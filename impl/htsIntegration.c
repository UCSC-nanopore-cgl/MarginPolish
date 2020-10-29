//
// Created by tpesout on 1/8/19.
//

#include <htslib/thread_pool.h>
#include <htslib/faidx.h>
#include "htsIntegration.h"
#include "margin.h"
#include "lp_lib.h"

//#include <htslib/vcf.h>
//#include <htslib/faidx.h>
//#include <htslib/bgzf.h>
#include "bedidx.h"

// This structure holds the bed information
// TODO rewrite the code to just use a void*
typedef struct samview_settings {
    void *bed;
} samview_settings_t;


/*
 * getAlignedReadLength computes the length of the read sequence which is aligned to the reference.  Hard-clipped bases
 * are never included in this calculation.  Soft-clipped bases are similarly not included, but will be returned via
 * the two parameters if the 2nd or 3rd version is invoked.  The boundaryAtMatch parameter handles cases where an
 * insertion or deletion is the first or last operation after or before clipping.  If set, these will be treated as
 * soft-clipping; otherwise they will included in the final return value.
 */
int64_t getAlignedReadLength(bam1_t *aln) {
    int64_t start_softclip = 0;
    int64_t end_softclip = 0;
    return getAlignedReadLength3(aln, &start_softclip, &end_softclip, TRUE);
}

int64_t getAlignedReadLength2(bam1_t *aln, int64_t *start_softclip, int64_t *end_softclip) {
    return getAlignedReadLength3(aln, start_softclip, end_softclip, TRUE);
}

int64_t getAlignedReadLength3(bam1_t *aln, int64_t *start_softclip, int64_t *end_softclip, bool boundaryAtMatch) {
    // start read needs to be init'd to 0 (mostly this is to avoid misuse)
    if (*start_softclip != 0) st_errAbort("getAlignedReadLength2 invoked with improper start_softclip parameter");
    if (*end_softclip != 0) st_errAbort("getAlignedReadLength2 invoked with improper end_softclip parameter");

    // get relevant cigar info
    int64_t len = aln->core.l_qseq;
    uint32_t *cigar = bam_get_cigar(aln);

    // data for tracking
    int64_t start_ref = 0;
    int64_t cig_idx = 0;

    // Find the correct starting locations on the read and reference sequence,
    // to deal with things like inserts / deletions / soft clipping
    while (cig_idx < aln->core.n_cigar) {
        int cigarOp = cigar[cig_idx] & BAM_CIGAR_MASK;
        int cigarNum = cigar[cig_idx] >> BAM_CIGAR_SHIFT;

        if (cigarOp == BAM_CMATCH || cigarOp == BAM_CEQUAL || cigarOp == BAM_CDIFF) {
            break;
        } else if (cigarOp == BAM_CDEL || cigarOp == BAM_CREF_SKIP) {
            if (!boundaryAtMatch) break;
            cig_idx++;
        } else if (cigarOp == BAM_CINS) {
            if (!boundaryAtMatch) break;
            *start_softclip += cigarNum;
            cig_idx++;
        } else if (cigarOp == BAM_CSOFT_CLIP) {
            *start_softclip += cigarNum;
            cig_idx++;
        } else if (cigarOp == BAM_CHARD_CLIP || cigarOp == BAM_CPAD) {
            cig_idx++;
        } else {
            st_errAbort("Unidentifiable cigar operation\n");
        }
    }

    // Check for soft clipping at the end
    cig_idx = aln->core.n_cigar - 1;
    while (cig_idx > 0) {
        int cigarOp = cigar[cig_idx] & BAM_CIGAR_MASK;
        int cigarNum = cigar[cig_idx] >> BAM_CIGAR_SHIFT;

        if (cigarOp == BAM_CMATCH || cigarOp == BAM_CEQUAL || cigarOp == BAM_CDIFF) {
            break;
        } else if (cigarOp == BAM_CDEL || cigarOp == BAM_CREF_SKIP) {
            if (!boundaryAtMatch) break;
            cig_idx--;
        } else if (cigarOp == BAM_CINS) {
            if (!boundaryAtMatch) break;
            *end_softclip += cigarNum;
            cig_idx--;
        } else if (cigarOp == BAM_CSOFT_CLIP) {
            *end_softclip += cigarNum;
            cig_idx--;
        } else if (cigarOp == BAM_CHARD_CLIP || cigarOp == BAM_CPAD) {
            cig_idx--;
        } else {
            st_errAbort("Unidentifiable cigar operation\n");
        }
    }

    // Count number of insertions & deletions in sequence
    int64_t numInsertions = 0;
    int64_t numDeletions = 0;
    countIndels(cigar, aln->core.n_cigar, &numInsertions, &numDeletions);
    int64_t trueLength = len - *start_softclip - *end_softclip + numDeletions - numInsertions;

    return trueLength;
}


/*
 * Counts the number of insertions and deletions in a read, given its cigar string.
 */
void countIndels(uint32_t *cigar, uint32_t ncigar, int64_t *numInsertions, int64_t *numDeletions) {
    for (uint32_t i = 0; i < ncigar; i++) {
        int cigarOp = cigar[i] & BAM_CIGAR_MASK;
        int cigarNum = cigar[i] >> BAM_CIGAR_SHIFT;
        if (cigarOp == BAM_CINS) *numInsertions += cigarNum;
        if (cigarOp == BAM_CDEL) *numDeletions += cigarNum;
    }
}

/*
 * Utility functions for BamChunk constructor
 */


int64_t getReadDepthInfoBucketSize(int64_t chunkSize) {
    int64_t depth = chunkSize / 32;
    if (depth == 0) depth = 1;
    return depth;
}

int64_t getEstimatedChunkDepth(stList *chunkDepths, int64_t contigStartPos, int64_t contigEndPos, int64_t chunkSize) {
    if (chunkDepths == NULL) return 0;
    int64_t bucketSize = getReadDepthInfoBucketSize(chunkSize);
    contigStartPos = contigStartPos/bucketSize;
    contigEndPos = contigEndPos/bucketSize;
    assert(stList_length(chunkDepths) > contigEndPos);
    int64_t totalSize = 0;
    for (int64_t pos = contigStartPos; pos < contigEndPos; pos++) {
        totalSize += (int64_t) stList_get(chunkDepths, pos);
    }
    int64_t chunkLength = contigEndPos - contigStartPos;
    if (chunkLength <= 0) chunkLength = 1;
    int64_t estimatedDepth = totalSize / chunkLength;
    return estimatedDepth;
}

int64_t saveContigChunks(stList *dest, BamChunker *parent, char *contig, int64_t contigStartPos, int64_t contigEndPos,
                         uint64_t chunkSize, uint64_t chunkMargin, stList *chunkDepths) {

    // whole contig case
    if (chunkSize == 0) {
        BamChunk *chunk = bamChunk_construct2(contig, stList_length(dest), contigStartPos, contigStartPos, contigEndPos,
                contigEndPos, getEstimatedChunkDepth(chunkDepths, contigStartPos, contigEndPos, chunkSize), parent);
        stList_append(dest, chunk);
        return 1;
    }

    // specific chunk size
    int64_t chunkCount = 0;
    for (int64_t i = contigStartPos; i < contigEndPos; i += chunkSize) {
        int64_t chunkEndPos = i + chunkSize;
        chunkEndPos = (chunkEndPos > contigEndPos ? contigEndPos : chunkEndPos);
        int64_t chunkMarginStartPos = i - chunkMargin;
        chunkMarginStartPos = (chunkMarginStartPos < contigStartPos ? contigStartPos : chunkMarginStartPos);
        int64_t chunkMarginEndPos = chunkEndPos + chunkMargin;
        chunkMarginEndPos = (chunkMarginEndPos > contigEndPos ? contigEndPos : chunkMarginEndPos);

        BamChunk *chunk = bamChunk_construct2(contig, stList_length(dest), chunkMarginStartPos, i, chunkEndPos,
                                              chunkMarginEndPos, getEstimatedChunkDepth(chunkDepths,
                                                      chunkMarginStartPos, chunkMarginEndPos, chunkSize), parent);
        stList_append(dest, chunk);
        chunkCount++;
    }
    return chunkCount;
}

void storeReadDepthInformation(stList *depthList, int64_t startPos, int64_t endPos, int64_t chunkSize) {
    int64_t bucketSize = getReadDepthInfoBucketSize(chunkSize);
    startPos = startPos/bucketSize;
    endPos = endPos/bucketSize;
    while (stList_length(depthList) <= endPos) {
        stList_append(depthList, (void*) 0);
    }
    for (int64_t pos = startPos; pos < endPos; pos++) {
        stList_set(depthList, pos, (void*) (((int64_t) stList_get(depthList, pos)) + 1));
    }
}


/*
 * These handle construction of the BamChunk object, by iterating through the bam (must be sorted), and finds the
 * first and last aligned location on each contig.  Then it generates a list of chunks based off of these positions,
 * with sizes determined by the parameters.
 */
BamChunker *bamChunker_construct(char *bamFile, PolishParams *params) {
    return bamChunker_construct2(bamFile, NULL, params, false);
}

BamChunker *bamChunker_construct2(char *bamFile, char *regionStr, PolishParams *params, bool recordFilteredReads) {

    // are we doing region filtering?
    bool filterByRegion = false;
    char regionContig[128] = "";
    int regionStart = 0;
    int regionEnd = 0;
    if (regionStr != NULL) {
        int scanRet = sscanf(regionStr, "%[^:]:%d-%d", regionContig, &regionStart, &regionEnd);
        if (scanRet != 3 && scanRet != 1) {
            st_errAbort("Region in unexpected format (expected %%s:%%d-%%d or %%s)): %s", regionStr);
        } else if (regionStart < 0 || regionEnd < 0 || regionEnd < regionStart) {
            st_errAbort("Start and end locations in region must be positive, start must be less than end: %s", regionStr);
        }
        filterByRegion = true;
    }

    // standard parameters
    uint64_t chunkSize = params->chunkSize;
    uint64_t chunkBoundary = params->chunkBoundary;
    bool includeSoftClip = params->includeSoftClipping;

    // the chunker we're building
    BamChunker *chunker = malloc(sizeof(BamChunker));
    chunker->bamFile = stString_copy(bamFile);
    chunker->chunkSize = chunkSize;
    chunker->chunkBoundary = chunkBoundary;
    chunker->includeSoftClip = includeSoftClip;
    chunker->params = params;
    chunker->chunks = stList_construct3(0, (void *) bamChunk_destruct);
    chunker->chunkCount = 0;

    // file initialization
    samFile *in = NULL;
    hts_idx_t *idx = NULL;
    hts_itr_multi_t *iter = NULL;
    // bam file
    if ((in = hts_open(bamFile, "r")) == 0) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamFile);
    }
    // bam index
    if ((idx = sam_index_load(in, bamFile)) == 0) {
        st_errAbort("ERROR: Cannot open index for bam file %s\n", bamFile);
    }
    // header
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    // read object
    bam1_t *aln = bam_init1();

    // thread stuff
    htsThreadPool threadPool = {NULL, 0};
    # ifdef _OPENMP
    int tc = omp_get_max_threads();
    # else
    int tc = 1;
    # endif
    if (!(threadPool.pool = hts_tpool_init(tc))) {
        fprintf(stderr, "Error creating thread pool\n");
    }
    hts_set_opt(in, HTS_OPT_THREAD_POOL, &threadPool);

    // prep for index (not entirely sure what all this does.  see samtools/sam_view.c
    int filter_state = ALL, filter_op = 0;
    samview_settings_t settings = {.bed = NULL};
    char *region[1] = {};
    int regcount = 0;
    hts_reglist_t *reglist = NULL;
    if (regionStr != NULL) {
        region[0] = stString_copy(regionStr);
        settings.bed = bed_hash_regions(settings.bed, region, 0, 1,
                                        &filter_op); //insert(1) or filter out(0) the regions from the command line in the same hash table as the bed file
        if (!filter_op) filter_state = FILTERED;
        reglist = bed_reglist(settings.bed, filter_state, &regcount);
        if (!reglist) {
            st_errAbort("ERROR: Could not create list of regions for read conversion");
        }
        if ((iter = sam_itr_regions(idx, bamHdr, reglist, regcount)) == 0) {
            st_errAbort("ERROR: Cannot open iterator for region %s for bam file %s\n", region[0], bamFile);
        }
    }

    // list of chunk boundaries
    char *currentContig = NULL;
    int64_t contigStartPos = 0;
    int64_t contigEndPos = 0;
    stList *estimatedChunkDepths = stList_construct();

    // get all reads
    // there is probably a better way (bai?) to find min and max aligned positions (which we need for chunk divisions)
    while ((regionStr == NULL ? sam_read1(in, bamHdr, aln) : sam_itr_multi_next(in, iter, aln)) > 0) {

        // basic filtering (no read length, no cigar)
        if (aln->core.l_qseq <= 0) continue;
        if (aln->core.n_cigar == 0) continue;
        if ((aln->core.flag & (uint16_t) 0x4) != 0)
            continue; //unaligned
        if (!params->includeSecondaryAlignments && (aln->core.flag & (uint16_t) 0x100) != 0)
            continue; //secondary
        if (!params->includeSupplementaryAlignments && (aln->core.flag & (uint16_t) 0x800) != 0)
            continue; //supplementary
        if (aln->core.qual < params->filterAlignmentsWithMapQBelowThisThreshold) { //low mapping quality
            if (!recordFilteredReads) continue;
        }

        //data
        char *chr = bamHdr->target_name[aln->core.tid];
        int64_t start_softclip = 0;
        int64_t end_softclip = 0;
        int64_t alnReadLength = getAlignedReadLength3(aln, &start_softclip, &end_softclip, FALSE);
        int64_t alnStartPos = aln->core.pos;
        int64_t alnEndPos = alnStartPos + alnReadLength;

        // does this belong in our chunk?
        if (alnReadLength <= 0) continue;

        // get start and stop position
        int64_t readStartPos = aln->core.pos;           // Left most position of alignment
        int64_t readEndPos = readStartPos + alnReadLength;

        // get contig
        char *contig = bamHdr->target_name[aln->core.tid];     // Contig name

        if (currentContig == NULL) {
            // first read
            currentContig = stString_copy(contig);
            contigStartPos = readStartPos;
            contigEndPos = readEndPos;
            storeReadDepthInformation(estimatedChunkDepths, readStartPos, readEndPos, chunkSize);
        } else if (stString_eq(currentContig, contig)) {
            // continue this contig's reads
            contigStartPos = readStartPos < contigStartPos ? readStartPos : contigStartPos;
            contigEndPos = readEndPos > contigEndPos ? readEndPos : contigEndPos;
            storeReadDepthInformation(estimatedChunkDepths, readStartPos, readEndPos, chunkSize);
        } else {
            // new contig (this should never happen if we're filtering by region)
            assert(!filterByRegion);
            int64_t savedChunkCount = saveContigChunks(chunker->chunks, chunker, currentContig,
                                                       contigStartPos, contigEndPos, chunkSize, chunkBoundary,
                                                       estimatedChunkDepths);
            chunker->chunkCount += savedChunkCount;
            free(currentContig);
            currentContig = stString_copy(contig);
            contigStartPos = readStartPos;
            contigEndPos = readEndPos;
            stList_destruct(estimatedChunkDepths);
            estimatedChunkDepths = stList_construct();
        }
    }
    // save last contig's chunks
    if (currentContig != NULL) {
        if (filterByRegion && regionStart != 0 && regionEnd != 0) {
            contigStartPos = (contigStartPos < regionStart ? regionStart : contigStartPos);
            contigEndPos = (contigEndPos > regionEnd ? regionEnd : contigEndPos);
        }
        int64_t savedChunkCount = saveContigChunks(chunker->chunks, chunker, currentContig,
                                                   contigStartPos, contigEndPos, chunkSize, chunkBoundary,
                                                   estimatedChunkDepths);
        chunker->chunkCount += savedChunkCount;
        free(currentContig);
    }

    // sanity check
    assert(stList_length(chunker->chunks) == chunker->chunkCount);

    // shut everything down
    if (regionStr != NULL) {
        hts_itr_multi_destroy(iter);
        free(region[0]);
        bed_destroy(settings.bed);
    }
    hts_idx_destroy(idx);
    stList_destruct(estimatedChunkDepths);
    bam_hdr_destroy(bamHdr);
    bam_destroy1(aln);
    sam_close(in);
    hts_tpool_destroy(threadPool.pool);

    return chunker;
}

BamChunker *bamChunker_constructFromFasta(char *fastaFile, char *bamFile, char *regionStr, PolishParams *params) {
    faidx_t *fai = fai_load_format(fastaFile, FAI_FASTA);
    if ( !fai ) {
        st_errAbort("[faidx] Could not load fai index of %s\n", fastaFile);
    }

    // standard parameters
    uint64_t chunkSize = params->chunkSize;
    uint64_t chunkBoundary = params->chunkBoundary;
    bool includeSoftClip = params->includeSoftClipping;

    // the chunker we're building
    BamChunker *chunker = malloc(sizeof(BamChunker));
    chunker->bamFile = stString_copy(bamFile);
    chunker->chunkSize = chunkSize;
    chunker->chunkBoundary = chunkBoundary;
    chunker->includeSoftClip = includeSoftClip;
    chunker->params = params;
    chunker->chunks = stList_construct3(0, (void *) bamChunk_destruct);
    chunker->chunkCount = 0;

    if (regionStr != NULL) {
        // prep parsing seq
        char regionContig[128] = "";
        int regionStart = 0;
        int regionEnd = 0;
        int scanRet = sscanf(regionStr, "%[^:]:%d-%d", regionContig, &regionStart, &regionEnd);
        if (scanRet != 3 && scanRet != 1) {
            st_errAbort("Region in unexpected format (expected %%s:%%d-%%d or %%s)): %s", regionStr);
        } else if (regionStart < 0 || regionEnd < 0 || regionEnd < regionStart) {
            st_errAbort("Start and end locations in region must be positive, start must be less than end: %s", regionStr);
        }
        // get actual values
        int seqLen;
        char *seq = fai_fetch(fai, regionStr, &seqLen);
        int64_t savedChunkCount = saveContigChunks(chunker->chunks, chunker, regionContig,
                                                   regionStart, regionStart + seqLen - 1, chunkSize, chunkBoundary,
                                                   NULL);
        chunker->chunkCount += savedChunkCount;
        free(seq);
    } else {
        char *faiFile = stString_print("%s.fai", fastaFile);
        FILE *fp = fopen(faiFile, "r");
        char *line = NULL;
        while ((line = stFile_getLineFromFile(fp)) != NULL) {
            stList *parts = stString_splitByString(line, "\t");
            if (stList_length(parts) < 2) st_errAbort("Unexpected fasta index file format: %s\n", faiFile);

            int64_t savedChunkCount = saveContigChunks(chunker->chunks, chunker, stList_get(parts, 0),
                                                       0, atol(stList_get(parts, 1)), chunkSize, chunkBoundary,
                                                       NULL);
            chunker->chunkCount += savedChunkCount;
        }
    }
    assert(chunker->chunkCount == stList_length(chunker->chunks));

    // cleanup and close
    fai_destroy(fai);
    return chunker;
}

BamChunker *bamChunker_copyConstruct(BamChunker *toCopy) {
    BamChunker *chunker = malloc(sizeof(BamChunker));
    chunker->bamFile = stString_copy(toCopy->bamFile);
    chunker->chunkSize = toCopy->chunkSize;
    chunker->chunkBoundary = toCopy->chunkBoundary;
    chunker->includeSoftClip = toCopy->includeSoftClip;
    chunker->params = toCopy->params;
    chunker->chunks = stList_construct3(0, (void *) bamChunk_destruct);
    chunker->chunkCount = 0;
    return chunker;
}

void bamChunker_destruct(BamChunker *bamChunker) {
    if (bamChunker->bamFile != NULL) free(bamChunker->bamFile);
    stList_destruct(bamChunker->chunks);
    free(bamChunker);
}

BamChunk *bamChunker_getChunk(BamChunker *bamChunker, int64_t chunkIdx) {
    BamChunk *chunk = stList_get(bamChunker->chunks, chunkIdx);
    return chunk;
}

BamChunk *bamChunk_construct() {
    return bamChunk_construct2(NULL, 0, 0, 0, 0, 0, 0, NULL);
}

BamChunk *bamChunk_construct2(char *refSeqName, int64_t chunkIndex, int64_t chunkOverlapStart, int64_t chunkStart, int64_t chunkEnd,
                              int64_t chunkOverlapEnd, int64_t depth, BamChunker *parent) {
    BamChunk *c = malloc(sizeof(BamChunk));
    c->chunkIdx = chunkIndex;
    c->refSeqName = stString_copy(refSeqName);
    c->chunkOverlapStart = chunkOverlapStart;
    c->chunkStart = chunkStart;
    c->chunkEnd = chunkEnd;
    c->chunkOverlapEnd = chunkOverlapEnd;
    c->estimatedDepth = depth;
    c->parent = parent;
    return c;
}

BamChunk *bamChunk_copyConstruct(BamChunk *toCopy) {
    BamChunk *c = malloc(sizeof(BamChunk));
    c->chunkIdx = toCopy->chunkIdx;
    c->refSeqName = stString_copy(toCopy->refSeqName);
    c->chunkOverlapStart = toCopy->chunkOverlapStart;
    c->chunkStart = toCopy->chunkStart;
    c->chunkEnd = toCopy->chunkEnd;
    c->chunkOverlapEnd = toCopy->chunkOverlapEnd;
    c->estimatedDepth = toCopy->estimatedDepth;
    c->parent = toCopy->parent;
    return c;
}

void bamChunk_destruct(BamChunk *bamChunk) {
    free(bamChunk->refSeqName);
    free(bamChunk);
}

int compareBamChunkDepthByIndexInList(const void *a, const void *b, const void *chunkList) {
    const BamChunk *chunk1 = stList_get((stList*)chunkList, stIntTuple_get((stIntTuple*) a, 0));
    const BamChunk *chunk2 = stList_get((stList*)chunkList, stIntTuple_get((stIntTuple*) b, 0));
    return chunk1->estimatedDepth < chunk2->estimatedDepth ? -1 : chunk1->estimatedDepth > chunk2->estimatedDepth ? 1 : 0;
}

/*
 * This generates a set of BamChunkReads (and alignments to the reference) from a BamChunk.  The BamChunk describes
 * positional information within the bam, from which the reads should be extracted.  The bam must be indexed.  Reads
 * which overlap the ends of the chunk are truncated.  A parameter in the BamChunk's parameters determines whether
 * softclipped portions of the reads should be included.
 */

uint32_t convertToReadsAndAlignmentsWithFiltered(BamChunk *bamChunk, RleString *reference, stList *reads,
        stList *alignments, stList *filteredReads, stList *filteredAlignments, PolishParams *polishParams) {

    // sanity check
    assert(stList_length(reads) == 0);
    assert(stList_length(alignments) == 0);

    uint64_t *ref_nonRleToRleCoordinateMap =
            reference == NULL ? NULL : rleString_getNonRleToRleCoordinateMap(reference);

    // prep
    int64_t chunkStart = bamChunk->chunkOverlapStart;
    int64_t chunkEnd = bamChunk->chunkOverlapEnd;
    bool includeSoftClip = bamChunk->parent->params->includeSoftClipping;
    char *bamFile = bamChunk->parent->bamFile;
    char *contig = bamChunk->refSeqName;
    uint32_t savedAlignments = 0;
    double randomDiscardChance = 1.0;
    //removing reads in filtering step is working, so taking this out for now
    //  (also removing the check in the loop below)
    /*if (bamChunk->estimatedDepth > polishParams->excessiveDepthThreshold) {
        char *logIdentifier = getLogIdentifier();
        randomDiscardChance = 1.0 * polishParams->excessiveDepthThreshold / bamChunk->estimatedDepth;
        st_logInfo(" %s Randomly removing reads from excessively deep (%"PRId64"/%"PRId64") chunk with chance %f\n",
                logIdentifier, bamChunk->estimatedDepth, polishParams->excessiveDepthThreshold, 1.0 - randomDiscardChance);
    }*/

    // prep for index (not entirely sure what all this does.  see samtools/sam_view.c
    int filter_state = ALL, filter_op = 0;
    int result;
    samview_settings_t settings = {.bed = NULL};
    char *region[1] = {};
    region[0] = stString_print("%s:%d-%d", bamChunk->refSeqName, bamChunk->chunkOverlapStart,
                               bamChunk->chunkOverlapEnd);
    settings.bed = bed_hash_regions(settings.bed, region, 0, 1,
                                    &filter_op); //insert(1) or filter out(0) the regions from the command line in the same hash table as the bed file
    if (!filter_op) filter_state = FILTERED;
    int regcount = 0;
    hts_reglist_t *reglist = bed_reglist(settings.bed, filter_state, &regcount);
    if (!reglist) {
        st_errAbort("ERROR: Could not create list of regions for read conversion");
    }

    // file initialization
    samFile *in = NULL;
    hts_idx_t *idx = NULL;
    hts_itr_multi_t *iter = NULL;
    // bam file
    if ((in = hts_open(bamFile, "r")) == 0) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamFile);
    }
    // bam index
    if ((idx = sam_index_load(in, bamFile)) == 0) {
        st_errAbort("ERROR: Cannot open index for bam file %s\n", bamFile);
    }
    // header  //todo samFile *in = hts_open(bamFile, "r");
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    // read object
    bam1_t *aln = bam_init1();
    // iterator for region
    if ((iter = sam_itr_regions(idx, bamHdr, reglist, regcount)) == 0) {
        st_errAbort("ERROR: Cannot open iterator for region %s for bam file %s\n", region[0], bamFile);
    }

    // fetch alignments //todo while(sam_read1(in,bamHdr,aln) > 0) {
    while ((result = sam_itr_multi_next(in, iter, aln)) >= 0) {
        bool filtered = FALSE;
        // basic filtering (no read length, no cigar)
        if (aln->core.l_qseq <= 0) continue;
        if (aln->core.n_cigar == 0) continue;
        if ((aln->core.flag & (uint16_t) 0x4) != 0)
            continue; //unaligned
        if (!polishParams->includeSecondaryAlignments && (aln->core.flag & (uint16_t) 0x100) != 0)
            continue; //secondary
        if (!polishParams->includeSupplementaryAlignments && (aln->core.flag & (uint16_t) 0x800) != 0)
            continue; //supplementary
        // see above, taking this out as removal in filtering step is working
        /*if (st_random() > randomDiscardChance)
            continue; // chunk is too deep*/
        if (aln->core.qual < polishParams->filterAlignmentsWithMapQBelowThisThreshold) { //low mapping quality
            if (filteredReads == NULL) continue;
            filtered = TRUE;
        }

        //data
        char *chr = bamHdr->target_name[aln->core.tid];
        int64_t start_softclip = 0;
        int64_t end_softclip = 0;
        int64_t alnReadLength = getAlignedReadLength3(aln, &start_softclip, &end_softclip, FALSE);
        if (alnReadLength <= 0) continue;
        int64_t alnStartPos = aln->core.pos;
        int64_t alnEndPos = alnStartPos + alnReadLength;

        // does this belong in our chunk?
        if (!stString_eq(contig, chr)) continue;
        if (alnStartPos >= chunkEnd) continue;
        if (alnEndPos <= chunkStart) continue;

        // get cigar and rep
        uint32_t *cigar = bam_get_cigar(aln);
        stList *cigRepr = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

        // Variables to keep track of position in sequence / cigar operations
        int64_t cig_idx = 0;
        int64_t currPosInOp = 0;
        int64_t cigarOp = -1;
        int64_t cigarNum = -1;
        int64_t cigarIdxInSeq = 0;
        int64_t cigarIdxInRef = alnStartPos;

        // positional modifications
        int64_t refCigarModification = -1 * chunkStart;

        // we need to calculate:
        //  a. where in the (potentially softclipped read) to start storing characters
        //  b. what the alignments are wrt those characters
        // so we track the first aligned character in the read (for a.) and what alignment modification to make (for b.)
        int64_t seqCigarModification;
        int64_t firstNonSoftclipAlignedReadIdxInChunk;

        // the handling changes based on softclip inclusion and where the chunk boundaries are
        if (includeSoftClip) {
            if (alnStartPos < chunkStart) {
                // alignment spans chunkStart (this will not be affected by softclipping)
                firstNonSoftclipAlignedReadIdxInChunk = -1; //need to find position of first alignment
                seqCigarModification = 0;
            } else if (alnStartPos - start_softclip <= chunkStart) {
                // softclipped bases span chunkStart
                firstNonSoftclipAlignedReadIdxInChunk = 0;
                int64_t includedSoftclippedBases = alnStartPos - chunkStart;
                seqCigarModification = includedSoftclippedBases;
                assert(includedSoftclippedBases >= 0);
                assert(start_softclip - includedSoftclippedBases >= 0);
            } else {
                // softclipped bases are after chunkStart
                firstNonSoftclipAlignedReadIdxInChunk = 0;
                seqCigarModification = start_softclip;
            }
        } else {
            if (alnStartPos < chunkStart) {
                // alignment spans chunkStart
                firstNonSoftclipAlignedReadIdxInChunk = -1;
                seqCigarModification = 0;
            } else {
                // alignment starts after chunkStart
                firstNonSoftclipAlignedReadIdxInChunk = 0;
                seqCigarModification = 0;
            }
        }

        // track number of characters in aligned portion (will inform softclipping at end of read)
        int64_t alignedReadLength = 0;

        // iterate over cigar operations
        for (uint32_t i = 0; i <= alnReadLength; i++) {
            // handles cases where last alignment is an insert or last is match
            if (cig_idx == aln->core.n_cigar) break;

            // do we need the next cigar operation?
            if (currPosInOp == 0) {
                cigarOp = cigar[cig_idx] & BAM_CIGAR_MASK;
                cigarNum = cigar[cig_idx] >> BAM_CIGAR_SHIFT;
            }

            // handle current character
            if (cigarOp == BAM_CMATCH || cigarOp == BAM_CEQUAL || cigarOp == BAM_CDIFF) {
                if (cigarIdxInRef >= chunkStart && cigarIdxInRef < chunkEnd) {
                    stList_append(cigRepr, stIntTuple_construct3(cigarIdxInRef + refCigarModification,
                                                                 cigarIdxInSeq + seqCigarModification,
                                                                 polishParams->p->diagonalExpansion));
//                                                                 polishParams != NULL && polishParams->p != NULL ?
//                                                                 polishParams->p->diagonalExpansion : 10)); //TODO: Tidy up so polish params is not optional
                    alignedReadLength++;
                }
                cigarIdxInSeq++;
                cigarIdxInRef++;
            } else if (cigarOp == BAM_CDEL || cigarOp == BAM_CREF_SKIP) {
                //delete
                cigarIdxInRef++;
            } else if (cigarOp == BAM_CINS) {
                //insert
                cigarIdxInSeq++;
                if (cigarIdxInRef >= chunkStart && cigarIdxInRef < chunkEnd) {
                    alignedReadLength++;
                }
                i--;
            } else if (cigarOp == BAM_CSOFT_CLIP || cigarOp == BAM_CHARD_CLIP || cigarOp == BAM_CPAD) {
                // nothing to do here. skip to next cigar operation
                currPosInOp = cigarNum - 1;
                i--;
            } else {
                st_logCritical("Unidentifiable cigar operation!\n");
            }

            // document read index in the chunk (for reads that span chunk boundary, used in read construction)
            if (firstNonSoftclipAlignedReadIdxInChunk < 0 && cigarIdxInRef >= chunkStart) {
                firstNonSoftclipAlignedReadIdxInChunk = cigarIdxInSeq;
                seqCigarModification = -1 * (firstNonSoftclipAlignedReadIdxInChunk + seqCigarModification);

            }

            // have we finished this last cigar
            currPosInOp++;
            if (currPosInOp == cigarNum) {
                cig_idx++;
                currPosInOp = 0;
            }
        }

        // get sequence positions
        int64_t seqLen = alignedReadLength;

        // modify start indices
        int64_t readStartIdxInChunk = firstNonSoftclipAlignedReadIdxInChunk;
        if (firstNonSoftclipAlignedReadIdxInChunk != 0) {
            // the aligned portion spans chunkStart, so no softclipped bases are included
            readStartIdxInChunk += start_softclip;
        } else if (!includeSoftClip) {
            // configured to not handle softclipped bases
            readStartIdxInChunk += start_softclip;
        } else if (alnStartPos - start_softclip <= chunkStart) {
            // configured to handle softclipped bases; softclipped bases span chunkStart
            int64_t includedSoftclippedBases = alnStartPos - chunkStart;
            seqLen += includedSoftclippedBases;
            readStartIdxInChunk += (start_softclip - includedSoftclippedBases);
        } else {
            // configured to handle softclipped bases; softclipped bases all occur after chunkStart
            seqLen += start_softclip;
            readStartIdxInChunk = 0;
        }

        // modify end indices
        int64_t readEndIdxInChunk = readStartIdxInChunk + seqLen;
        if (alnEndPos < chunkEnd && includeSoftClip) {
            // all other cases mean we don't need to handle softclip (by config or aln extends past chunk end)
            if (alnEndPos + end_softclip <= chunkEnd) {
                // all softclipped bases fit in chunk
                readEndIdxInChunk += end_softclip;
                seqLen += end_softclip;
            } else {
                // softclipping spands chunkEnd
                int64_t includedSoftclippedBases = chunkEnd - alnEndPos;
                seqLen += includedSoftclippedBases;
                readEndIdxInChunk += includedSoftclippedBases;
            }
        }

        // get sequence - all data we need is encoded in readStartIdxInChunk (start), readEnd idx, and seqLen
        char *seq = st_calloc(seqLen + 1, sizeof(char));
        uint8_t *seqBits = bam_get_seq(aln);
        int64_t idxInOutputSeq = 0;
        int64_t idxInBamRead = readStartIdxInChunk;
        while (idxInBamRead < readEndIdxInChunk) {
            seq[idxInOutputSeq] = seq_nt16_str[bam_seqi(seqBits, idxInBamRead)];
            idxInBamRead++;
            idxInOutputSeq++;
        }
        seq[seqLen] = '\0';

        // get sequence qualities (if exists)
        char *readName = stString_copy(bam_get_qname(aln));
        uint8_t *qualBits = bam_get_qual(aln);
        uint8_t *qual = NULL;
        if (qualBits[0] != 0xff) { //inital score of 255 means qual scores are unavailable
            idxInOutputSeq = 0;
            idxInBamRead = readStartIdxInChunk;
            qual = st_calloc(seqLen, sizeof(uint8_t));
            while (idxInBamRead < readEndIdxInChunk) {
                qual[idxInOutputSeq] = qualBits[idxInBamRead];
                idxInBamRead++;
                idxInOutputSeq++;

            }
            assert(idxInOutputSeq == strlen(seq));
        };

        // failure case
        if (stList_length(cigRepr) == 0 || strlen(seq) == 0) {
            stList_destruct(cigRepr);
            free(readName);
            free(seq);
            if (qual != NULL) free(qual);
            continue;
        }

        // sanity check
        assert(stIntTuple_get((stIntTuple *) stList_peek(cigRepr), 1) < strlen(seq));

        // save to read
        bool forwardStrand = !bam_is_rev(aln);
        BamChunkRead *chunkRead = bamChunkRead_construct3(readName, seq, qual, forwardStrand, aln->l_data,
                                                          polishParams->useRunLengthEncoding);
        stList_append(filtered ? filteredReads: reads, chunkRead);

        // save alignment
        if (polishParams->useRunLengthEncoding) {
            // ref_nonRleToRleCoordinateMap should only be null w/ RLE in tests
            if (ref_nonRleToRleCoordinateMap != NULL) {
                // rle the alignment and save it
                uint64_t *read_nonRleToRleCoordinateMap = rleString_getNonRleToRleCoordinateMap(chunkRead->rleRead);
                stList_append(filtered ? filteredAlignments : alignments,
                        runLengthEncodeAlignment(cigRepr, ref_nonRleToRleCoordinateMap, read_nonRleToRleCoordinateMap));
                stList_destruct(cigRepr);
                free(read_nonRleToRleCoordinateMap);
            }
        } else {
            stList_append(filtered ? filteredAlignments : alignments, cigRepr);
        }
        savedAlignments++;

        // cleanup
        free(readName);
        free(seq);
        if (qual != NULL) free(qual);
    }
    // the status from "get reads from iterator"
    if (result < -1) {
        st_errAbort("ERROR: Retrieval of region %d failed due to truncated file or corrupt BAM index file\n",
                    iter->curr_tid);
    }

    // close it all down


    hts_itr_multi_destroy(iter);
    hts_idx_destroy(idx);
    free(region[0]);
    bed_destroy(settings.bed);
    bam_hdr_destroy(bamHdr);
    bam_destroy1(aln);
    sam_close(in);
    if (ref_nonRleToRleCoordinateMap != NULL)
        free(ref_nonRleToRleCoordinateMap);
    return savedAlignments;
}

uint32_t convertToReadsAndAlignments(BamChunk *bamChunk, RleString *reference, stList *reads, stList *alignments,
                                     PolishParams *polishParams) {
    return convertToReadsAndAlignmentsWithFiltered(bamChunk, reference, reads, alignments, NULL, NULL, polishParams);
}

bool downsampleViaReadLikelihood(int64_t intendedDepth, BamChunk *bamChunk, stList *inputReads, stList *inputAlignments,
                                 stList *maintainedReads, stList *maintainedAlignments, stList *discardedReads,
                                 stList *discardedAlignments) {

    // calculate depth
    int64_t totalNucleotides = 0;
    for (int64_t i = 0; i < stList_length(inputReads); i++) {
        BamChunkRead *bcr = stList_get(inputReads, i);
        totalNucleotides += bcr->rleRead->length;
    }
    int64_t chunkSize = bamChunk->chunkOverlapEnd - bamChunk->chunkOverlapStart;
    double averageDepth = 1.0 * totalNucleotides / (chunkSize);

    // do we need to downsample?
    if (averageDepth < intendedDepth) {
        return FALSE;
    }

    // we do need to downsample
    char *logIdentifier = getLogIdentifier();
    st_logInfo(" %s Downsampling chunk with average depth %.2fx to %dx \n", logIdentifier, averageDepth, intendedDepth);
    free(logIdentifier);

    // keep some ratio of reads
    double ratioToKeep = intendedDepth / averageDepth;
    for (int64_t i = 0; i < stList_length(inputReads); i++) {
        if (st_random() < ratioToKeep) {
            stList_append(maintainedReads, stList_get(inputReads, i));
            stList_append(maintainedAlignments, stList_get(inputAlignments, i));
        } else {
            stList_append(discardedReads, stList_get(inputReads, i));
            stList_append(discardedAlignments, stList_get(inputAlignments, i));
        }
    }

    return TRUE;
}

#define HET_SITE_SCALE 2
int64_t countVcfSitesInAlignmentRange(stList *vcfEntries, stList *alignment) {
    if (stList_length(alignment) <= 1) return 0;
    int64_t startPos = stIntTuple_get(stList_get(alignment, 0), 0);
    int64_t endPos = stIntTuple_get(stList_get(alignment, stList_length(alignment) - 1), 0);
    int64_t sites = 0;
    for (int64_t i = 0; i < stList_length(vcfEntries); i++) {
        VcfEntry *vcfEntry = stList_get(vcfEntries, i);
        if (vcfEntry->refPos < startPos) continue;
        if (vcfEntry->refPos >= endPos) break;
        sites++;
    }
    return sites;
}

/* Copywrite Jordan "Big Brain" Eizinga
 * maximize    \sum_i p_i * h_i
 * subject to  \sum_i l_i * p_i = C * L
 *                          p_i ≥ 0       i = 1...N
 *                          p_i ≤ 1       i = 1...N
 */
double* computeReadProbsByLengthAndSecondMetric(int *read_lengths, int *read_metric, int num_reads,
                                                double target_coverage, int region_length) {
    lprec* lp = make_lp(0, num_reads);
    assert(lp != NULL);
    // this makes adding constraints more efficient
    set_add_rowmode(lp, true);
    // make the equality constraint to control expected coverage
    double total_cov = target_coverage * region_length;
    REAL* lens = (REAL*) malloc(sizeof(REAL) * (num_reads + 1));
    lens[0] = 0; // why do these libraries have this dumb 1-based shit?
    for (int i = 1; i <= num_reads; ++i) {
        lens[i] = (REAL) read_lengths[i - 1];
    }
    add_constraint(lp, lens, EQ, total_cov);
    free(lens);
    // constrain each probability to [0, 1] (these are sparse constraints)
    int* colno = (int*) malloc(sizeof(int));
    REAL* coef = (REAL*) malloc(sizeof(REAL));
    *coef = 1.0;
    for (int i = 1; i <= num_reads; ++i) {
        *colno = i;
        add_constraintex(lp, 1, coef, colno, LE, 1.0);
        add_constraintex(lp, 1, coef, colno, GE, 0.0);
    }
    free(colno);
    free(coef);
    // we're done adding constraints
    set_add_rowmode(lp, false);
    // set the objective function
    REAL* hets = (REAL*) malloc(sizeof(REAL) * (num_reads + 1));
    for (int i = 1; i <= num_reads; ++i) {
        hets[i] = (REAL) read_metric[i - 1];
    }
    set_obj_fn(lp, hets);
    free(hets);
    // we want to maximize the objective
    set_maxim(lp);
    // only print important logging info
    set_verbose(lp, IMPORTANT);
    // solve the problem
    int return_code = solve(lp);
    assert(return_code == OPTIMAL);
    // retrieve the results
    REAL* probs_real = (REAL*) malloc(sizeof(REAL) * num_reads);
    get_variables(lp, probs_real);
    // maybe convert to double
    double* probs = (double*) malloc(sizeof(double) * num_reads);
    for (int i = 0; i < num_reads; ++i) {
        probs[i] = (double) probs_real[i];
    }
    free(probs_real);
    // housekeeping
    delete_lp(lp);
    return probs;
}

bool downsampleViaHetSpanLikelihood(int64_t intendedDepth, BamChunk *bamChunk, stList *vcfEntries,
                                    stList *inputReads, stList *inputAlignments, stList *maintainedReads, stList *maintainedAlignments,
                                    stList *discardedReads, stList *discardedAlignments) {

    // calculate depth, generate bcrwhs list
    int64_t totalNucleotides = 0;
    int64_t totalHetSiteCount = 0;
    int64_t totalScaledHetSiteCount = 0;

    int *readLengths = st_calloc(stList_length(inputReads), sizeof(int));
    int *readHets = st_calloc(stList_length(inputReads), sizeof(int));

    for (int64_t i = 0; i < stList_length(inputReads); i++) {
        BamChunkRead *bcr = stList_get(inputReads, i);
        totalNucleotides += bcr->rleRead->length;
        readLengths[i] = (int) bcr->rleRead->length;
        int64_t hetCount = countVcfSitesInAlignmentRange(vcfEntries, stList_get(inputAlignments, i));
        totalHetSiteCount += hetCount;
        // give reads a pseudocount
        hetCount = HET_SITE_SCALE * hetCount + 1;
        readHets[i] = (int) hetCount;
        totalScaledHetSiteCount += hetCount;
    }
    int64_t chunkSize = bamChunk->chunkOverlapEnd - bamChunk->chunkOverlapStart;
    double averageDepth = 1.0 * totalNucleotides / (chunkSize);

    // do we need to downsample?
    if (averageDepth < intendedDepth) {
        free(readLengths);
        free(readHets);
        return FALSE;
    }

    // we do need to downsample
    char *logIdentifier = getLogIdentifier();
    st_logInfo(" %s Downsampling chunk via het sites with average depth %.2fx to %dx over %"PRId64" total spanned HET sites\n",
               logIdentifier, averageDepth, intendedDepth, totalHetSiteCount);

    // get likelihood of keeping each read
    double *probs = computeReadProbsByLengthAndSecondMetric(readLengths, readHets, stList_length(inputReads),
                                                            intendedDepth,
                                                            chunkSize);

    // keep some ratio of reads
    int64_t totalKeptNucleotides = 0;
    for (int64_t i = 0; i < stList_length(inputReads); i++) {
        BamChunkRead *bcr = stList_get(inputReads, i);
        if (st_random() < probs[i]) {
            stList_append(maintainedReads, stList_get(inputReads, i));
            stList_append(maintainedAlignments, stList_get(inputAlignments, i));
            totalKeptNucleotides += bcr->rleRead->length;
        } else {
            stList_append(discardedReads, stList_get(inputReads, i));
            stList_append(discardedAlignments, stList_get(inputAlignments, i));
        }
    }
    st_logInfo(" %s Downsampled chunk via het sites to average depth %.2fx (expected %dx)\n",
               logIdentifier, 1.0 * totalKeptNucleotides / chunkSize, intendedDepth);

    free(probs);
    free(readLengths);
    free(readHets);
    free(logIdentifier);
    return TRUE;
}

bool downsampleViaFullReadLengthLikelihood(int64_t intendedDepth, BamChunk *bamChunk, stList *inputReads,
        stList *inputAlignments, stList *maintainedReads, stList *maintainedAlignments,
        stList *discardedReads, stList *discardedAlignments) {

    // calculate depth, generate bcrwhs list
    int64_t totalNucleotides = 0;

    int *readLengths = st_calloc(stList_length(inputReads), sizeof(int));
    int *readFullLengths = st_calloc(stList_length(inputReads), sizeof(int));

    for (int64_t i = 0; i < stList_length(inputReads); i++) {
        BamChunkRead *bcr = stList_get(inputReads, i);
        totalNucleotides += bcr->rleRead->length;
        readLengths[i] = (int) bcr->rleRead->length;
        // give reads a pseudocount
        readFullLengths[i] = (int) bcr->fullReadLength;
    }
    int64_t chunkSize = bamChunk->chunkOverlapEnd - bamChunk->chunkOverlapStart;
    double averageDepth = 1.0 * totalNucleotides / (chunkSize);

    // do we need to downsample?
    if (averageDepth < intendedDepth) {
        free(readLengths);
        free(readFullLengths);
        return FALSE;
    }

    // we do need to downsample
    char *logIdentifier = getLogIdentifier();
    st_logInfo(" %s Downsampling chunk via full read length with average depth %.2fx to %dx.\n",
               logIdentifier, averageDepth, intendedDepth);

    // get likelihood of keeping each read
    double *probs = computeReadProbsByLengthAndSecondMetric(readLengths, readFullLengths,
                                                            (int) stList_length(inputReads),
                                                            intendedDepth, (int) chunkSize);

    // keep some ratio of reads
    int64_t totalKeptNucleotides = 0;
    for (int64_t i = 0; i < stList_length(inputReads); i++) {
        BamChunkRead *bcr = stList_get(inputReads, i);
        if (st_random() < probs[i]) {
            stList_append(maintainedReads, stList_get(inputReads, i));
            stList_append(maintainedAlignments, stList_get(inputAlignments, i));
            totalKeptNucleotides += bcr->rleRead->length;
        } else {
            stList_append(discardedReads, stList_get(inputReads, i));
            stList_append(discardedAlignments, stList_get(inputAlignments, i));
        }
    }
    st_logInfo(" %s Downsampled chunk via full read length to average depth %.2fx (expected %dx)\n",
               logIdentifier, 1.0 * totalKeptNucleotides / chunkSize, intendedDepth);

    free(probs);
    free(readLengths);
    free(readFullLengths);
    free(logIdentifier);
    return TRUE;
}



bool downsampleBamChunkReadWithVcfEntrySubstringsViaFullReadLengthLikelihood(int64_t intendedDepth,
                                                                             stList *chunkVcfEntries,
                                                                             stList *inputReads,
                                                                             stList *maintainedReads,
                                                                             stList *discardedReads) {

    // calculate depth, generate bcrwhs list
    int64_t totalVcfEntries = 0;

    int *readLengths = st_calloc(stList_length(inputReads), sizeof(int));
    int *readFullLengths = st_calloc(stList_length(inputReads), sizeof(int));

    for (int64_t i = 0; i < stList_length(inputReads); i++) {
        BamChunkRead *bcr = stList_get(inputReads, i);
        BamChunkReadVcfEntrySubstrings *bcrves = bcr->bamChunkReadVcfEntrySubstrings;
        assert(bcrves != NULL);
        totalVcfEntries += (int) stList_length(bcrves->vcfEntries);
        readLengths[i] = (int) stList_length(bcrves->vcfEntries);
        // give reads a pseudocount
        readFullLengths[i] = (int) bcr->fullReadLength;
    }
    int64_t chunkSize = stList_length(chunkVcfEntries);
    double averageDepth = 1.0 * totalVcfEntries / chunkSize;

    // do we need to downsample?
    if (averageDepth < intendedDepth) {
        free(readLengths);
        free(readFullLengths);
        return FALSE;
    }
    char *logIdentifier = getLogIdentifier();

    // is there something wrong with this chunk?
    if (chunkSize == 0 || totalVcfEntries == 0) {
        st_logInfo(" %s Downsampling all reads in chunk with %"PRId64" reads (%"PRId64" incoming filtered), as it has %"PRId64" spanned variants (chunk has %"PRId64")\n",
                   logIdentifier, stList_length(inputReads), stList_length(discardedReads), totalVcfEntries, chunkSize);

        // "filter" everything and return
        for (int64_t i = 0; i < stList_length(inputReads); i++) {
            stList_append(discardedReads, stList_get(inputReads, i));
        }
        free(readLengths);
        free(readFullLengths);
        free(logIdentifier);
        return TRUE;
    }

    // we do need to downsample
    st_logInfo(" %s Downsampling chunk with %"PRId64" reads via spanned variant count with average depth %.2fx to %dx.\n",
               logIdentifier, stList_length(inputReads), averageDepth, intendedDepth);

    // get likelihood of keeping each read
    double *probs = computeReadProbsByLengthAndSecondMetric(readLengths, readFullLengths,
                                                            (int) stList_length(inputReads),
                                                            intendedDepth, (int) chunkSize);

    // keep some ratio of reads
    int64_t totalKeptVariants = 0;
    for (int64_t i = 0; i < stList_length(inputReads); i++) {
        BamChunkRead *bcr = stList_get(inputReads, i);
        if (st_random() < probs[i]) {
            stList_append(maintainedReads, stList_get(inputReads, i));
            totalKeptVariants += stList_length(bcr->bamChunkReadVcfEntrySubstrings->readSubstrings);
        } else {
            stList_append(discardedReads, stList_get(inputReads, i));
        }
    }
    st_logInfo(" %s Downsampled chunk via spanned variant count with to average depth %.2fx (expected %dx)\n",
               logIdentifier, 1.0 * totalKeptVariants / chunkSize, intendedDepth);

    free(probs);
    free(readLengths);
    free(readFullLengths);
    free(logIdentifier);
    return TRUE;
}




void writeHaplotaggedBam(BamChunk *bamChunk, char *inputBamLocation, char *outputBamFileBase,
                         stSet *readsInH1, stSet *readsInH2, Params *params, char *logIdentifier) {
    /*
     * Write out sam files with reads in each split based on which haplotype partition they are in.
     */

    // Prep

    char *chunkIdentifier = NULL;
    if (bamChunk == NULL) {
        chunkIdentifier = stString_print("");
    } else if (bamChunk->chunkIdx == -1) {
        chunkIdentifier = stString_print(".%s-%"PRId64"-%"PRId64, bamChunk->refSeqName,
            bamChunk->chunkOverlapStart, bamChunk->chunkOverlapEnd);
    } else {
        chunkIdentifier = stString_print(".C%05"PRId64".%s-%"PRId64"-%"PRId64, bamChunk->chunkIdx, bamChunk->refSeqName,
                                         bamChunk->chunkOverlapStart, bamChunk->chunkOverlapEnd);
    }
    char *haplotaggedBamOutFile = stString_print("%s%s.haplotagged.bam", outputBamFileBase, chunkIdentifier);
    st_logInfo(" %s Writing BAM haplotype output to: %s \n", logIdentifier, haplotaggedBamOutFile);

    int64_t h1Count = 0;
    int64_t h2Count = 0;
    int64_t h0Count = 0;

    // input file management
    samFile *in = hts_open(inputBamLocation, "r");
    hts_idx_t *idx = NULL;
    hts_itr_multi_t *iter = NULL;
    if (in == NULL) {
        st_errAbort("ERROR: Cannot open bam file %s\n", inputBamLocation);
    }
    // bam index
    if ((idx = sam_index_load(in, inputBamLocation)) == 0) {
        st_errAbort("ERROR: Cannot open index for bam file %s\n", inputBamLocation);
    }
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();
    int r;

    // output file mangement
    samFile *out = NULL;
    if (outputBamFileBase != NULL) {
        out = hts_open(haplotaggedBamOutFile, "wb");
        r = sam_hdr_write(out, bamHdr);
    }
    FILE *readOut = NULL;

    // thread stuff
    htsThreadPool threadPool = {NULL, 0};
    # ifdef _OPENMP
    int tc = omp_get_max_threads();
    # else
    int tc = 1;
    # endif
    if (!(threadPool.pool = hts_tpool_init(tc))) {
        fprintf(stderr, "Error creating thread pool\n");
    }
    hts_set_opt(in, HTS_OPT_THREAD_POOL, &threadPool);
    hts_set_opt(out, HTS_OPT_THREAD_POOL, &threadPool);

    // prep for index (not entirely sure what all this does.  see samtools/sam_view.c
    int filter_state = ALL, filter_op = 0;
    int result;
    samview_settings_t settings = { .bed = NULL };
    char* region[1] = {};

    if (bamChunk != NULL) {
        region[0] = stString_print("%s:%d-%d", bamChunk->refSeqName, bamChunk->chunkOverlapStart,
                                   bamChunk->chunkOverlapEnd);
        settings.bed = bed_hash_regions(settings.bed, region, 0, 1,
                                        &filter_op); //insert(1) or filter out(0) the regions from the command line in the same hash table as the bed file
        if (!filter_op) filter_state = FILTERED;
        int regcount = 0;
        hts_reglist_t *reglist = bed_reglist(settings.bed, filter_state, &regcount);
        if (!reglist) {
            st_errAbort("ERROR: Could not create list of regions for read conversion");
        }
        if ((iter = sam_itr_regions(idx, bamHdr, reglist, regcount)) == 0) {
            st_errAbort("ERROR: Cannot open iterator for region %s for bam file %s\n", region[0], inputBamLocation);
        }
    }

    // fetch alignments (either all reads or without reads)
    while ((bamChunk == NULL ? sam_read1(in,bamHdr,aln) : sam_itr_multi_next(in, iter, aln)) >= 0) {
        // basic filtering (no read length, no cigar)
        if (aln->core.l_qseq <= 0) continue;
        if (aln->core.n_cigar == 0) continue;
        if ((aln->core.flag & (uint16_t) 0x4) != 0)
            continue; //unaligned
        if (!params->polishParams->includeSecondaryAlignments && (aln->core.flag & (uint16_t) 0x100) != 0)
            continue; //secondary
        if (!params->polishParams->includeSupplementaryAlignments && (aln->core.flag & (uint16_t) 0x800) != 0)
            continue; //supplementary

        char *readName = bam_get_qname(aln);
        bool has_tag = bam_aux_get(aln, "HP") != NULL;
        bool inH1 = stSet_search(readsInH1, readName);
        bool inH2 = stSet_search(readsInH2, readName);
        if (inH1 & !inH2) {
            int32_t ht = 1;
            if (has_tag) {
                bam_aux_update_int(aln, "HP", ht);
            } else {
                bam_aux_append(aln, "HP", 'i', sizeof(ht), (uint8_t*) &ht);
            }
            h1Count++;
        } else if (!inH1 & inH2) {
            int32_t ht = 2;
            if (has_tag) {
                bam_aux_update_int(aln, "HP", ht);
            } else {
                bam_aux_append(aln, "HP", 'i', sizeof(ht), (uint8_t*) &ht);
            }
            h2Count++;
        } else {
            int32_t ht = 0;
            if (has_tag) {
                bam_aux_update_int(aln, "HP", ht);
            } else {
                bam_aux_append(aln, "HP", 'i', sizeof(ht), (uint8_t*) &ht);
            }
            h0Count++;
        }
        r = sam_write1(out, bamHdr, aln);
    }
    st_logInfo(" %s Separated reads with divisions: H1 %"PRId64", H2 %"PRId64", and H0 %"PRId64"\n", logIdentifier,
            h1Count, h2Count, h0Count);

    // Cleanup
    if (bamChunk != NULL) {
        hts_itr_multi_destroy(iter);
        free(region[0]);
        bed_destroy(settings.bed);
    }
    hts_idx_destroy(idx);
    bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
    sam_close(in);
    sam_close(out);
    hts_tpool_destroy(threadPool.pool);
    free(chunkIdentifier);
    free(haplotaggedBamOutFile);
}


void poa_writeSupplementalChunkInformation2(char *outputBase, char *haplotypeIdentifier, int64_t chunkIdx,
                                            BamChunk *bamChunk, Poa *poa, stList *reads, Params *params,
                                            bool outputPoaDOT, bool outputPoaCSV, bool outputRepeatCounts) {

    if (outputPoaDOT) {
        char *outputPoaDotFilename = stString_print("%s.poa.C%05"PRId64".%s-%"PRId64"-%"PRId64"%s.dot",
                                                    outputBase, chunkIdx, bamChunk->refSeqName, bamChunk->chunkOverlapStart, bamChunk->chunkOverlapEnd,
                                                    haplotypeIdentifier);
        FILE *outputPoaTsvFileHandle = safe_fopen(outputPoaDotFilename, "w");
        poa_printDOT(poa, outputPoaTsvFileHandle, reads);
        fclose(outputPoaTsvFileHandle);
        free(outputPoaDotFilename);
    }
    if (outputPoaCSV) {
        char *outputPoaCsvFilename = stString_print("%s.poa.C%05"PRId64".%s-%"PRId64"-%"PRId64"%s.csv",
                                                    outputBase, chunkIdx, bamChunk->refSeqName, bamChunk->chunkOverlapStart, bamChunk->chunkOverlapEnd,
                                                    haplotypeIdentifier);
        FILE *outputPoaCsvFileHandle = safe_fopen(outputPoaCsvFilename, "w");
        poa_printCSV(poa, outputPoaCsvFileHandle, reads, params->polishParams->repeatSubMatrix, 5);
        fclose(outputPoaCsvFileHandle);
        free(outputPoaCsvFilename);
    }
    if (outputRepeatCounts) {
        char *outputRepeatCountFilename = stString_print("%s.repeatCount.C%05"PRId64".%s-%"PRId64"-%"PRId64"%s.csv",
                                                         outputBase, chunkIdx, bamChunk->refSeqName, bamChunk->chunkOverlapStart, bamChunk->chunkOverlapEnd,
                                                         haplotypeIdentifier);
        FILE *outputRepeatCountFileHandle = safe_fopen(outputRepeatCountFilename, "w");
        poa_printRepeatCountsCSV(poa, outputRepeatCountFileHandle, reads);
        fclose(outputRepeatCountFileHandle);
        free(outputRepeatCountFilename);
    }
}

void poa_writeSupplementalChunkInformation(char *outputBase, int64_t chunkIdx,
                                           BamChunk *bamChunk, Poa *poa, stList *reads, Params *params,
                                           bool outputPoaDOT, bool outputPoaCSV, bool outputRepeatCounts) {
    poa_writeSupplementalChunkInformation2(outputBase, "", chunkIdx, bamChunk, poa, reads, params,
                                           outputPoaDOT, outputPoaCSV, outputRepeatCounts);
}

void poa_writeSupplementalChunkInformationDiploid(char *outputBase, int64_t chunkIdx,
                                                  BamChunk *bamChunk, stGenomeFragment *genomeFragment, Poa *poaH1, Poa *poaH2, stList *bamChunkReads,
                                                  stSet *readsInHap1, stSet *readsInHap2, Params *params, bool outputPoaDOT, bool outputPoaCSV,
                                                  bool outputRepeatCounts, bool outputHaplotypedReadIdCsv, bool outputHaplotypedBam, char *logIdentifier) {

    poa_writeSupplementalChunkInformation2(outputBase, ".hap1", chunkIdx, bamChunk, poaH1, bamChunkReads,
                                           params, outputPoaDOT, outputPoaCSV, outputRepeatCounts);
    poa_writeSupplementalChunkInformation2(outputBase, ".hap2", chunkIdx, bamChunk, poaH2, bamChunkReads,
                                           params, outputPoaDOT, outputPoaCSV, outputRepeatCounts);

    if (outputHaplotypedReadIdCsv) {
        char *readIdsHap1Filename = stString_print("%s.readIds.C%05"PRId64".%s-%"PRId64"-%"PRId64".hap1.csv",
                                                   outputBase, chunkIdx, bamChunk->refSeqName, bamChunk->chunkOverlapStart, bamChunk->chunkOverlapEnd);
        FILE *readIdsHap1File = safe_fopen(readIdsHap1Filename, "w");
        stGenomeFragment_printPartitionAsCSV(genomeFragment, readIdsHap1File, params->phaseParams, TRUE, NULL);
        fclose(readIdsHap1File);

        char *readIdsHap2Filename = stString_print("%s.readIds.C%05"PRId64".%s-%"PRId64"-%"PRId64".hap2.csv",
                                                   outputBase, chunkIdx, bamChunk->refSeqName, bamChunk->chunkOverlapStart, bamChunk->chunkOverlapEnd);
        FILE *readIdsHap2File = safe_fopen(readIdsHap2Filename, "w");
        stGenomeFragment_printPartitionAsCSV(genomeFragment, readIdsHap2File, params->phaseParams, FALSE, NULL);
        fclose(readIdsHap2File);

        free(readIdsHap1Filename);
        free(readIdsHap2Filename);
    }

    if (outputHaplotypedBam) {
        // setup
        stSet *readIdsInHap1 = bamChunkRead_to_readName(readsInHap1);
        stSet *readIdsInHap2 = bamChunkRead_to_readName(readsInHap2);

        // write it
        writeHaplotaggedBam(bamChunk, bamChunk->parent->bamFile, outputBase, readIdsInHap1, readIdsInHap2, params,
                            logIdentifier);

        // cleanup
        stSet_destruct(readIdsInHap1);
        stSet_destruct(readIdsInHap2);
    }

}

void saveStartingVcfEntries(stList *vcfEntries, stHash *currentVcfEntries, int64_t *nextVcfEntriesIndex, int64_t cigarIdxInRef,
        int64_t refCigarModification, int64_t firstNonSoftclipAlignedReadIdxInChunk, int64_t cigarIdxInSeq,
        int64_t startSoftclipAmount) {

    VcfEntry *nextVcfEntry = *nextVcfEntriesIndex < stList_length(vcfEntries) ?
                             stList_get(vcfEntries, *nextVcfEntriesIndex) : NULL;
    while (nextVcfEntry != NULL && nextVcfEntry->refAlnStart <= cigarIdxInRef + refCigarModification) {
        // in theory, we should always land exactly on the right position, never skip it
        assert( firstNonSoftclipAlignedReadIdxInChunk == cigarIdxInSeq ||
                nextVcfEntry->refAlnStart == cigarIdxInRef + refCigarModification);
        // save in list of vcf entries spanning current pos
        stHash_insert(currentVcfEntries, nextVcfEntry, (void*)cigarIdxInSeq + startSoftclipAmount);
        // increment (it is possible to have mulitple entries at same pos)
        (*nextVcfEntriesIndex)++;
        nextVcfEntry = *nextVcfEntriesIndex < stList_length(vcfEntries) ?
                       stList_get(vcfEntries, *nextVcfEntriesIndex) : NULL;
    }
}


void saveFinishedVcfEntries(stHash *currentVcfEntries, int64_t currentAlnRefPos, int64_t cigarIdxInSeq,
        int64_t startSoftclipAmount, const uint8_t *seqBits, const uint8_t *qualBits, BamChunkReadVcfEntrySubstrings *bcrves,
        bool endOfRead) {

    // remove old vcf entries
    stList *noLongerCurrentVcfEntries = stList_construct();
    stHashIterator *currVcfEntryItor = stHash_getIterator(currentVcfEntries);
    VcfEntry *currVcfEntry = NULL;
    while ((currVcfEntry = stHash_getNext(currVcfEntryItor)) != NULL) {
        if (endOfRead || currVcfEntry->refAlnStopIncl <= currentAlnRefPos) {
            // in theory we should always land exactly on the end pos
            assert(endOfRead || currVcfEntry->refAlnStopIncl == currentAlnRefPos);

            // read start and stop pos
            int64_t seqStartPos = (int64_t) stHash_search(currentVcfEntries, currVcfEntry);
            int64_t seqEndPos = cigarIdxInSeq + startSoftclipAmount;

            // unsaved termination cases: delete over full variant, or end of read and haven't gotten to variant pos
            int64_t seqLen = seqEndPos - seqStartPos;
            if (seqLen == 0 || (endOfRead && currentAlnRefPos < currVcfEntry->refPos)) {
                stList_append(noLongerCurrentVcfEntries, currVcfEntry);
                continue;
            }

            // get sequence
            char *seq = st_calloc(seqLen + 1, sizeof(char));
            int64_t idxInOutputSeq = 0;
            int64_t idxInBamRead = seqStartPos;
            while (idxInBamRead < seqEndPos) {
                seq[idxInOutputSeq] = seq_nt16_str[bam_seqi(seqBits, idxInBamRead)];
                idxInBamRead++;
                idxInOutputSeq++;
            }
            seq[seqLen] = '\0';

            // get sequence qualities (if exists)
            uint8_t *qual = NULL;
            if (qualBits[0] != 0xff) { //inital score of 255 means qual scores are unavailable
                idxInOutputSeq = 0;
                idxInBamRead = seqStartPos;
                qual = st_calloc(seqLen, sizeof(uint8_t));
                while (idxInBamRead < seqEndPos) {
                    qual[idxInOutputSeq] = qualBits[idxInBamRead];
                    idxInBamRead++;
                    idxInOutputSeq++;
                }
                assert(idxInOutputSeq == strlen(seq));
            }

            // save it
            bamChunkReadVcfEntrySubstrings_saveSubstring(bcrves, seq, qual, currVcfEntry);
            stList_append(noLongerCurrentVcfEntries, currVcfEntry);
        }
    }

    // remove any vcf entries we're finished with
    for (int64_t r = 0; r < stList_length(noLongerCurrentVcfEntries); r++) {
        stHash_remove(currentVcfEntries, stList_get(noLongerCurrentVcfEntries, r));
    }

    // cleanup
    stHash_destructIterator(currVcfEntryItor);
    stList_destruct(noLongerCurrentVcfEntries);
}



uint32_t extractReadSubstringsAtVariantPositions(BamChunk *bamChunk, stList *vcfEntries, stList *reads,
        stList *filteredReads, PolishParams *polishParams) {

    // sanity check
    assert(stList_length(reads) == 0);

    // prep
    int64_t chunkOverlapStart = bamChunk->chunkOverlapStart;
    int64_t chunkOverlapEnd = bamChunk->chunkOverlapEnd;
    char *bamFile = bamChunk->parent->bamFile;
    char *contig = bamChunk->refSeqName;
    uint32_t savedAlignments = 0;

    // prep for index (not entirely sure what all this does.  see samtools/sam_view.c
    int filter_state = ALL, filter_op = 0;
    int result;
    samview_settings_t settings = {.bed = NULL};
    char *region[1] = {};
    region[0] = stString_print("%s:%d-%d", bamChunk->refSeqName, bamChunk->chunkOverlapStart,
                               bamChunk->chunkOverlapEnd);
    settings.bed = bed_hash_regions(settings.bed, region, 0, 1,
                                    &filter_op); //insert(1) or filter out(0) the regions from the command line in the same hash table as the bed file
    if (!filter_op) filter_state = FILTERED;
    int regcount = 0;
    hts_reglist_t *reglist = bed_reglist(settings.bed, filter_state, &regcount);
    if (!reglist) {
        st_errAbort("ERROR: Could not create list of regions for read conversion");
    }

    // file initialization
    samFile *in = NULL;
    hts_idx_t *idx = NULL;
    hts_itr_multi_t *iter = NULL;
    // bam file
    if ((in = hts_open(bamFile, "r")) == 0) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamFile);
    }
    // bam index
    if ((idx = sam_index_load(in, bamFile)) == 0) {
        st_errAbort("ERROR: Cannot open index for bam file %s\n", bamFile);
    }
    // header  //todo samFile *in = hts_open(bamFile, "r");
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    // read object
    bam1_t *aln = bam_init1();
    // iterator for region
    if ((iter = sam_itr_regions(idx, bamHdr, reglist, regcount)) == 0) {
        st_errAbort("ERROR: Cannot open iterator for region %s for bam file %s\n", region[0], bamFile);
    }

    // fetch alignments //todo while(sam_read1(in,bamHdr,aln) > 0) {
    while ((result = sam_itr_multi_next(in, iter, aln)) >= 0) {
        bool filtered = FALSE;
        // basic filtering (no read length, no cigar)
        if (aln->core.l_qseq <= 0) continue;
        if (aln->core.n_cigar == 0) continue;
        if ((aln->core.flag & (uint16_t) 0x4) != 0)
            continue; //unaligned
        if (!polishParams->includeSecondaryAlignments && (aln->core.flag & (uint16_t) 0x100) != 0)
            continue; //secondary
        if (!polishParams->includeSupplementaryAlignments && (aln->core.flag & (uint16_t) 0x800) != 0)
            continue; //supplementary
        if (aln->core.qual < polishParams->filterAlignmentsWithMapQBelowThisThreshold) { //low mapping quality
            if (filteredReads == NULL) continue;
            filtered = TRUE;
        }

        // data
        char *chr = bamHdr->target_name[aln->core.tid];
        int64_t start_softclip = 0;
        int64_t end_softclip = 0;
        int64_t alnReadLength = getAlignedReadLength3(aln, &start_softclip, &end_softclip, FALSE);
        if (alnReadLength <= 0) continue;
        int64_t alnStartPos = aln->core.pos;
        int64_t alnEndPos = alnStartPos + alnReadLength;

        // does this belong in our chunk?
        assert(stString_eq(contig, chr));
        // excludes reads starting before nominal chunk start, and after chunk end
        if (alnStartPos >= bamChunk->chunkEnd) continue;
        if (alnEndPos <= bamChunk->chunkStart) continue;

        // get read data
        uint8_t *seqBits = bam_get_seq(aln);
        char *readName = bam_get_qname(aln);
        uint8_t *qualBits = bam_get_qual(aln);
        bool forwardStrand = !bam_is_rev(aln);

        // for vcf tracking at local level
        // +1 because vcf->refPos is in 1-based space, and this reference position is in 0-based
        int64_t nextVcfEntriesIndex = binarySearchVcfListForFirstIndexAfterRefPos(vcfEntries, alnStartPos - chunkOverlapStart + 1);
        if (nextVcfEntriesIndex == -1) continue; // all vcf entries are before this read's start
        stHash *currentVcfEntries = stHash_construct(); // current vcf entries to read start pos
        // the bcrves we will be populating with vcf substrings
        BamChunkReadVcfEntrySubstrings *bcrves = bamChunkReadVcfEntrySubstrings_construct();
        BamChunkRead *bcr = bamChunkRead_constructWithVcfEntrySubstrings(readName, forwardStrand, alnReadLength, bcrves);

        // get cigar and rep
        uint32_t *cigar = bam_get_cigar(aln);

        // Variables to keep track of position in sequence / cigar operations
        int64_t cig_idx = 0;
        int64_t currPosInOp = 0;
        int64_t cigarOp = -1;
        int64_t cigarNum = -1;
        int64_t cigarIdxInSeq = 0;
        int64_t cigarIdxInRef = alnStartPos;

        // positional modifications
        int64_t refCigarModification = -1 * chunkOverlapStart;

        // we need to calculate:
        //  a. where in the (potentially softclipped read) to start storing characters
        //  b. what the alignments are wrt those characters
        // so we track the first aligned character in the read (for a.) and what alignment modification to make (for b.)
        int64_t firstNonSoftclipAlignedReadIdxInChunk;

        // the handling changes based on softclip inclusion and where the chunk boundaries are
        if (alnStartPos < chunkOverlapStart) {
            // alignment spans chunkStart
            firstNonSoftclipAlignedReadIdxInChunk = -1;
        } else {
            // alignment starts after chunkStart
            firstNonSoftclipAlignedReadIdxInChunk = 0;
        }

        // start with any vcf entries that may coincide with the beginning of the read
        if (start_softclip == 0) {
            saveStartingVcfEntries(vcfEntries, currentVcfEntries, &nextVcfEntriesIndex, cigarIdxInRef,
                                   refCigarModification, firstNonSoftclipAlignedReadIdxInChunk,
                                   cigarIdxInSeq, start_softclip);
        }

        // track number of characters in aligned portion (will inform softclipping at end of read)
        int64_t alignedReadLength = 0;

        // iterate over cigar operations
        for (uint32_t i = 0; i <= alnReadLength; i++) {
            // handles cases where last alignment is an insert or last is match
            if (cig_idx == aln->core.n_cigar) break;

            // do we need the next cigar operation?
            if (currPosInOp == 0) {
                cigarOp = cigar[cig_idx] & BAM_CIGAR_MASK;
                cigarNum = cigar[cig_idx] >> BAM_CIGAR_SHIFT;
            }

            // handle current character
            if (cigarOp == BAM_CMATCH || cigarOp == BAM_CEQUAL || cigarOp == BAM_CDIFF) {
                if (cigarIdxInRef >= chunkOverlapStart && cigarIdxInRef < chunkOverlapEnd) {
                    alignedReadLength++;
                }
                cigarIdxInSeq++;
                cigarIdxInRef++;
            } else if (cigarOp == BAM_CDEL || cigarOp == BAM_CREF_SKIP) {
                //delete
                cigarIdxInRef++;
            } else if (cigarOp == BAM_CINS) {
                //insert
                cigarIdxInSeq++;
                if (cigarIdxInRef >= chunkOverlapStart && cigarIdxInRef < chunkOverlapEnd) {
                    alignedReadLength++;
                }
                i--;
            } else if (cigarOp == BAM_CSOFT_CLIP) {
                // nothing to do here. skip to next cigar operation
                currPosInOp = cigarNum - 1;
                i--;
            } else if (cigarOp == BAM_CHARD_CLIP || cigarOp == BAM_CPAD) {
                // nothing to do here. skip to next cigar operation
                currPosInOp = cigarNum - 1;
                i--;
            } else {
                st_logCritical("Unidentifiable cigar operation!\n");
            }

            // document read index in the chunk (for reads that span chunk boundary, used in read construction)
            if (firstNonSoftclipAlignedReadIdxInChunk < 0 && cigarIdxInRef >= chunkOverlapStart) {
                firstNonSoftclipAlignedReadIdxInChunk = cigarIdxInSeq;
            }

            //  add new vcf entries
            saveStartingVcfEntries(vcfEntries, currentVcfEntries, &nextVcfEntriesIndex, cigarIdxInRef,
                    refCigarModification, firstNonSoftclipAlignedReadIdxInChunk, cigarIdxInSeq, start_softclip);
            // remove old vcf entries
            saveFinishedVcfEntries(currentVcfEntries, cigarIdxInRef + refCigarModification,
                    cigarIdxInSeq, start_softclip, seqBits, qualBits, bcrves, FALSE);

            // have we finished this last cigar
            currPosInOp++;
            if (currPosInOp == cigarNum) {
                cig_idx++;
                currPosInOp = 0;
            }
        }

        // finish final vcf stuff
        saveFinishedVcfEntries(currentVcfEntries, cigarIdxInRef + refCigarModification, cigarIdxInSeq, start_softclip,
                seqBits, qualBits, bcrves, TRUE);
        assert(stHash_size(currentVcfEntries) == 0);

        // save
        stList_append(filtered ? filteredReads: reads, bcr);

        savedAlignments++;

        // cleanup
        stHash_destruct(currentVcfEntries);
    }
    // the status from "get reads from iterator"
    if (result < -1) {
        st_errAbort("ERROR: Retrieval of region %d failed due to truncated file or corrupt BAM index file\n",
                    iter->curr_tid);
    }

    // close it all down
    hts_itr_multi_destroy(iter);
    hts_idx_destroy(idx);
    free(region[0]);
    bed_destroy(settings.bed);
    bam_hdr_destroy(bamHdr);
    bam_destroy1(aln);
    sam_close(in);

    return savedAlignments;
}


char *getSequenceFromReference(char *fastaFile, char *contig, int64_t startPos, int64_t endPosExcl) {
    faidx_t *fai = fai_load_format(fastaFile, FAI_FASTA);
    if ( !fai ) {
        st_errAbort("[faidx] Could not load fai index of %s\n", fastaFile);
    }
    // the faidx api is 1 based, tolerates use of idx 0 (ignores it), and is inclusive
    // this effectively makes it 0-based and exclusive
    startPos++;
    char *regionStr = stString_print("%s:%"PRId64"-%"PRId64, contig, startPos, endPosExcl);
    int seqLen;
    char *seq = fai_fetch(fai, regionStr, &seqLen);

    // sanity check
    assert(seqLen == endPosExcl - startPos + 1); // this sequence is 1-based (but tolerates idx 0)
    assert(seqLen == strlen(seq)); // seq len returns size of seq w/ \0-termination

    // convert to upper
    for (int i = 0; i < seqLen; i++) {
        seq[i] = (char) toupper(seq[i]);
    }

    // close and return
    free(regionStr);
    fai_destroy(fai);
    return seq;
}
