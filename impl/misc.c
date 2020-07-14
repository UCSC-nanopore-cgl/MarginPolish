/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"
#include <sys/stat.h>
#include <sonLibListPrivate.h>
#include <helenFeatures.h>
#include <htsIntegration.h>

char *getTimeDescriptorFromSeconds(int64_t seconds) {
    int64_t minutes = (int64_t) (seconds / 60);
    int64_t hours = (int64_t) (minutes / 60);
    char *timeDescriptor;

    if (hours > 0) {
        timeDescriptor = stString_print("%"PRId64"h %"PRId64"m", hours,
                minutes - (hours * 60));
    } else if (minutes > 0) {
        timeDescriptor = stString_print("%"PRId64"m %"PRId64"s", minutes,
                seconds - (minutes * 60));
    } else {
        timeDescriptor = stString_print("%"PRId64"s", seconds);
    }
    return timeDescriptor;
}

stHash *parseReferenceSequences(char *referenceFastaFile) {
    /*
     * Get hash of reference sequence names in fasta to their sequences, doing some munging on the sequence names.
     */
    st_logCritical("> Parsing reference sequences from file: %s\n", referenceFastaFile);
    FILE *fh = fopen(referenceFastaFile, "r");
    stHash *referenceSequences = fastaReadToMap(fh);  //valgrind says blocks from this allocation are "still reachable"
    fclose(fh);
    // log names and transform (if necessary)
    stList *refSeqNames = stHash_getKeys(referenceSequences);
    int64_t origRefSeqLen = stList_length(refSeqNames);
    st_logDebug("\tReference contigs: \n");
    for (int64_t i = 0; i < origRefSeqLen; ++i) {
        char *fullRefSeqName = (char *) stList_get(refSeqNames, i);
        st_logDebug("\t\t%s\n", fullRefSeqName);
        char refSeqName[128] = "";
        if (sscanf(fullRefSeqName, "%s", refSeqName) == 1 && !stString_eq(fullRefSeqName, refSeqName)) {
            // this transformation is necessary for cases where the reference has metadata after the contig name:
            // >contig001 length=1000 date=1999-12-31
            char *newKey = stString_copy(refSeqName);
            char *refSeq = stHash_search(referenceSequences, fullRefSeqName);
            stHash_insert(referenceSequences, newKey, refSeq);
            stHash_removeAndFreeKey(referenceSequences, fullRefSeqName);
            st_logDebug("\t\t\t-> %s\n", newKey);
        }
    }
    stList_destruct(refSeqNames);

    return referenceSequences;
}

char *getFileBase(char *base, char *defawlt) {
    struct stat fileStat;
    int64_t rc = stat(base, &fileStat);
    if (S_ISDIR(fileStat.st_mode)) {
        if (optarg[strlen(base) - 1] == '/') optarg[strlen(base) - 1] = '\0';
        return stString_print("%s/%s", base, defawlt);
    } else {
        return stString_copy(base);
    }
}

RleString *bamChunk_getReferenceSubstring(BamChunk *bamChunk, stHash *referenceSequences, Params *params) {
    /*
     * Get corresponding substring of the reference for a given bamChunk.
     */
    char *fullReferenceString = stHash_search(referenceSequences, bamChunk->refSeqName);
    if (fullReferenceString == NULL) {
        st_logCritical("> ERROR: Reference sequence missing from reference map: %s \n", bamChunk->refSeqName);
        return NULL;
    }
    int64_t refLen = strlen(fullReferenceString);
    char *referenceString = stString_getSubString(fullReferenceString, bamChunk->chunkBoundaryStart,
                                                  (refLen < bamChunk->chunkBoundaryEnd ? refLen : bamChunk->chunkBoundaryEnd) - bamChunk->chunkBoundaryStart);

    int64_t subRefLen = strlen(referenceString);
    // TODO: Decide where the proper place to do this is
    for (int64_t i = 0; i < subRefLen; i++) {
        referenceString[i] = (char) toupper(referenceString[i]);
    }

    RleString *rleRef = params->polishParams->useRunLengthEncoding ?
                        rleString_construct(referenceString) : rleString_construct_no_rle(referenceString);

    free(referenceString);

    return rleRef;
}


uint64_t *getPaddedHaplotypeString(uint64_t *hap, stGenomeFragment *gf, BubbleGraph *bg, Params *params) {
    /*
     * Pads a haplotype string from the genome fragment to account for any missing prefix or suffix.
     */
    uint64_t *paddedHap = bubbleGraph_getConsensusPath(bg, params->polishParams);

    for(uint64_t i=0; i<gf->length; i++) {
        paddedHap[i+gf->refStart] = hap[i];
    }

    return paddedHap;
}

stSet *bamChunkRead_to_readName(stSet *bamChunkReads) {
    stSet *readNames = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    stSetIterator *itor = stSet_getIterator(bamChunkReads);
    BamChunkRead *bcr = NULL;
    while ((bcr = stSet_getNext(itor)) != NULL) {
        stSet_insert(readNames, stString_copy(bcr->readName));
    }
    stSet_destructIterator(itor);
    return readNames;
}

stList *copyListOfIntTuples(stList *toCopy) {
    stList *copy = stList_construct3(0, toCopy->destructElement);
    for (int64_t j = 0; j < stList_length(toCopy); j++) {
        stIntTuple *tupleFrom = stList_get(toCopy, j);
        int64_t n = stIntTuple_length(tupleFrom);
        stIntTuple *tupleTo = stIntTuple_constructN(n, &(tupleFrom[1]));
        stList_append(copy, tupleTo);

        //todo delete
        assert(stIntTuple_length(tupleTo) == n);
        for (int64_t k = 0; k < n; k++) {
            assert(stIntTuple_get(tupleFrom, k) == stIntTuple_get(tupleTo, k));
        }
    }
    return copy;
}

double toPhred(double prob) {
    return -10 * log10(prob);
}

double fromPhred(double phred) {
    return pow(10, phred / -10.0);
}


void removeReadsStartingAfterChunkEnd(BamChunk *bamChunk, stList *reads, stList *alignments, char *logIdentifier) {
    int64_t chunkEnd = bamChunk->chunkEnd - bamChunk->chunkBoundaryStart;
    stList *indicesToRemove = stList_construct();
    for (int64_t i = stList_length(reads) - 1; i >= 0; i--) {
        stList *alignment = stList_get(alignments, i);
        if (stList_length(alignment) == 0) continue;
        int64_t firstAlignPos = stIntTuple_get(stList_get(alignment, 0), 0);
        if (firstAlignPos >= chunkEnd) {
            stList_append(indicesToRemove, (void*) i);
        }
    }
    st_logInfo(" %s Removing %"PRIu64" of %"PRIu64" reads for starting after chunkEnd\n", logIdentifier,
            stList_length(indicesToRemove), stList_length(reads));
    for (int64_t i = 0; i < stList_length(indicesToRemove); i++) {
        int64_t indexToRemove = (int64_t) stList_get(indicesToRemove, i);
        bamChunkRead_destruct(stList_get(reads, indexToRemove));
        stList_destruct(stList_get(alignments, indexToRemove));
        stList_remove(reads, indexToRemove);
        stList_remove(alignments, indexToRemove);
    }
    stList_destruct(indicesToRemove);

}

void writePhasedReadInfoJSON(BamChunk *bamChunk, stList *primaryReads, stList *primaryAlignments, stList *filteredReads,
        stList *filteredAlignments, stSet *readsInHap1, stSet *readsInHap2, uint64_t *reference_rleToNonRleCoordMap,
        FILE *out) {


    fprintf(out, ",\n \"reads\": [");

    // get scores and save to appropriate sets
    for (int64_t i = 0; i < stList_length(primaryReads); i++) {
        BamChunkRead *read = stList_get(primaryReads, i);
        stList *alignment = stList_get(primaryAlignments, i);
        stIntTuple *firstAlign = stList_get(alignment, 0);
        stIntTuple *lastAlign = stList_get(alignment, stList_length(alignment) - 1);
        int hap = 0;
        if (stSet_search(readsInHap1, read)) {
            hap = 1;
        } else if (stSet_search(readsInHap2, read)) {
            hap = 2;
        }

        //write
        if (i != 0) {
            fprintf(out, ",");
        }
        fprintf(out, "\n  {\n");
        fprintf(out, "     \"name\": \"%s\",\n", read->readName);
        fprintf(out, "     \"strand\": \"%s\",\n", read->forwardStrand ? "+" : "-");
        fprintf(out, "     \"startPos\": %"PRId64",\n", bamChunk->chunkBoundaryStart + reference_rleToNonRleCoordMap[stIntTuple_get(firstAlign, 0)]);
        fprintf(out, "     \"endPos\": %"PRId64",\n", bamChunk->chunkBoundaryStart + reference_rleToNonRleCoordMap[stIntTuple_get(lastAlign, 0)]);
        fprintf(out, "     \"hap\": %d\n", hap);
        fprintf(out, "  }");
    }

    // filtered reads
    for (int64_t i = 0; i < stList_length(filteredReads); i++) {
        BamChunkRead *read = stList_get(filteredReads, i);
        stList *alignment = stList_get(filteredAlignments, i);
        stIntTuple *firstAlign = stList_get(alignment, 0);
        stIntTuple *lastAlign = stList_get(alignment, stList_length(alignment) - 1);
        int hap = 0;
        if (stSet_search(readsInHap1, read)) {
            hap = 1;
        } else if (stSet_search(readsInHap2, read)) {
            hap = 2;
        }

        //write
        fprintf(out, ",\n  {\n");
        fprintf(out, "     \"name\": \"%s\",\n", read->readName);
        fprintf(out, "     \"strand\": \"%s\",\n", read->forwardStrand ? "+" : "-");
        fprintf(out, "     \"startPos\": %"PRId64",\n", bamChunk->chunkBoundaryStart + reference_rleToNonRleCoordMap[stIntTuple_get(firstAlign, 0)]);
        fprintf(out, "     \"endPos\": %"PRId64",\n", bamChunk->chunkBoundaryStart + reference_rleToNonRleCoordMap[stIntTuple_get(lastAlign, 0)]);
        fprintf(out, "     \"hap\": %d\n", hap);
        fprintf(out, "  }");
    }

    // write
    if (out != NULL) {
        fprintf(out, "\n ]");
    }
}


stList *produceVcfEntriesFromBubbleGraph(BamChunk *bamChunk, BubbleGraph *bg, stHash *readsToPSeqs,
                                         stGenomeFragment *gF, double strandSkewThreshold,
                                         double readSkewThreshold) {

    stList *vcfEntries = stList_construct3(0, (void(*)(void*))vcfEntry_destruct);

    char *logIdentifier = getLogIdentifier();
    int64_t passes = 0;
    int64_t total = 0;
    int64_t failsStrandSkew = 0;
    int64_t failsReadSkew = 0;
    for (uint64_t i = 0; i < gF->length; i++) {

        // bubble and hap info
        Bubble *b = &bg->bubbles[gF->refStart + i];
        RleString *hap1 = b->alleles[gF->haplotypeString1[i]];
        RleString *hap2 = b->alleles[gF->haplotypeString2[i]];

        // we only care about hets
        if (hap1 == hap2) continue;

        // read info
        int64_t hap1AlleleNo = gF->haplotypeString1[i];
        int64_t hap2AlleleNo = gF->haplotypeString2[i];
        int64_t totalReads = 0;
        int64_t hap1Reads = 0;
        for (uint64_t j = 0; j < b->readNo; j++) {

            double readHap1Support = b->alleleReadSupports[hap1AlleleNo * b->readNo + j];
            double readHap2Support = b->alleleReadSupports[hap2AlleleNo * b->readNo + j];

            if (readHap1Support != readHap2Support) {
                totalReads++;
                if (readHap1Support > readHap2Support) {
                    hap1Reads++;
                }
            }
        }

        // bubble skew info
        double strandSkew = bubble_phasedStrandSkew(b, readsToPSeqs, gF);
        double readSupportSkew = binomialPValue(totalReads, hap1Reads); // 1-bpv = prob of having a less extreme

        bool pass = TRUE;
        if (strandSkew < strandSkewThreshold) {
            failsStrandSkew++;
            pass = FALSE;
        }
        if (readSupportSkew < readSkewThreshold) {
            failsReadSkew++;
            pass = FALSE;
        }

        // maybe save to vcf entry
        uint64_t bubblePos = (uint64_t ) (b->refStart + b->bubbleLength / 2);
        if (pass) {
            stList_append(vcfEntries, vcfEntry_construct(bamChunk->refSeqName, bubblePos, -1, -1, rleString_copy(hap1),
                                                         rleString_copy(hap2)));
            passes++;
        }
        total++;

        //loggit
        st_logInfo(" %s Bubble at %"PRIu64" has strand skew %.5f%s and read support skew %.5f%s (%2"PRId64":%-2"PRId64" ): %s\n",
                logIdentifier, bubblePos, strandSkew, strandSkew < strandSkewThreshold ? "*" : " ",
                readSupportSkew, readSupportSkew < readSkewThreshold ? "*" : " ", hap1Reads, totalReads-hap1Reads,
                pass ? "PASS" : "FAIL");
    }

    //loggit
    st_logInfo(" %s Kept %"PRId64" of %"PRId64" bubbles after quality filtering. Failures: %"PRId64" strand, %"PRId64" reads\n",
            logIdentifier, passes, total, failsStrandSkew, failsReadSkew);
    free(logIdentifier);
    return vcfEntries;
}