/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <memory.h>
#include <hashTableC.h>
#include <unistd.h>
#include <time.h>
#include "marginVersion.h"

#include "margin.h"
#include "htsIntegration.h"
#include "helenFeatures.h"


/*
 * Main functions
 */

void usage() {
    fprintf(stderr, "usage: marginVVcfPhase <ALIGN_BAM> <REFERENCE_FASTA> <VARIANT_VCF> <PARAMS> [options]\n");
    fprintf(stderr, "Version: %s \n\n", MARGIN_POLISH_VERSION_H);
    fprintf(stderr, "Tags reads in ALIGN_BAM using variants in VARIANT_VCF.\n");

    fprintf(stderr, "\nRequired arguments:\n");
    fprintf(stderr, "    BAM_FILE is the alignment of reads to the assembly (or reference).\n");
    fprintf(stderr, "    ASSEMBLY_FASTA is the reference sequence BAM file in fasta format.\n");
    fprintf(stderr, "    PARAMS is the file with marginPolish parameters.\n");

    fprintf(stderr, "\nDefault options:\n");
    fprintf(stderr, "    -h --help                : Print this help screen\n");
    fprintf(stderr, "    -a --logLevel            : Set the log level [default = info]\n");
# ifdef _OPENMP
    fprintf(stderr, "    -t --threads             : Set number of concurrent threads [default = 1]\n");
#endif
    fprintf(stderr, "    -o --outputBase          : Name to use for output files [default = 'output']\n");
    fprintf(stderr, "    -r --region              : If set, will only compute for given chromosomal region\n");
    fprintf(stderr, "                                 Format: chr:start_pos-end_pos (chr3:2000-3000)\n");
    fprintf(stderr, "    -p --depth               : Will override the downsampling depth set in PARAMS\n");
    fprintf(stderr, "    -k --tempFilesToDisk     : Write temporary files to disk (for --diploid or supplementary output)\n");

    fprintf(stderr, "\n");
}


int main(int argc, char *argv[]) {

    // Parameters / arguments
    char *logLevelString = stString_copy("critical");
    char *bamInFile = NULL;
    char *paramsFile = NULL;
    char *referenceFastaFile = NULL;
    char *outputBase = stString_copy("output");
    char *regionStr = NULL;
    char *vcfFile = NULL;
    int numThreads = 1;
    int64_t maxDepth = -1;
    bool inMemory = TRUE;

    if (argc < 4) {
        free(outputBase);
        free(logLevelString);
        usage();
        return 0;
    }

    bamInFile = stString_copy(argv[1]);
    referenceFastaFile = stString_copy(argv[2]);
    vcfFile = stString_copy(argv[3]);
    paramsFile = stString_copy(argv[4]);

    // Parse the options
    while (1) {
        static struct option long_options[] = {
                { "help", no_argument, 0, 'h' },
                { "logLevel", required_argument, 0, 'a' },
# ifdef _OPENMP
                { "threads", required_argument, 0, 't'},
#endif
                { "outputBase", required_argument, 0, 'o'},
                { "region", required_argument, 0, 'r'},
                { "depth", required_argument, 0, 'p'},
                { "tempFilesToDisk", no_argument, 0, 'k'},
                { 0, 0, 0, 0 } };

        int option_index = 0;
        int key = getopt_long(argc-2, &argv[2], "ha:o:p:t:r:k", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
        case 'a':
            free(logLevelString);
            logLevelString = stString_copy(optarg);
            break;
        case 'h':
            usage();
            return 0;
        case 'o':
            free(outputBase);
            outputBase = getFileBase(optarg, "output");
            break;
        case 'r':
            regionStr = stString_copy(optarg);
            break;
        case 'p':
            maxDepth = atoi(optarg);
            if (maxDepth < 0) {
                st_errAbort("Invalid maxDepth: %s", optarg);
            }
            break;
        case 't':
            numThreads = atoi(optarg);
            if (numThreads <= 0) {
                st_errAbort("Invalid thread count: %d", numThreads);
            }
            break;
        case 'k':
            inMemory = FALSE;
            break;
        default:
            usage();
            free(outputBase);
            free(logLevelString);
            free(bamInFile);
            free(referenceFastaFile);
            free(paramsFile);
            return 0;
        }
    }

    // sanity check (verify files exist)
    if (access(bamInFile, R_OK) != 0) {
        st_errAbort("Could not read from input bam file: %s\n", bamInFile);
        char *idx = stString_print("%s.bai", bamInFile);
        if (access(idx, R_OK) != 0) {
            st_errAbort("BAM does not appear to be indexed: %s\n", bamInFile);
        }
        free(idx);
    }
    if (access(referenceFastaFile, R_OK) != 0) {
        st_errAbort("Could not read from reference fastafile: %s\n", referenceFastaFile);
    }
    if (access(paramsFile, R_OK) != 0) {
        st_errAbort("Could not read from vcf file: %s\n", vcfFile);
    }
    if (access(paramsFile, R_OK) != 0) {
        st_errAbort("Could not read from params file: %s\n", paramsFile);
    }

    // Initialization from arguments
    time_t startTime = time(NULL);
    st_setLogLevelFromString(logLevelString);
    free(logLevelString);
    if (st_getLogLevel() >= info) {
        st_setCallocDebug(true);
    }
# ifdef _OPENMP
    if (numThreads <= 0) {
        numThreads = 1;
    }
    omp_set_num_threads(numThreads);
    st_logCritical("Running OpenMP with %d threads.\n", omp_get_max_threads());
    # endif

    // Parse parameters
    st_logCritical("> Parsing model parameters from file: %s\n", paramsFile);
    Params *params = params_readParams(paramsFile);

    // update depth (if set)
    if (maxDepth >= 0) {
        st_logCritical("> Changing maxDepth parameter from %"PRId64" to %"PRId64"\n", params->polishParams->maxDepth,
                       maxDepth);
        params->polishParams->maxDepth = (uint64_t) maxDepth;
    }

    // Print a report of the parsed parameters
    if (st_getLogLevel() == debug) {
        params_printParameters(params, stderr);
    }

    // get vcf entries (if set)
    stHash *vcfEntries = NULL;
    if (vcfFile != NULL) {
        vcfEntries = parseVcf2(vcfFile, regionStr, params);
    }

    // get chunker for bam.  if regionStr is NULL, it will be ignored
    time_t chunkingStart = time(NULL);
    BamChunker *bamChunker = bamChunker_construct2(bamInFile, regionStr, params->polishParams, TRUE);
    st_logCritical(
            "> Set up bam chunker in %"PRId64"s with chunk size %i and overlap %i (for region=%s), resulting in %i total chunks\n",
            time(NULL) - chunkingStart, (int) bamChunker->chunkSize, (int) bamChunker->chunkBoundary,
            regionStr == NULL ? "all" : regionStr, bamChunker->chunkCount);
    if (bamChunker->chunkCount == 0) {
        st_errAbort("> Found no valid reads!\n");
    }

    // print chunk info
    char *outputChunksFile = stString_print("%s.chunks.csv", outputBase);
    FILE *chunksOut = safe_fopen(outputChunksFile, "w");
    for (int64_t i = 0; i < bamChunker->chunkCount; i++) {
        BamChunk *c = stList_get(bamChunker->chunks, i);
        fprintf(chunksOut, "%s,%"PRId64",%"PRId64",%"PRId64",%"PRId64"\n", c->refSeqName, c->chunkOverlapStart,
                c->chunkOverlapEnd, c->chunkStart, c->chunkEnd);
    }
    fclose(chunksOut);
    free(outputChunksFile);

    // output chunker tracks intermediate output files
    OutputChunkers *outputChunkers = outputChunkers_construct(numThreads, params, NULL, NULL, NULL, NULL,
            ".hap1", ".hap2", inMemory);

    // (may) need to shuffle chunks
    stList *chunkOrder = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < bamChunker->chunkCount; i++) {
        stList_append(chunkOrder, stIntTuple_construct1(i));
    }
    if (params->polishParams->shuffleChunks) {
        switch (params->polishParams->shuffleChunksMethod) {
            case SCM_SIZE_DESC:
                st_logCritical("> Ordering chunks by estimated depth\n");
                stList_sort2(chunkOrder, compareBamChunkDepthByIndexInList, bamChunker->chunks);
                stList_reverse(chunkOrder);
                break;
            case SCM_RANDOM:
                st_logCritical("> Randomly shuffling chunks\n");
                stList_shuffle(chunkOrder);
                break;
        }
    }

    // multiproccess the chunks, save to results
    st_logCritical("> Setup complete, beginning run\n");
    int64_t lastReportedPercentage = 0;
    time_t polishStartTime = time(NULL);

    # ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1)
    # endif
    for (int64_t i = 0; i < bamChunker->chunkCount; i++) {
        int64_t chunkIdx = stIntTuple_get(stList_get(chunkOrder, i), 0);
        // Time all chunks
        time_t chunkStartTime = time(NULL);

        // Get chunk
        BamChunk *bamChunk = bamChunker_getChunk(bamChunker, chunkIdx);

        // logging
        char *logIdentifier;
        bool logProgress = FALSE;
        int64_t currentPercentage = (int64_t) (100 * i / bamChunker->chunkCount);
        # ifdef _OPENMP
        int64_t threadIdx = omp_get_thread_num();
        logIdentifier = stString_print(" T%02d_C%05"PRId64, threadIdx, chunkIdx);
        if (threadIdx == 0) {
            if (currentPercentage != lastReportedPercentage) {
                logProgress = TRUE;
                lastReportedPercentage = currentPercentage;
            }
        }
        # else
        int64_t threadIdx = 0;
        logIdentifier = stString_copy("");
        if (currentPercentage != lastReportedPercentage) {
            logProgress = TRUE;
            lastReportedPercentage = currentPercentage;
        }
        # endif

        // prints percentage complete and estimated time remaining
        if (logProgress) {
            // log progress
            int64_t timeTaken = (int64_t) (time(NULL) - polishStartTime);
            int64_t secondsRemaining = (int64_t) floor(1.0 * timeTaken / currentPercentage * (100 - currentPercentage));
            char *timeDescriptor = (secondsRemaining == 0 && currentPercentage <= 50 ?
                                    stString_print("unknown") : getTimeDescriptorFromSeconds(secondsRemaining));
            st_logCritical("> Polishing %2"PRId64"%% complete (%"PRId64"/%"PRId64").  Estimated time remaining: %s\n",
                           currentPercentage, i, bamChunker->chunkCount, timeDescriptor);
            free(timeDescriptor);
        }

        // Get reference string for chunk of alignment
        char *chunkReference = getSequenceFromReference(referenceFastaFile, bamChunk->refSeqName,
                bamChunk->chunkOverlapStart, bamChunk->chunkOverlapEnd);
        st_logInfo(">%s Going to process a chunk for reference sequence: %s, starting at: %i and ending at: %i\n",
                   logIdentifier, bamChunk->refSeqName, (int) bamChunk->chunkOverlapStart, bamChunk->chunkOverlapEnd);

        // get VCF string
        stList *chunkVcfEntries = getVcfEntriesForRegion(vcfEntries, NULL, bamChunk->refSeqName,
                                                 bamChunk->chunkOverlapStart,  bamChunk->chunkOverlapEnd);
        updateVcfEntriesWithSubstringsAndPositions(chunkVcfEntries, chunkReference, strlen(chunkReference),
                FALSE, params);

        // Convert bam lines into corresponding reads and alignments
        st_logInfo(" %s Parsing input reads from file: %s\n", logIdentifier, bamInFile);
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *filteredReads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        extractReadSubstringsAtVariantPositions(bamChunk, chunkVcfEntries, reads, filteredReads, params->polishParams);

        // do downsampling if appropriate
        if (params->polishParams->maxDepth > 0) {
            // get downsampling structures
            stList *maintainedReads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);

            bool didDownsample = downsampleBamChunkReadWithVcfEntrySubstringsViaFullReadLengthLikelihood(
                    params->polishParams->maxDepth, chunkVcfEntries, reads, maintainedReads, filteredReads);

            // we need to destroy the discarded reads and structures
            if (didDownsample) {
                st_logInfo(" %s Downsampled from %"PRId64" to %"PRId64" reads\n", logIdentifier,
                           stList_length(reads), stList_length(maintainedReads));
                // still has all the old reads, need to not free these
                stList_setDestructor(reads, NULL);
                stList_destruct(reads);
                // and keep the filtered reads
                reads = maintainedReads;
            }
                // no downsampling, we just need to free the (empty) objects
            else {
                assert(stList_length(maintainedReads) == 0);
                stList_destruct(maintainedReads);
            }
        }

        time_t primaryPhasingStart = time(NULL);

        // iteratively find bubbles
        int64_t bubbleFindingIteration = 0;
        BubbleGraph *bg = NULL;
        stHash *readsToPSeqs = NULL;
        stSet *readsBelongingToHap1 = NULL, *readsBelongingToHap2 = NULL;
        stGenomeFragment *gf = NULL;
        stReference *ref = NULL;
        stList *vcfEntriesToBubbles = NULL;

        // Get the bubble graph representation
        bg =  bubbleGraph_constructFromVCFAndBamChunkReadVcfEntrySubstrings(reads, chunkVcfEntries, params,
                &vcfEntriesToBubbles);

        // Now make a POA for each of the haplotypes
        ref = bubbleGraph_getReference(bg, bamChunk->refSeqName, params);
        gf = bubbleGraph_phaseBubbleGraph(bg, ref, reads, params, &readsToPSeqs);

        stGenomeFragment_phaseBamChunkReads(gf, readsToPSeqs, reads, &readsBelongingToHap1, &readsBelongingToHap2,
                                            params->phaseParams);
        st_logInfo(" %s After phasing, of %i reads got %i reads partitioned into hap1 and %i reads partitioned "
                   "into hap2 (%i unphased)\n", logIdentifier, (int) stList_length(reads),
                   (int) stSet_size(readsBelongingToHap1), (int) stSet_size(readsBelongingToHap2),
                   (int) (stList_length(reads) - stSet_size(readsBelongingToHap1) -
                          stSet_size(readsBelongingToHap2)));


        st_logInfo(" %s Phased primary reads in %d sec\n", logIdentifier, time(NULL) - primaryPhasingStart);

        // should included filtered reads in output
        // get reads
        for (int64_t bcrIdx = 0; bcrIdx < stList_length(reads); bcrIdx++) {
            BamChunkRead *bcr = stList_get(reads, bcrIdx);
            if (!stSet_search(readsBelongingToHap1, bcr) && !stSet_search(readsBelongingToHap2, bcr)) {
                // was filtered in some form
                stList_append(filteredReads, bamChunkRead_constructCopy(bcr));
            }
        }
        st_logInfo(" %s Assigning %"PRId64" filtered reads to haplotypes\n", logIdentifier, stList_length(filteredReads));

        time_t filteredPhasingStart = time(NULL);

        bubbleGraph_partitionFilteredReadsFromVcfEntries(filteredReads, gf, bg, vcfEntriesToBubbles, readsBelongingToHap1,
                readsBelongingToHap2, params, logIdentifier);
        st_logInfo(" %s Partitioned filtered reads in %d sec.\n", logIdentifier, time(NULL) - filteredPhasingStart);


        // Output
        outputChunkers_processChunkSequencePhased(outputChunkers, threadIdx, chunkIdx, bamChunk->refSeqName,
                                                  NULL, NULL, reads, readsBelongingToHap1, readsBelongingToHap2, gf,
                                                  params);

        // Cleanup
        if (chunkVcfEntries != NULL) stList_destruct(chunkVcfEntries);
        stSet_destruct(readsBelongingToHap1);
        stSet_destruct(readsBelongingToHap2);
        bubbleGraph_destruct(bg);
        stGenomeFragment_destruct(gf);
        stReference_destruct(ref);
        stHash_destruct(readsToPSeqs);
        stList_destruct(vcfEntriesToBubbles);
        free(chunkReference);

        // report timing
        if (st_getLogLevel() >= info) {
            st_logInfo(">%s Chunk with ~%"PRId64" reads processed in %d sec\n",
                       logIdentifier, stList_length(reads) + stList_length(filteredReads), (int) (time(NULL) - chunkStartTime));
        }

        // final post-completion logging cleanup
        stList_destruct(reads);
        stList_destruct(filteredReads);
        free(logIdentifier);
    }

    // for writing haplotyped chunks
    stList *allReadIdsHap1 = stList_construct3(0, free);
    stList *allReadIdsHap2 = stList_construct3(0, free);

    // merge chunks
    time_t mergeStartTime = time(NULL);
    st_logCritical("> Starting merge\n");
    outputChunkers_stitchAndTrackReadIds(outputChunkers, TRUE, bamChunker->chunkCount, allReadIdsHap1, allReadIdsHap2);
    time_t mergeEndTime = time(NULL);
    char *tds = getTimeDescriptorFromSeconds((int) mergeEndTime - mergeStartTime);
    st_logCritical("> Merging took %s\n", tds);
    outputChunkers_destruct(outputChunkers);
    free(tds);
    tds = getTimeDescriptorFromSeconds((int) time(NULL) - mergeEndTime);
    st_logCritical("> Merge cleanup took %s\n", tds);
    free(tds);

    // maybe write final haplotyped bams
    time_t hapBamStart = time(NULL);
    st_logInfo("> Writing final haplotyped BAMs\n");

    stSet *allReadIdsForHaplotypingHap1 = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, NULL);
    stSet *allReadIdsForHaplotypingHap2 = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, NULL);

    for(int64_t i = 0; i < stList_length(allReadIdsHap1); i++) {
        stSet_insert(allReadIdsForHaplotypingHap1, stList_get(allReadIdsHap1, i));
    }
    for(int64_t i = 0; i < stList_length(allReadIdsHap2); i++) {
        stSet_insert(allReadIdsForHaplotypingHap2, stList_get(allReadIdsHap2, i));
    }

    // write it
    BamChunk *whbBamChunk = NULL;
    if (regionStr != NULL) {
        char regionContig[128] = "";
        int regionStart = 0;
        int regionEnd = 0;
        int scanRet = sscanf(regionStr, "%[^:]:%d-%d", regionContig, &regionStart, &regionEnd);
        if (scanRet != 3) {
            regionEnd = (int)((BamChunk*)stList_get(bamChunker->chunks, bamChunker->chunkCount -1))->chunkOverlapEnd;
        }
        whbBamChunk = bamChunk_construct2(regionContig, -1, regionStart, regionStart, regionEnd,
                regionEnd, 0, bamChunker);
    }
    writeHaplotaggedBam(whbBamChunk, bamChunker->bamFile, outputBase,
                        allReadIdsForHaplotypingHap1, allReadIdsForHaplotypingHap2, params, "");
    if (whbBamChunk != NULL) {
        bamChunk_destruct(whbBamChunk);
    }

    char *hapBamTDS = getTimeDescriptorFromSeconds(time(NULL) - hapBamStart);
    st_logCritical("> Wrote haplotyped bams in %s\n", hapBamTDS);

    // cleanup
    stSet_destruct(allReadIdsForHaplotypingHap1);
    stSet_destruct(allReadIdsForHaplotypingHap2);
    free(hapBamTDS);

    bamChunker_destruct(bamChunker);
    params_destruct(params);
    if (regionStr != NULL) free(regionStr);
    stList_destruct(chunkOrder);
    free(vcfFile);
    stHash_destruct(vcfEntries);
    if (allReadIdsHap1 != NULL) stList_destruct(allReadIdsHap1);
    if (allReadIdsHap2 != NULL) stList_destruct(allReadIdsHap2);
    free(outputBase);
    free(bamInFile);
    free(referenceFastaFile);
    free(paramsFile);

    // log completion
    char *timeDescriptor = getTimeDescriptorFromSeconds(time(NULL) - startTime);
    st_logCritical("> Finished phasing in %s.\n", timeDescriptor);
    free(timeDescriptor);

//    while(1); // Use this for testing for memory leaks

    return 0;
}

