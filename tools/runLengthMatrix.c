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
    fprintf(stderr, "usage: runLengthMatrix <ALIGN_BAM> <REFERENCE_FASTA> <PARAMS> [options]\n");
    fprintf(stderr, "Version: %s \n\n", MARGIN_POLISH_VERSION_H);
    fprintf(stderr, "Produces a run length matrix of reads in ALIGN_BAM to REFERENCE_FASTA.\n");

    fprintf(stderr, "\nRequired arguments:\n");
    fprintf(stderr, "    ALIGN_BAM is the alignment of reads to the reference.\n");
    fprintf(stderr, "    REFERENCE_FASTA is the reference sequence BAM file in fasta format.\n");
//    fprintf(stderr, "    VARIANT_VCF is the set of variants to use for phasing.\n");
    fprintf(stderr, "    PARAMS is the file with margin parameters.\n");
//
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
    fprintf(stderr, "    -l --maxRunLength        : Maximum run length (default 50)\n");

    fprintf(stderr, "\n");
}

int64_t charToNuclIdx(char nucl, bool forward) {
    switch (nucl) {
        case 'a':
        case 'A':
            return forward ? 0 : 3;
        case 'c':
        case 'C':
            return forward ? 1 : 2;
        case 'g':
        case 'G':
            return forward ? 2 : 1;
        case 't':
        case 'T':
            return forward ? 3 : 0;
        default:
            return -1;
    }
}


int64_t getRunLengthArrayIndex(int threadIdx, int64_t nuclIdx, uint64_t refRL, uint64_t readRL, int64_t maxRL) {
    // bad thread
    assert(threadIdx >= 0);
    assert(nuclIdx < 4);

    int64_t threadPos = threadIdx * 4 * maxRL * maxRL;
    int64_t nuclPos = nuclIdx * maxRL * maxRL;
    int64_t refRlPos = (refRL < maxRL ? refRL : maxRL - 1) * maxRL;
    int64_t readRlPos = (readRL < maxRL ? readRL : maxRL - 1);

    // bad nucl
    if (nuclPos < 0) return -1;

    return threadPos + nuclPos + refRlPos + readRlPos;
}


int64_t testRunLengthConstruction() {
    int64_t threadCount = 10;
    int64_t maxRunLenght = 10;
    int64_t nuclCount = 4;
    int64_t maxArraySize = threadCount * nuclCount * maxRunLenght * maxRunLenght;

    uint64_t *myArray = st_calloc(maxArraySize, sizeof(uint64_t));
    for (int thread = 0 ; thread < threadCount; thread++) {
        for (uint64_t refRl = 0; refRl < maxRunLenght; refRl++) {
            for (uint64_t readRl = 0; readRl < maxRunLenght; readRl++) {
                for (int64_t nucl = 0; nucl < nuclCount; nucl++) {
                    char nuc = (nucl==0 ? 'A' : (nucl==1 ? 'C' : (nucl==2 ? 'G' : 'T')));
                    for (int64_t strand = 0; strand < 2; strand++) {
                        int64_t idx = getRunLengthArrayIndex(thread, charToNuclIdx(nuc, strand == 0), refRl, readRl, maxRunLenght);
                        assert(idx < maxArraySize);
                        myArray[idx] += 1;
                    }
                }
            }
        }
    }

    for (int64_t i = 0; i < maxArraySize; i++) {
        assert(myArray[i] == 2);
    }

    free(myArray);
}


int main(int argc, char *argv[]) {

    // Parameters / arguments
    char *logLevelString = stString_copy("critical");
    char *bamInFile = NULL;
    char *referenceFastaFile = NULL;
    char *paramsFile = NULL;
    char *outputBase = stString_copy("output");
    char *regionStr = NULL;
    int numThreads = 1;
    int64_t maxDepth = -1;
    int64_t maxRunLengthExcl = 51;

    if (argc < 3) {
        free(outputBase);
        free(logLevelString);
        usage();
        return 0;
    }

    bamInFile = stString_copy(argv[1]);
    referenceFastaFile = stString_copy(argv[2]);
    paramsFile = stString_copy(argv[3]);

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
                { "maxRunLength", no_argument, 0, 'l'},
                { 0, 0, 0, 0 } };

        int option_index = 0;
        int key = getopt_long(argc-2, &argv[2], "ha:o:p:t:r:l:", long_options, &option_index);

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
        case 't':
            numThreads = atoi(optarg);
            if (numThreads <= 0) {
                st_errAbort("Invalid thread count: %d", numThreads);
            }
            break;
        case 'l':
            maxRunLengthExcl = atoi(optarg) +1;
            if (maxRunLengthExcl < 1) {
                st_errAbort("Invalid max run length: %s", optarg);
            }
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

    //testing
    testRunLengthConstruction();

    // Parse parameters
    st_logCritical("> Parsing model parameters from file: %s\n", paramsFile);
    Params *params = params_readParams(paramsFile);

    // parameter updates
    st_logInfo("  Setting chunkBoundary to 0\n");
    params->polishParams->chunkBoundary = 0;

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

    // this is the run length data we want
    int64_t totalSize = numThreads * 4 * maxRunLengthExcl * maxRunLengthExcl;
    uint64_t *runLengthDataForAllThreads = st_calloc(totalSize, sizeof(uint64_t));

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

        RleString *rleReference = bamChunk_getReferenceSubstring(bamChunk, referenceFastaFile, params);
        st_logInfo(">%s Going to process a chunk for reference sequence: %s, starting at: %i and ending at: %i\n",
                   logIdentifier, bamChunk->refSeqName, (int) bamChunk->chunkOverlapStart, bamChunk->chunkOverlapEnd);

        // Convert bam lines into corresponding reads and alignments
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *alignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
        stList *filteredReads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *filteredAlignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
        convertToReadsAndAlignments(bamChunk, rleReference, reads, alignments, params->polishParams);

        // do downsampling if appropriate
        if (params->polishParams->maxDepth > 0) {
            // get downsampling structures
            stList *maintainedReads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
            stList *maintainedAlignments = stList_construct3(0, (void (*)(void *)) stList_destruct);

            bool didDownsample = downsampleViaReadLikelihood(params->polishParams->maxDepth, bamChunk, reads,
                                                             alignments, maintainedReads, maintainedAlignments,
                                                             filteredReads, filteredAlignments);

            // we need to destroy the discarded reads and structures
            if (didDownsample) {
                st_logInfo(" %s Downsampled from %"PRId64" to %"PRId64" reads\n", logIdentifier,
                           stList_length(reads), stList_length(maintainedReads));
                // still has all the old reads, need to not free these
                stList_setDestructor(reads, NULL);
                stList_setDestructor(alignments, NULL);
                stList_destruct(reads);
                stList_destruct(alignments);
                // and keep the filtered reads
                reads = maintainedReads;
                alignments = maintainedAlignments;
            }
                // no downsampling, we just need to free the (empty) objects
            else {
                assert(stList_length(maintainedReads) == 0);
                assert(stList_length(maintainedAlignments) == 0);
                stList_destruct(maintainedReads);
                stList_destruct(maintainedAlignments);
            }
        }

        // prep for polishing
        Poa *poa = NULL; // The poa alignment

        // Generate partial order alignment (POA) (destroys rleAlignments in the process)
        poa = poa_realignOnlyAnchorAlignments(reads, alignments, rleReference, params->polishParams);

        for (int64_t pos = 1; pos < stList_length(poa->nodes); pos++) {
            PoaNode *node = stList_get(poa->nodes, pos);
            char refNucl = node->base;
            uint64_t refRL = node->repeatCount;
            for (int64_t o = 0; o < stList_length(node->observations); o++) {
                PoaBaseObservation *obs = stList_get(node->observations, o);
                BamChunkRead *read = stList_get(reads, obs->readNo);
                char readNucl = read->rleRead->rleString[obs->offset];
                uint64_t readRL = read->rleRead->repeatCounts[obs->offset];

                if (readNucl == refNucl) {
                    int64_t idx = getRunLengthArrayIndex(threadIdx, charToNuclIdx(readNucl, read->forwardStrand),
                            refRL, readRL, maxRunLengthExcl);
                    if (idx < 0) {
                        continue;
                    }
                    assert(idx < totalSize);
                    runLengthDataForAllThreads[idx] += 1;
                }
            }
        }


        // report timing
        if (st_getLogLevel() >= info) {
            st_logInfo(">%s Chunk with ~%"PRId64" reads processed in %d sec\n",
                       logIdentifier, stList_length(reads), (int) (time(NULL) - chunkStartTime));
        }

        // final post-completion logging cleanup
        poa_destruct(poa);
        rleString_destruct(rleReference);
        stList_destruct(reads);
        stList_destruct(alignments);
        stList_destruct(filteredReads);
        stList_destruct(filteredAlignments);
        free(logIdentifier);
    }

    st_logCritical("> Consolidating all run lengths\n");

    // condense all values
    uint64_t *condensedRunLengthArray = st_calloc(4 * maxRunLengthExcl * maxRunLengthExcl, sizeof(uint64_t));
    for (int t = 0; t < numThreads; t++) {
        for (int64_t nucl = 0; nucl < 4; nucl++) {
            for (uint64_t refRL = 1; refRL < maxRunLengthExcl; refRL++) {
                for (uint64_t readRL = 1; readRL < maxRunLengthExcl; readRL++) {
                    int64_t fullDataPos = getRunLengthArrayIndex(t, nucl, refRL, readRL, maxRunLengthExcl);
                    int64_t condensedPos = getRunLengthArrayIndex(0, nucl, refRL, readRL, maxRunLengthExcl);
                    assert(fullDataPos >= 0);
                    assert(condensedPos >= 0);
                    condensedRunLengthArray[condensedPos] += runLengthDataForAllThreads[fullDataPos];
                }
            }
        }
    }

    // printit
    char *countFilenameA = stString_print("%s.run_lengths.A.tsv", outputBase);
    FILE *countFileA = fopen(countFilenameA, "w");
    char *countFilenameC = stString_print("%s.run_lengths.C.tsv", outputBase);
    FILE *countFileC = fopen(countFilenameC, "w");
    char *countFilenameG = stString_print("%s.run_lengths.G.tsv", outputBase);
    FILE *countFileG = fopen(countFilenameG, "w");
    char *countFilenameT = stString_print("%s.run_lengths.T.tsv", outputBase);
    FILE *countFileT = fopen(countFilenameT, "w");
    if (countFileA == NULL || countFileC == NULL || countFileG == NULL || countFileT == NULL) {
        st_errAbort("Could not open output files for writing!", countFilenameA);
    } else {
        st_logCritical("> Writing counts to %s, %s, %s %s\n", countFilenameA, countFilenameC, countFilenameG, countFilenameT);
    }

    for (uint64_t refRL = 0; refRL < maxRunLengthExcl; refRL++) {
        for (uint64_t readRL = 0; readRL < maxRunLengthExcl; readRL++) {
            if (refRL == 0) {
                // header
                if (readRL == 0) {
                    fprintf(countFileA, "#ref_rl");
                    fprintf(countFileC, "#ref_rl");
                    fprintf(countFileG, "#ref_rl");
                    fprintf(countFileT, "#ref_rl");
                } else {
                    fprintf(countFileA, "read_%"PRId64"%s", readRL, readRL == maxRunLengthExcl - 1 ? "+" : "");
                    fprintf(countFileC, "read_%"PRId64"%s", readRL, readRL == maxRunLengthExcl - 1 ? "+" : "");
                    fprintf(countFileG, "read_%"PRId64"%s", readRL, readRL == maxRunLengthExcl - 1 ? "+" : "");
                    fprintf(countFileT, "read_%"PRId64"%s", readRL, readRL == maxRunLengthExcl - 1 ? "+" : "");
                }
            } else {
                if (readRL == 0) {
                    // header (ish)
                    fprintf(countFileA, "%"PRIu64, refRL);
                    fprintf(countFileC, "%"PRIu64, refRL);
                    fprintf(countFileG, "%"PRIu64, refRL);
                    fprintf(countFileT, "%"PRIu64, refRL);
                } else {
                    // data
                    int64_t condensedPosA = getRunLengthArrayIndex(0, charToNuclIdx('A', TRUE), refRL, readRL, maxRunLengthExcl);
                    int64_t condensedPosC = getRunLengthArrayIndex(0, charToNuclIdx('C', TRUE), refRL, readRL, maxRunLengthExcl);
                    int64_t condensedPosG = getRunLengthArrayIndex(0, charToNuclIdx('G', TRUE), refRL, readRL, maxRunLengthExcl);
                    int64_t condensedPosT = getRunLengthArrayIndex(0, charToNuclIdx('T', TRUE), refRL, readRL, maxRunLengthExcl);

                    uint64_t countA = condensedRunLengthArray[condensedPosA];
                    uint64_t countC = condensedRunLengthArray[condensedPosC];
                    uint64_t countG = condensedRunLengthArray[condensedPosG];
                    uint64_t countT = condensedRunLengthArray[condensedPosT];

                    fprintf(countFileA, "%"PRIu64, countA);
                    fprintf(countFileC, "%"PRIu64, countC);
                    fprintf(countFileG, "%"PRIu64, countG);
                    fprintf(countFileT, "%"PRIu64, countT);
                }
            }

            // increment
            if (readRL == maxRunLengthExcl - 1) {
                fprintf(countFileA, "\n");
                fprintf(countFileC, "\n");
                fprintf(countFileG, "\n");
                fprintf(countFileT, "\n");
            } else {
                fprintf(countFileA, "\t");
                fprintf(countFileC, "\t");
                fprintf(countFileG, "\t");
                fprintf(countFileT, "\t");
            }
        }
    }


    // close files
    fclose(countFileA);
    fclose(countFileC);
    fclose(countFileG);
    fclose(countFileT);

    // cleanup
    free(countFilenameA);
    free(countFilenameC);
    free(countFilenameG);
    free(countFilenameT);
    free(condensedRunLengthArray);
    free(runLengthDataForAllThreads);
    bamChunker_destruct(bamChunker);
    params_destruct(params);
    if (regionStr != NULL) free(regionStr);
    stList_destruct(chunkOrder);
    free(outputBase);
    free(bamInFile);
    free(referenceFastaFile);
    free(paramsFile);

    // log completion
    char *timeDescriptor = getTimeDescriptorFromSeconds(time(NULL) - startTime);
    st_logCritical("> Finished generating run length matrix in %s.\n", timeDescriptor);
    free(timeDescriptor);

//    while(1); // Use this for testing for memory leaks

    return 0;
}

