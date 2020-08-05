/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"
#include "bedidx.h"

#define DEFAULT_ALIGNMENT_SCORE 10

BamChunkRead *bamChunkRead_construct() {
    return bamChunkRead_construct2(NULL, NULL, NULL, TRUE, NULL);
}

BamChunkRead *bamChunkRead_construct2(char *readName, char *nucleotides,
                                      uint8_t *qualities, bool forwardStrand, bool useRunLengthEncoding) {
    return bamChunkRead_construct3(readName, nucleotides, qualities, forwardStrand, 0, useRunLengthEncoding);
}

BamChunkRead *bamChunkRead_construct3(char *readName, char *nucleotides, uint8_t *qualities, bool forwardStrand,
                                      int64_t fullReadLength, bool useRunLengthEncoding) {
    BamChunkRead *r = calloc(1, sizeof(BamChunkRead));
    r->readName = stString_copy(readName);
    r->forwardStrand = forwardStrand;
    r->fullReadLength = fullReadLength;
    assert(nucleotides != NULL);
    r->rleRead = useRunLengthEncoding ? rleString_construct(nucleotides) : rleString_construct_no_rle(nucleotides);
    if (qualities != NULL) {
        r->qualities = rleString_rleQualities(r->rleRead, qualities);
    }

    return r;
}

BamChunkRead *bamChunkRead_constructCopy(BamChunkRead *copy) {
    BamChunkRead *r = calloc(1, sizeof(BamChunkRead));
    r->readName = stString_copy(copy->readName);
    r->forwardStrand = copy->forwardStrand;
    r->fullReadLength = copy->fullReadLength;
    r->rleRead = rleString_copy(copy->rleRead);
    if (copy->qualities != NULL) {
        r->qualities = st_calloc(r->rleRead->length, sizeof(uint8_t));
        for (int64_t i = 0; i < r->rleRead->length; i++) {
            r->qualities[i] = copy->qualities[i];
        }
    }

    return r;
}

void bamChunkRead_destruct(BamChunkRead *r) {
    if (r->readName != NULL) free(r->readName);
    if (r->rleRead != NULL) rleString_destruct(r->rleRead);
    if (r->qualities != NULL) free(r->qualities);
    free(r);
}


//TODO remove?
//char *getMissingChunkSpacerSequence(int64_t overlapSize) {
//    //overlap size is 2x chunk boundary, we want 3x because the chunck boundary (2x) will all be overwritten)
//    int64_t spacerSize = (overlapSize == 0 ? 50 : (int64_t) ceil(overlapSize * 1.5));
//    char *missingChunkSpacer = st_calloc(spacerSize + 1, sizeof(char));
//    for (int64_t i = 0; i < spacerSize; i++) {
//        missingChunkSpacer[i] = 'N';
//    }
//    missingChunkSpacer[spacerSize] = '\0';
//}
//
//void alignAndSaveChunkConsensusSequence(stList *polishedReferenceStrings, char *currentChunk, Params *params, int64_t chunkIdx) {
//
//    char *previousChunk = stList_peek(polishedReferenceStrings);
//
//    // Trim the currrent and previous polished reference strings to remove overlap
//    int64_t prefixStringCropEnd, suffixStringCropStart;
//    int64_t overlapMatchWeight = removeOverlap(previousChunk, strlen(previousChunk), currentChunk, strlen(currentChunk),
//                                               params->polishParams->chunkBoundary * 2, params->polishParams,
//                                               &prefixStringCropEnd, &suffixStringCropStart);
//
//    // we have an overlap
//    if (overlapMatchWeight > 0) {
//        st_logInfo(
//                "    Removed overlap between neighbouring chunks at %"PRId64". "
//                "Approx overlap size: %i, overlap-match weight: %f, "
//                "left-trim: %i, right-trim: %i:\n", chunkIdx, params->polishParams->chunkBoundary * 2,
//                (float) overlapMatchWeight / PAIR_ALIGNMENT_PROB_1,
//                strlen(previousChunk) - prefixStringCropEnd, suffixStringCropStart);
//
//        // Crop the suffix of the previous chunk's polished reference string
//        previousChunk[prefixStringCropEnd] = '\0';
//
//        // Crop the the prefix of the current chunk's polished reference string
//        currentChunk = stString_copy(&(currentChunk[suffixStringCropStart]));
//
//    // no good alignment, likely missing chunks but have to be able to handle freaky situations also
//    } else {
//        if (strlen(currentChunk) == 0) {
//            // missing chunk
//            st_logInfo("    No overlap found for empty chunk at %"PRId64". Filling empty chunk with Ns.\n", chunkIdx);
//            currentChunk = getMissingChunkSpacerSequence(params->polishParams->chunkBoundary);
//        } else if (params->polishParams->chunkBoundary == 0) {
//            // poorly configured but could be done (freaky)
//            st_logInfo("    No overlap configured with non-empty (len %"PRId64") chunk at %"PRId64". \n",
//                       strlen(currentChunk), chunkIdx);
//            currentChunk = stString_copy(currentChunk);
//        } else {
//            // couldn't find an overlap (freaky)
//            st_logInfo("    No overlap found at %"PRId64". Filling Ns in stitch position.\n", chunkIdx);
//            stList_append(polishedReferenceStrings, stString_copy("NNNNNNNNNN"));
//            currentChunk = stString_copy(currentChunk);
//        }
//    }
//
//    // save it
//    stList_append(polishedReferenceStrings, currentChunk);
//}
//
//char *mergeContigChunks(char **chunks, int64_t startIdx, int64_t endIdxExclusive, Params *params) {
//
//    // merge chunks
//    stList *polishedReferenceStrings = stList_construct3(0, free); // The polished reference strings, one for each chunk
//
//    for (int64_t chunkIdx = startIdx; chunkIdx < endIdxExclusive; chunkIdx++) {
//        // Get chunk and polished
//        char *currentChunk = chunks[chunkIdx];
//
//        if (stList_length(polishedReferenceStrings) == 0) {
//            // we must copy the first one because the original and merged strings are freed separately
//            stList_append(polishedReferenceStrings, stString_copy(currentChunk));
//        } else {
//            alignAndSaveChunkConsensusSequence(polishedReferenceStrings, currentChunk, params, chunkIdx);
//            /*char *previousChunk = stList_peek(polishedReferenceStrings);
//
//            // Trim the currrent and previous polished reference strings to remove overlap
//            int64_t prefixStringCropEnd, suffixStringCropStart;
//            int64_t overlapMatchWeight = removeOverlap(previousChunk, strlen(previousChunk), currentChunk,
//                                                       strlen(currentChunk),
//                                                       overlap, params->polishParams,
//                                                       &prefixStringCropEnd, &suffixStringCropStart);
//
//            // we have an overlap
//            if (overlapMatchWeight > 0) {
//                st_logInfo(
//                        "    Removed overlap between neighbouring chunks at %"PRId64". "
//                        "Approx overlap size: %i, overlap-match weight: %f, "
//                        "left-trim: %i, right-trim: %i:\n", chunkIdx, overlap,
//                        (float) overlapMatchWeight / PAIR_ALIGNMENT_PROB_1,
//                        strlen(previousChunk) - prefixStringCropEnd, suffixStringCropStart);
//
//                // Crop the suffix of the previous chunk's polished reference string
//                previousChunk[prefixStringCropEnd] = '\0';
//
//                // Crop the the prefix of the current chunk's polished reference string
//                currentChunk = stString_copy(&(currentChunk[suffixStringCropStart]));
//
//                // no good alignment, likely missing chunks but have to be able to handle freaky situations also
//            } else {
//                if (strlen(currentChunk) == 0) {
//                    // missing chunk
//                    st_logInfo("    No overlap found for empty chunk at %"PRId64". Filling empty chunk with Ns.\n",
//                               chunkIdx);
//                    currentChunk = stString_copy(missingChunkSpacer);
//                } else if (overlap == 0) {
//                    // poorly configured but could be done (freaky)
//                    st_logInfo("    No overlap configured with non-empty (len %"PRId64") chunk at %"PRId64". \n",
//                               strlen(currentChunk), chunkIdx);
//                    currentChunk = stString_copy(currentChunk);
//                } else {
//                    // couldn't find an overlap (freaky)
//                    st_logInfo("    No overlap found at %"PRId64". Filling Ns in stitch position.\n", chunkIdx);
//                    stList_append(polishedReferenceStrings, stString_copy("NNNNNNNNNN"));
//                    currentChunk = stString_copy(currentChunk);
//                }
//            }*/
//        }
//    }
//
//    // finish
//    char *merged = stString_join2("", polishedReferenceStrings);
//    stList_destruct(polishedReferenceStrings);
//    return merged;
//}
//
//
//char *mergeContigChunksThreaded(char **chunks, int64_t startIdx, int64_t endIdxExclusive, int64_t numThreads,
//                                Params *params, char *referenceSequenceName) {
//
//    // special unthreaded case
//    if (numThreads == 1) return mergeContigChunks(chunks, startIdx, endIdxExclusive, params);
//
//    // divide into chunks
//    int64_t totalChunks = endIdxExclusive - startIdx;
//    int64_t chunksPerThread = (int64_t) ceil(1.0 * totalChunks / numThreads);
//    while (startIdx + chunksPerThread * (numThreads - 1) >= endIdxExclusive) { numThreads--; }
//    char **outputChunks = st_calloc(numThreads, sizeof(char *));
//
//    // multithread loop
//    st_logInfo(
//            "  Merging chunks for %s from (%"PRId64", %"PRId64"] with %"PRId64" chunks per thread on %"PRId64" threads \n",
//            referenceSequenceName, startIdx, endIdxExclusive, chunksPerThread, numThreads);
//#pragma omp parallel for schedule(static, 1)
//    for (int64_t thread = 0; thread < numThreads; thread++) {
//        int64_t threadedStartIdx = startIdx + chunksPerThread * thread;
//        int64_t threadedEndIdxExclusive = threadedStartIdx + chunksPerThread;
//        if (endIdxExclusive < threadedEndIdxExclusive) threadedEndIdxExclusive = endIdxExclusive;
//
//        outputChunks[thread] = mergeContigChunks(chunks, threadedStartIdx, threadedEndIdxExclusive, params);
//    }
//
//    // finish
//    char *contig = mergeContigChunks(outputChunks, 0, numThreads, params);
//    for (int64_t i = 0; i < numThreads; i++) {
//        free(outputChunks[i]);
//    }
//    free(outputChunks);
//    return contig;
//}
//
//
//char **mergeContigChunksDiploid(char **chunksH1, char **chunksH2, stSet **readsH1, stSet **readsH2,
//        stSet **lastReadsH1, stSet **lastReadsH2, int64_t startIdx, int64_t endIdxExclusive, Params *params) {
//
//    // merge chunks
//    stList *polishedReferenceStringsH1 = stList_construct3(0, free); // The polished reference strings, one for each chunk
//    stList *polishedReferenceStringsH2 = stList_construct3(0, free); // The polished reference strings, one for each chunk
//    stSet *previousReadsH1 = NULL;
//    stSet *previousReadsH2 = NULL;
//
//    for (int64_t chunkIdx = startIdx; chunkIdx < endIdxExclusive; chunkIdx++) {
//        // Get chunk and polished
//        char* currentChunkH1 = chunksH1[chunkIdx];
//        char* currentChunkH2 = chunksH2[chunkIdx];
//        stSet *currentReadsH1 = readsH1[chunkIdx];
//        stSet *currentReadsH2 = readsH2[chunkIdx];
//
//        if (stList_length(polishedReferenceStringsH1) == 0) {
//            // we must copy the first one because the original and merged strings are freed separately
//            stList_append(polishedReferenceStringsH1, stString_copy(currentChunkH1));
//            stList_append(polishedReferenceStringsH2, stString_copy(currentChunkH2));
//        } else {
//            char *previousChunkH1 = stList_peek(polishedReferenceStringsH1);
//            char *previousChunkH2 = stList_peek(polishedReferenceStringsH2);
//
//            //TODO
//            bool prevH1_to_currH1 = st_randomInt(0,2) == 0;
//
//            // just swap reads, always do h1 to h1 after this point
//            if (!prevH1_to_currH1) {
//                stSet *tmpReads = currentReadsH1;
//                currentReadsH1 = currentReadsH2;
//                currentReadsH2 = tmpReads;
//                char *tmpChunk = currentChunkH1;
//                currentChunkH1 = currentChunkH2;
//                currentChunkH2 = tmpChunk;
//            }
//
//            // save each segment separately
//            alignAndSaveChunkConsensusSequence(polishedReferenceStringsH1, currentChunkH1, params, chunkIdx);
//            alignAndSaveChunkConsensusSequence(polishedReferenceStringsH2, currentChunkH2, params, chunkIdx);
//        }
//
//        // remember previous reads per chunk
//        previousReadsH1 = currentReadsH1;
//        previousReadsH2 = currentReadsH2;
//    }
//
//    // sanity check
//    assert(previousReadsH1 != NULL);
//    assert(previousReadsH2 != NULL);
//
//    // finish
//    char **mergedHaps = st_calloc(2, sizeof(char*));
//    mergedHaps[0] = stString_join2("", polishedReferenceStringsH1);
//    mergedHaps[1] = stString_join2("", polishedReferenceStringsH2);
//    stList_destruct(polishedReferenceStringsH1);
//    stList_destruct(polishedReferenceStringsH2);
//    *lastReadsH1 = previousReadsH1;
//    *lastReadsH2 = previousReadsH2;
//    return mergedHaps;
//}
//
//
//char **mergeContigChunksDiploidThreaded(char **chunksH1, char **chunksH2, stSet **readsH1, stSet **readsH2,
//        int64_t startIdx, int64_t endIdxExclusive, int64_t numThreads, Params *params, char *referenceSequenceName) {
//
//    // special unthreaded case
//    stSet *unused1;
//    stSet *unused2;
//    if (numThreads == 1) return mergeContigChunksDiploid(chunksH1, chunksH2, readsH1, readsH2, &unused1, &unused2,
//            startIdx, endIdxExclusive, params);
//
//    // divide into chunks
//    int64_t totalChunks = endIdxExclusive - startIdx;
//    int64_t chunksPerThread = (int64_t) ceil(1.0 * totalChunks / numThreads);
//    while (startIdx + chunksPerThread * (numThreads - 1) >= endIdxExclusive) {numThreads--;}
//    char **mergedChunksH1 = st_calloc(numThreads, sizeof(char*));
//    char **mergedChunksH2 = st_calloc(numThreads, sizeof(char*));
//    stSet **mergedReadsH1 = st_calloc(numThreads, sizeof(stSet*));
//    stSet **mergedReadsH2 = st_calloc(numThreads, sizeof(stSet*));
//
//    // multithread loop
//    st_logInfo("  Merging chunks for %s from (%"PRId64", %"PRId64"] with %"PRId64" chunks per thread on %"PRId64" threads \n",
//               referenceSequenceName, startIdx, endIdxExclusive, chunksPerThread, numThreads);
//    #pragma omp parallel for schedule(static,1)
//    for (int64_t thread = 0; thread < numThreads; thread++) {
//        int64_t threadedStartIdx = startIdx + chunksPerThread * thread;
//        int64_t threadedEndIdxExclusive = threadedStartIdx + chunksPerThread;
//        if (endIdxExclusive < threadedEndIdxExclusive) threadedEndIdxExclusive = endIdxExclusive;
//
//        char **mergedThreadContigs = mergeContigChunksDiploid(chunksH1, chunksH2, readsH1, readsH2,
//                &(mergedReadsH1[thread]), &(mergedReadsH2[thread]), threadedStartIdx, threadedEndIdxExclusive, params);
//        mergedChunksH1[thread] = mergedThreadContigs[0];
//        mergedChunksH2[thread] = mergedThreadContigs[1];
//        free(mergedThreadContigs);
//    }
//
//    // final merging
//    char **finalContigHaplotypes = mergeContigChunksDiploid(mergedChunksH1, mergedChunksH2, mergedReadsH1, mergedReadsH2,
//            &unused1, &unused2, 0, numThreads, params);
//
//    // cleanup
//    for (int64_t i = 0; i < numThreads; i++) {
//        free(mergedChunksH1[i]);
//        free(mergedChunksH2[i]);
//    }
//    free(mergedChunksH1);
//    free(mergedChunksH1);
//    return finalContigHaplotypes;
//}
