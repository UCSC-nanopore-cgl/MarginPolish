//
// Created by Benedict Paten on 3/14/20.
//
// Code for stitching together "chunks" of inferred sequence
//

#include "margin.h"
#include "htsIntegration.h"



int chunkToStitch_cmp(ChunkToStitch *chunk1, ChunkToStitch *chunk2) {
    /*
     * Compares two chunks by their ordinal in the output order.
     */
    return chunk1->chunkOrdinal < chunk2->chunkOrdinal ? -1 : (chunk1->chunkOrdinal > chunk2->chunkOrdinal ? 1 : 0);
}

ChunkToStitch *chunkToStitch_construct(char *seqName, int64_t chunkOrdinal, bool phased,
        bool initRepeatCounts, bool initPoa) {
    ChunkToStitch *chunk = calloc(1, sizeof(ChunkToStitch));
    chunk->seqName = seqName;
    chunk->chunkOrdinal = chunkOrdinal;
    if (phased) {
        chunk->readsHap1Lines = stList_construct3(0, free);
        chunk->readsHap2Lines = stList_construct3(0, free);
    }
    if (initRepeatCounts) {
        chunk->repeatCountLinesHap1 = stList_construct3(0, free);
        if (phased) chunk->repeatCountLinesHap2 = stList_construct3(0, free);
    }
    if (initPoa) {
        chunk->poaHap1StringsLines = stList_construct3(0, free);
        if (phased) chunk->poaHap2StringsLines = stList_construct3(0, free);
    }

    return chunk;
}

void chunkToStitch_destruct(ChunkToStitch *chunkToStitch) {
    if (chunkToStitch->seqName != NULL) free(chunkToStitch->seqName);
    if (chunkToStitch->seqHap1 != NULL) free(chunkToStitch->seqHap1);

    // Second sequence and remaining fields are optional
    if (chunkToStitch->seqHap2 != NULL) {
        free(chunkToStitch->seqHap2);
    }

    if (chunkToStitch->poaHap1StringsLines != NULL) {
        stList_destruct(chunkToStitch->poaHap1StringsLines);
    }
    if (chunkToStitch->poaHap2StringsLines != NULL) {
        stList_destruct(chunkToStitch->poaHap2StringsLines);
    }

    if (chunkToStitch->repeatCountLinesHap1 != NULL) {
        stList_destruct(chunkToStitch->repeatCountLinesHap1);
    }
    if (chunkToStitch->repeatCountLinesHap2 != NULL) {
        stList_destruct(chunkToStitch->repeatCountLinesHap2);
    }

    if (chunkToStitch->readsHap1Lines != NULL) {
        stList_destruct(chunkToStitch->readsHap1Lines);
    }
    if (chunkToStitch->readsHap2Lines != NULL) {
        stList_destruct(chunkToStitch->readsHap2Lines);
    }

    free(chunkToStitch);
}

stList *readChunk(FILE *fh, char **seqName, int64_t *chunkOrdinal) {
    /*
     * Reads "chunks" from a file. Each chunk starts with a line formatted as:
     * SEQ_NAME,CHUNK_ORDINAL,LINE_NUMBER\n
     * where SEQ_NAME is the name of the string, containing any characters other than white space,
     * CHUNK_ORDINAL is the order of the chunk in the output,
     * and LINE_NUMBER is an integer >= 0 that gives the remaining number of lines in the chunk to read.
     * Returns the chunk as a list of lines, in order, and initializes the seqName argument to be a string representing
     * seqName.
     *
     * If the EOF is reached will return NULL, set seqName to NULL and chunkOrdinal to -1.
     *
     * Newlines characters are omitted from the ends of each line string.
     */
    char *headerLine = stFile_getLineFromFile(fh);
    if (headerLine == NULL) {
        *chunkOrdinal = -1;
        *seqName = NULL; // Seq seqName to NULL;
        return NULL;
    }

    stList *tokens = stString_splitByString(headerLine, ",");

    if (stList_length(tokens) != 3) {
        st_errAbort("Expected three tokens in header line, got %" PRIi64 "\n", stList_length(tokens));
    }

    *seqName = stList_removeFirst(tokens); // Set seqName

    *chunkOrdinal = strtol(stList_get(tokens, 0), NULL, 10); // Get chunk ordinal

    int64_t lineNo = strtol(stList_peek(tokens), NULL, 10); // Get line number

    stList *lines = stList_construct3(0, free);
    for (int64_t i = 0; i < lineNo; i++) {
        char *bodyLine = stFile_getLineFromFile(fh);
        if (bodyLine == NULL) {
            st_errAbort("Failed to read body line from chunk, line %" PRIi64 " of %" PRIi64 " lines\n", i, lineNo);
        }
        stList_append(lines, bodyLine);
    }

    // Cleanup
    stList_destruct(tokens);
    free(headerLine);

    return lines;
}

stList *readChunk2(FILE *fh, char *expectedSequenceName, int64_t expectedChunkOrdinal) {
    /*
     * As readChunk, but creates an error if no chunk is found or the read seqName is not equal to expectedSequenceName.
     */
    char *seqName;
    int64_t chunkOrdinal;
    stList *lines = readChunk(fh, &seqName, &chunkOrdinal);
    if (lines == NULL) {
        return NULL;
    }
    if (!stString_eq(seqName, expectedSequenceName)) {
        st_errAbort("Got an unexpected sequence name: %s in reading chunk (expected: %s)\n", seqName,
                    expectedSequenceName);
    }
    if (expectedChunkOrdinal != chunkOrdinal) {
        st_errAbort("Got an unexpected chunk ordinal (%" PRIi64 ") in reading chunk (expected: %" PRIi64 ")\n",
                    chunkOrdinal, expectedChunkOrdinal);
    }
    free(seqName);
    return lines;
}

bool chunkToStitch_readSequenceChunk(FILE *fh, ChunkToStitch *chunk, bool phased) {
    /*
     * Reads a "chunk" from a sequence containing file, adding it to chunk. Returns non-zero if no more chunks remain.
     */

    // Read the next set of lines from the file
    stList *lines = readChunk(fh, &chunk->seqName, &chunk->chunkOrdinal);

    // If we get nothing we have exhausted the file
    if (lines == NULL) {
        return 0;
    }

    chunk->seqHap1 = stString_join2("", lines);   // Concatenate the lines to make the sequence
    stList_destruct(lines); // Cleanup

    if (phased) {
        int64_t i;
        char *name;
        lines = readChunk(fh, &name, &i);
        if (lines == NULL) {
            st_errAbort("Error trying get alt sequence from chunk");
        }
        if (i != chunk->chunkOrdinal) {
            st_errAbort(
                    "Got an unexpected chunk ordinal (%" PRIi64 ") in reading second haplotype chunk (expected: %" PRIi64 ")\n",
                    i, chunk->chunkOrdinal);
        }
        if (!stString_eq(chunk->seqName, name)) {
            st_errAbort("Got an unexpected hap2 sequence name: %s in reading chunk (expected: %s)\n", name,
                        chunk->seqName);
        }
        chunk->seqHap2 = stString_join2("", lines);
        stList_destruct(lines); // Cleanup
        free(name);
    }

    return 1;
}

void chunkToStitch_readPoaChunk(FILE *fh, ChunkToStitch *chunk, bool phased) {
    /*
     * Reads a chunk from a POA containing file.
     */
    chunk->poaHap1StringsLines = readChunk2(fh, chunk->seqName, chunk->chunkOrdinal);
    if (phased) {
        chunk->poaHap2StringsLines = readChunk2(fh, chunk->seqName, chunk->chunkOrdinal);
    }
}

bool chunkToStitch_readReadPhasingChunk(FILE *fh, ChunkToStitch *chunk) {
    /*
     * Reads a read phasing chunk.  Returns true if chunk was found, returns FALSE if end of chunk
     */
    assert(chunk->readsHap1Lines == NULL);
    assert(chunk->readsHap2Lines == NULL);
    if (chunk->seqName == NULL) { // case where we are skipping fasta output for improved speed
        chunk->readsHap1Lines = readChunk(fh, &chunk->seqName, &chunk->chunkOrdinal);
        int64_t i;
        char *name;
        chunk->readsHap2Lines = readChunk(fh, &name, &i);
        if (chunk->readsHap2Lines != NULL) {
            if (i != chunk->chunkOrdinal) {
                st_errAbort("Got an unexpected chunk ordinal (%"PRIi64") in reading sequence lines (expected: %"PRIi64")\n",
                        i, chunk->chunkOrdinal);
            }
            if (!stString_eq(chunk->seqName, name)) {
                st_errAbort("Got an unexpected hap2 sequence name: %s in reading sequence lines (expected: %s)\n", name,
                            chunk->seqName);
            }
            free(name);
        }
    } else { // standard case, chunk is initialized from sequence file
        chunk->readsHap1Lines = readChunk2(fh, chunk->seqName, chunk->chunkOrdinal);
        chunk->readsHap2Lines = readChunk2(fh, chunk->seqName, chunk->chunkOrdinal);
    }

    if (chunk->readsHap1Lines == NULL ^ chunk->readsHap2Lines == NULL) {
        st_errAbort("Got reads for one chunk but not another! Expected chunk %"PRId64"\n", chunk->chunkOrdinal);
    }
    return chunk->readsHap1Lines != NULL;
}

void chunkToStitch_readRepeatCountChunk(FILE *fh, ChunkToStitch *chunk, bool phased) {
    /*
     * Reads repeat counts chunk
     */
    chunk->repeatCountLinesHap1 = readChunk2(fh, chunk->seqName, chunk->chunkOrdinal);
    if (phased) {
        chunk->repeatCountLinesHap2 = readChunk2(fh, chunk->seqName, chunk->chunkOrdinal);
    }
}

static void swap(void **a, void **b) {
    void *c = *a;
    *a = *b;
    *b = c;
}

static void addToHapReadsSeen(stHash *hapReads, stHash *otherHapReads, stHash *readsToAdd) {
    /*
     * Adds read names / probs from readsToAdd to hapReads that are not in otherHapReads. Cleans up readsToAdd.
     */
    stHashIterator *it = stHash_getIterator(readsToAdd);
    char *readName;
    while ((readName = stHash_getNext(it)) != NULL) {
        double *prob = stHash_search(readsToAdd, readName);

        /*
         * Check if in reads for other haplotype.
         * If it is, then if the prob of being in this haplotype is greater then
         * remove it from the other haplotype so it can be added to this haplotype,
         * otherwise do not add it to this haplotype.
         */
        double *pProb;
        if ((pProb = stHash_search(otherHapReads, readName)) != NULL) {
            if (*prob > *pProb) {
                free(stHash_removeAndFreeKey(otherHapReads, readName)); // Remove from otherHapReads
            } else {
                free(prob);
                continue;
            }
        }

        /*
         * Now add the read to this haplotype
         */
        if ((pProb = stHash_search(hapReads, readName)) == NULL) {
            stHash_insert(hapReads, stString_copy(readName), prob);
        } else if (*prob > *pProb) {
            free(stHash_removeAndFreeKey(hapReads, readName)); // Cleanup the prior entries
            stHash_insert(hapReads, stString_copy(readName), prob);
        } else {
            free(prob);
        }
    }
    stHash_destructIterator(it);
    stHash_setDestructValues(readsToAdd, NULL);
    stHash_destruct(readsToAdd);
}

stHash *getReadNames(stList *readPartitionLines) {
    /*
     * Parse the names of the reads from the lines of output representing the relative read phasing and return as a set
     * of strings.
     */
    stHash *readNames = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    for (int64_t i = 1; i < stList_length(readPartitionLines); i++) {
        char *line = stList_get(readPartitionLines, i);
        stList *tokens = stString_splitByString(line, ",");
        assert(stHash_search(readNames, stList_get(tokens, 0)) ==
               NULL); // Sanity check that read name is not present twice
        double *prob = st_calloc(1, sizeof(double));
        *prob = strtof(stList_peek(tokens), NULL); // Get the log prob of the read being in the partition
        stHash_insert(readNames, stList_removeFirst(tokens), prob); // First field is the read name
        stList_destruct(tokens);
    }
    return readNames;
}

int64_t sizeOfIntersection(stHash *pSet, stHash *nSet) {
    /*
     * Returns the number of keys in nSet also in pSet.
     */
    stHashIterator *it = stHash_getIterator(nSet);
    int64_t i = 0;
    char *readName;
    while ((readName = stHash_getNext(it)) != NULL) {
        if ((stHash_search(pSet, readName)) != NULL) {
            i++;
        }
    }
    stHash_destructIterator(it);
    return i;
}

int64_t sizeOfIntersectionWithNonNegativeValues(stHash *pSet, stHash *nSet) {
    /*
     * Returns the number of keys in nSet also in pSet.
     */
    stHashIterator *it = stHash_getIterator(nSet);
    int64_t i = 0;
    char *readName;
    while ((readName = stHash_getNext(it)) != NULL) {
        double nLikelihood = *((double*) stHash_search(nSet, readName));
        if (nLikelihood < 0)
            continue;
        if ((stHash_search(pSet, readName)) != NULL) {
            double pLikelihood = *((double*) stHash_search(pSet, readName));
            if (pLikelihood < 0)
                continue;
            i++;
        }
    }
    stHash_destructIterator(it);
    return i;
}

void chunkToStitch_phaseAdjacentChunks(ChunkToStitch *chunk, stHash *readsInHap1, stHash *readsInHap2) {
    /*
     * Phases chunk so that hap1 in chunk corresponds to hap1 in the prior chunks (as best as we can tell).
     */

    // Get the names of the reads in the different read sets
    stHash *chunkHap1Reads = getReadNames(chunk->readsHap1Lines);
    stHash *chunkHap2Reads = getReadNames(chunk->readsHap2Lines);

    // Calculate the intersection between reads shared between the chunks
    int64_t i = sizeOfIntersectionWithNonNegativeValues(readsInHap1, chunkHap1Reads);
    int64_t j = sizeOfIntersectionWithNonNegativeValues(readsInHap1, chunkHap2Reads);
    int64_t k = sizeOfIntersectionWithNonNegativeValues(readsInHap2, chunkHap1Reads);
    int64_t l = sizeOfIntersectionWithNonNegativeValues(readsInHap2, chunkHap2Reads);

    // Calculate support for the cis (keeping the current relative phasing) and the trans (switching the phasing) configurations
    int64_t cisPhase = i + l; // Number of reads consistently phased in cis configuration
    int64_t transPhase = j + k; // Number of reads consistently phased in trans configuration

    // Log the support for the phasing
    char *logIdentifier = getLogIdentifier();
    st_logInfo(" %s In phasing chunk %"PRId64" got %"
               PRIi64 " hap1 and %" PRIi64 " hap2 reads\n",
               logIdentifier, chunk->chunkOrdinal, stHash_size(chunkHap1Reads), stHash_size(chunkHap2Reads));
    st_logInfo(
            " %s Support for phasing cis-configuration,   Total: %" PRIi64 " (%f), %" PRIi64 " (%f%%) in h11-h21 intersection, %" PRIi64 " (%f%%) in h12-h22 intersection\n",
            logIdentifier, i + l, (double) (i + l) / (double) (i + j + k + l), i, (double) i / (double) (i + j), l,
            (double) l / (double) (i + j));
    st_logInfo(
            " %s Support for phasing trans-configuration, Total: %" PRIi64 " (%f), %" PRIi64 " (%f%%) in h11-h22 intersection, %" PRIi64 " (%f%%) in h12-h21 intersection\n",
            logIdentifier, j + k, (double) (j + k) / (double) (i + j + k + l), j, (double) j / (double) (k + l), k,
            (double) k / (double) (k + l));

    // Switch the relative phasing if the trans phase configuration has more support
    if (cisPhase < transPhase) {
        st_logInfo(" %s Flipping phase of chunk\n", logIdentifier);
        swap((void *) &chunk->seqHap1, (void *) &chunk->seqHap2);
        swap((void *) &chunk->poaHap1StringsLines, (void *) &chunk->poaHap2StringsLines);
        swap((void *) &chunk->readsHap1Lines, (void *) &chunk->readsHap2Lines);
        swap((void *) &chunk->repeatCountLinesHap1, (void *) &chunk->repeatCountLinesHap2);
        swap((void *) &chunkHap1Reads, (void *) &chunkHap2Reads);
    }

    //Remove duplicated reads from output
    addToHapReadsSeen(readsInHap1, readsInHap2, chunkHap1Reads);
    addToHapReadsSeen(readsInHap2, readsInHap1, chunkHap2Reads);

    // Cleanup
    free(logIdentifier);
}

int64_t removeOverlap(char *prefixString, int64_t prefixStringLength, char *suffixString, int64_t suffixStringLength,
                      int64_t approxOverlap, PolishParams *polishParams,
                      int64_t *prefixStringCropEnd, int64_t *suffixStringCropStart) {
    char *logIdentifier = getLogIdentifier();
    // Align the overlapping suffix of the prefixString and the prefix of the suffix string

    // Get coordinates of substrings to be aligned
    int64_t i = (prefixStringLength - approxOverlap) < 0 ? 0 : prefixStringLength - approxOverlap;
    int64_t j = approxOverlap < suffixStringLength ? approxOverlap : suffixStringLength;

    // Crop suffix
    char c = suffixString[j];
    suffixString[j] = '\0';

    // Symbol strings
    SymbolString sX = symbolString_construct(&(prefixString[i]), 0, strlen(&(prefixString[i])), polishParams->alphabet);
    SymbolString sY = symbolString_construct(suffixString, 0, strlen(suffixString), polishParams->alphabet);

    // Use default state machine for alignment
    StateMachine *sM = stateMachine3_constructNucleotide(threeState);

    // Get quick and dirty anchor pairs
    stList *anchorPairs = getKmerAlignmentAnchors(sX, sY, polishParams->p->diagonalExpansion);

    // Run the alignment
    stList *alignedPairs = getAlignedPairsUsingAnchors(sM, sX, sY, anchorPairs, polishParams->p, 1, 1);

    st_logInfo(" %s Removing overlap: got %"PRId64" anchor pairs and %"PRId64" aligned pairs for sequences of length x:%"PRId64", y:%"PRId64"\n",
            logIdentifier, stList_length(anchorPairs), stList_length(alignedPairs), sX.length, sY.length);

    // Cleanup
    symbolString_destruct(sX);
    symbolString_destruct(sY);
    stateMachine_destruct(sM);
    stList_destruct(anchorPairs);

    if (stList_length(alignedPairs) == 0 && st_getLogLevel() >= info) {
        st_logInfo(" %s Failed to find good overlap. Prefix-string: %s\n", logIdentifier, &(prefixString[i]));
        st_logInfo(" %s Failed to find good overlap. Suffix-string: %s\n", logIdentifier, suffixString);
    }

    // Remove the suffix crop
    suffixString[j] = c;

    // Pick the median point
    stIntTuple *maxPair = NULL;
    int64_t maxPairIdx = -1;
    for (int64_t k = 0; k < stList_length(alignedPairs); k++) {
        stIntTuple *aPair = stList_get(alignedPairs, k);
        if (maxPair == NULL || stIntTuple_get(aPair, 0) > stIntTuple_get(maxPair, 0)) {
            maxPairIdx = k;
            maxPair = aPair;
        }
    }
    if (maxPair == NULL) {
        st_logCritical(" %s Failed to find any aligned pairs between overlapping strings, not "
                       "doing any trimming (approx overlap: %i, len x: %i, len y: %i)\n", logIdentifier, approxOverlap,
                       prefixStringLength, suffixStringLength);
        *prefixStringCropEnd = prefixStringLength;
        *suffixStringCropStart = 0;
    } else {
        *prefixStringCropEnd = stIntTuple_get(maxPair, 1) + i; // Exclusive
        *suffixStringCropStart = stIntTuple_get(maxPair, 2);  // Inclusive
        st_logInfo(" %s Selecting best aligned pair at index %"PRId64" with pos p:%"PRId64"+%"PRId64", s:%"PRId64" with weight %"PRId64"\n",
                logIdentifier, maxPairIdx, stIntTuple_get(maxPair, 1), i, stIntTuple_get(maxPair, 2),
                stIntTuple_get(maxPair, 0));
    }

    int64_t overlapWeight = maxPair == NULL ? 0 : stIntTuple_get(maxPair, 0);

    stList_destruct(alignedPairs);
    free(logIdentifier);

    return overlapWeight;
}

void renumberCSVLines(stList *csvLines, int64_t index) {
    /*
     * Renumber the CSV lines in the output so that they are all sequential
     */
    for (int64_t i = 0; i < stList_length(csvLines); i++) { // Start from 0 as header line is already gone
        char *line = stList_get(csvLines, i);
        stList *tokens = stString_splitByString(line, ",");
        free(stList_get(tokens, 0)); // Cleanup the old index
        stList_set(tokens, 0, stString_print("%" PRIi64 "", index++));
        stList_set(csvLines, i, stString_join2(",", tokens));
        stList_destruct(tokens);
        free(line);
    }
}

void chunkToStitch_trimAdjacentChunks2(char **pSeq, char **seq,
                                       stList *pPoa, stList *poa, stList *pRepeatCounts, stList *repeatCounts,
                                       Params *params, int64_t *lengthOfSequenceOutputSoFar) {

    // for very fast case, where we don't write sequences
    if (*pSeq == NULL && *seq == NULL) {
        return;
    }

    // Convert to RLE space
    RleString *pSeqRle = params->polishParams->useRunLengthEncoding ?
                         rleString_construct(*pSeq) : rleString_construct_no_rle(*pSeq);
    RleString *seqRle = params->polishParams->useRunLengthEncoding ?
                        rleString_construct(*seq) : rleString_construct_no_rle(*seq);

    // Get the trim factor
    int64_t pSeqCropEnd = -1, seqCropStart = -1;
    int64_t overlapMatchWeight = removeOverlap(pSeqRle->rleString, pSeqRle->length, seqRle->rleString, seqRle->length,
                                               params->polishParams->chunkBoundary * 2,
                                               params->polishParams, &pSeqCropEnd, &seqCropStart);

    // Log
    char *logIdentifier = getLogIdentifier();
    st_logInfo(
            " %s Removing overlap between neighbouring chunks (in RLE space). Approx overlap size: %i, "
            "overlap-match weight: %f, left-trim: %i, right-trim: %i:\n", logIdentifier,
            (int) params->polishParams->chunkBoundary * 2,
            (float) overlapMatchWeight / PAIR_ALIGNMENT_PROB_1, pSeqRle->length - pSeqCropEnd, seqCropStart);

    // debug logging
    if (st_getLogLevel() >= info) {
        char *tmpSeq;
        if (pSeqRle->length > 17) {
            char ch = pSeqRle->rleString[8];
            pSeqRle->rleString[8] = '\0';
            tmpSeq = stString_print("%s...%s", pSeqRle->rleString, &(pSeqRle->rleString[pSeqRle->length - 8]));
            pSeqRle->rleString[8] = ch;
        } else {
            tmpSeq = stString_copy(pSeqRle->rleString);
        }
        st_logInfo(" %s pSeq:  pSeqCropEnd:%7"PRId64", LenRLE:%7"PRId64", LenRAW:%7"PRId64", seq: %s\n",
                   logIdentifier, pSeqCropEnd, pSeqRle->length, pSeqRle->nonRleLength, tmpSeq);
        free(tmpSeq);

        if (seqRle->length > 17) {
            char ch = seqRle->rleString[8];
            seqRle->rleString[8] = '\0';
            tmpSeq = stString_print("%s...%s", seqRle->rleString, &(seqRle->rleString[seqRle->length - 8]));
            seqRle->rleString[8] = ch;
        } else {
            tmpSeq = stString_copy(seqRle->rleString);
        }
        st_logInfo(" %s  seq: seqCropStart:%7"PRId64", LenRLE:%7"PRId64", LenRAW:%7"PRId64", seq: %s\n",
                   logIdentifier, seqCropStart, seqRle->length, seqRle->nonRleLength, tmpSeq);
        free(tmpSeq);
    }

    // Trim the sequences

    // Crop the suffix of the previous sequence
    RleString *pSeqRleCropped = rleString_copySubstring(pSeqRle, 0, pSeqCropEnd);
    free(*pSeq);
    *pSeq = rleString_expand(pSeqRleCropped);

    // Crop the the prefix of the current chunk's sequence
    RleString *seqRleCropped = rleString_copySubstring(seqRle, seqCropStart, seqRle->length - seqCropStart);
    free(*seq);
    *seq = rleString_expand(seqRleCropped);

    // Trim the remaining stuff
    int64_t suffixRleTrimLength = pSeqRle->length - pSeqCropEnd;
    *lengthOfSequenceOutputSoFar += pSeqCropEnd;

    // Poa
    if (poa != NULL) {
        stList_removeInterval(pPoa, stList_length(pPoa) - suffixRleTrimLength, suffixRleTrimLength);
        stList_removeInterval(poa, 0, seqCropStart + 2);
        renumberCSVLines(poa, 1 + *lengthOfSequenceOutputSoFar);
    }

    // Repeat counts
    if (repeatCounts != NULL) {
        stList_removeInterval(pRepeatCounts, stList_length(pRepeatCounts) - suffixRleTrimLength, suffixRleTrimLength);
        stList_removeInterval(repeatCounts, 0, seqCropStart + 2);
        renumberCSVLines(repeatCounts, 1 + *lengthOfSequenceOutputSoFar);
    }

    // Cleanup
    free(logIdentifier);
    rleString_destruct(pSeqRle);
    rleString_destruct(pSeqRleCropped);
    rleString_destruct(seqRle);
    rleString_destruct(seqRleCropped);
}

void chunkToStitch_trimAdjacentChunks(ChunkToStitch *pChunk, ChunkToStitch *chunk, Params *params,
                                      int64_t *lengthOfSequenceOutputSoFarHap1,
                                      int64_t *lengthOfSequenceOutputSoFarHap2) {
    /*
     * Trims the right end of pChunk and the left end of chunk so that they do not overlap, but are directly contiguous.
     */
    // Checks that they are part of the same sequence
    assert(stString_eq(pChunk->seqName, chunk->seqName));

    // Trim haplotype 1 sequences
    chunkToStitch_trimAdjacentChunks2(&pChunk->seqHap1, &chunk->seqHap1,
                                      pChunk->poaHap1StringsLines, chunk->poaHap1StringsLines,
                                      pChunk->repeatCountLinesHap1, chunk->repeatCountLinesHap1,
                                      params, lengthOfSequenceOutputSoFarHap1);

    // Trim haplotype 2 sequences, if it exists
    if (chunk->seqHap2 != NULL) {
        chunkToStitch_trimAdjacentChunks2(&pChunk->seqHap2, &chunk->seqHap2,
                                          pChunk->poaHap2StringsLines, chunk->poaHap2StringsLines,
                                          pChunk->repeatCountLinesHap2, chunk->repeatCountLinesHap2,
                                          params, lengthOfSequenceOutputSoFarHap2);
    }
}

/*
 * OutputChunker
 */

typedef struct _outputChunker {
    /*
     * Object for managing the output of a polished sequence.
     */
    bool useMemoryBuffers;
    // Sequence file
    bool outputSequence;
    char *outputSequenceFile; // This is either the name of the file or the location in memory of the file buffer
    FILE *outputSequenceFileHandle;
    size_t outputSequenceFileBufferSize; // Used if in memory
    // Poa file
    bool outputPoa;
    char *outputPoaFile;
    FILE *outputPoaFileHandle;
    size_t outputPoaFileBufferSize;
    // Repeat count file
    bool outputRepeatCounts;
    char *outputRepeatCountFile;
    FILE *outputRepeatCountFileHandle;
    size_t outputRepeatCountFileBufferSize;
    // Read partition file - this must be specified if phasing is to be performed, as the
    // information is needed for stitching
    bool outputReadPartition;
    char *outputReadPartitionFile;
    FILE *outputReadPartitionFileHandle;
    size_t outputRepeatPartitionFileBufferSize;

    Params *params;
} OutputChunker;

static FILE *open(bool output, char **file, size_t *outputBufferSize, bool inMemory, char *openStr) {
    if (!output) {
        return NULL;
    }
    if (inMemory) {
        if (stString_eq(openStr, "w")) {
            return open_memstream(file, outputBufferSize);
        }
        assert(stString_eq(openStr, "r"));
        return fmemopen(*file, (*outputBufferSize) + 1, "r");
    }
    return safe_fopen(*file, openStr);
}

void outputChunker_open(OutputChunker *outputChunker, char *openStr) {
    /*
     * Open the files.
     */
    outputChunker->outputSequenceFileHandle = open(outputChunker->outputSequence, &(outputChunker->outputSequenceFile),
                                                   &(outputChunker->outputSequenceFileBufferSize),
                                                   outputChunker->useMemoryBuffers, openStr);
    outputChunker->outputPoaFileHandle = open(outputChunker->outputPoa, &(outputChunker->outputPoaFile),
                                              &(outputChunker->outputPoaFileBufferSize),
                                              outputChunker->useMemoryBuffers, openStr);
    outputChunker->outputRepeatCountFileHandle = open(outputChunker->outputRepeatCounts,
                                                      &(outputChunker->outputRepeatCountFile),
                                                      &(outputChunker->outputRepeatCountFileBufferSize),
                                                      outputChunker->useMemoryBuffers, openStr);
    outputChunker->outputReadPartitionFileHandle = open(outputChunker->outputReadPartition,
                                                        &(outputChunker->outputReadPartitionFile),
                                                        &(outputChunker->outputRepeatPartitionFileBufferSize),
                                                        outputChunker->useMemoryBuffers, openStr);
}

OutputChunker *
outputChunker_construct(Params *params, char *outputSequenceFile, char *outputPoaFile, char *outputReadPartitionFile,
                        char *outputRepeatCountFile) {
    /*
     * Create an OutputChunker object, ready to write chunks of output to the given output files.
     */
    OutputChunker *outputChunker = st_calloc(1, sizeof(OutputChunker));

    // Initialize variables
    outputChunker->useMemoryBuffers = 0;
    outputChunker->outputSequence = outputSequenceFile != NULL;
    outputChunker->outputSequenceFile = outputSequenceFile;
    outputChunker->outputPoaFile = outputPoaFile;
    outputChunker->outputPoa = outputPoaFile != NULL;
    outputChunker->outputRepeatCountFile = outputRepeatCountFile;
    outputChunker->outputRepeatCounts = outputRepeatCountFile != NULL;
    outputChunker->outputReadPartitionFile = outputReadPartitionFile;
    outputChunker->outputReadPartition = outputReadPartitionFile != NULL;
    outputChunker->params = params;

    // Open files for writing
    outputChunker_open(outputChunker, "w");

    return outputChunker;
}

OutputChunker *
outputChunker_constructInMemory(Params *params, bool outputSequence,  bool outputPoaFile, bool outputReadPartitionFile,
                                bool outputRepeatCountFile) {
    /*
     * Create an OutputChunker object, ready to write chunks of output to the given output files.
     */
    OutputChunker *outputChunker = st_calloc(1, sizeof(OutputChunker));

    // Initialize variables
    outputChunker->useMemoryBuffers = 1;
    outputChunker->outputSequence = outputSequence;
    outputChunker->outputPoa = outputPoaFile;
    outputChunker->outputReadPartition = outputReadPartitionFile;
    outputChunker->outputRepeatCounts = outputRepeatCountFile;
    outputChunker->params = params;

    // Open files for writing
    outputChunker_open(outputChunker, "w");

    return outputChunker;
}

void
outputChunker_processChunkSequence(OutputChunker *outputChunker, int64_t chunkOrdinal, char *sequenceName, Poa *poa,
                                   stList *reads) {
    // Create chunk name
    char *headerLinePrefix = stString_print("%s,%" PRIi64 ",", sequenceName, chunkOrdinal);

    // Sequence
    if (outputChunker->outputSequenceFileHandle != NULL) {
        // Do run-length decoding
        char *outputSequence = rleString_expand(poa->refString);

        // Output the sequence, putting the sequence all on one line
        fprintf(outputChunker->outputSequenceFileHandle, "%s1\n%s\n", headerLinePrefix, outputSequence);

        free(outputSequence);
    }

    // Poa
    if (outputChunker->outputPoaFileHandle != NULL) {
        fprintf(outputChunker->outputPoaFileHandle, "%s%" PRIi64 "\n", headerLinePrefix, stList_length(poa->nodes) + 1);
        poa_printCSV(poa, outputChunker->outputPoaFileHandle, reads,
                     outputChunker->params->polishParams->repeatSubMatrix, 5);
    }

    // Now repeat counts
    if (outputChunker->outputRepeatCountFileHandle != NULL) {
        fprintf(outputChunker->outputRepeatCountFileHandle, "%s%" PRIi64 "\n", headerLinePrefix,
                stList_length(poa->nodes) + 1);
        poa_printRepeatCountsCSV(poa, outputChunker->outputRepeatCountFileHandle, reads);
    }

    // Cleanup
    free(headerLinePrefix);
}

void
outputChunker_processChunkSequencePhased2(OutputChunker *outputChunker, char *headerLinePrefix,
                                          Poa *poa, stList *reads, stSet *readsBelongingToHap1,
                                          stSet *readsBelongingToHap2) {
    // Output the sequence
    if (outputChunker->outputSequenceFileHandle != NULL) {
        char *outputSequence = rleString_expand(poa->refString);  // Do run-length decoding
        fprintf(outputChunker->outputSequenceFileHandle, "%s1\n%s\n", headerLinePrefix, outputSequence);
        free(outputSequence); // Cleanup
    }

    // Poa
    if (outputChunker->outputPoaFileHandle != NULL) {
        fprintf(outputChunker->outputPoaFileHandle, "%s%" PRIi64 "\n", headerLinePrefix, stList_length(poa->nodes) + 1);
        poa_printPhasedCSV(poa, outputChunker->outputPoaFileHandle, reads, readsBelongingToHap1, readsBelongingToHap2,
                           outputChunker->params->polishParams->repeatSubMatrix, 5);
    }

    // Output the repeat counts
    if (outputChunker->outputRepeatCountFileHandle != NULL) {
        fprintf(outputChunker->outputRepeatCountFileHandle, "%s%" PRIi64 "\n", headerLinePrefix,
                stList_length(poa->nodes) + 1);
        poa_printRepeatCountsCSV(poa, outputChunker->outputRepeatCountFileHandle, reads);
    }
}

void outputChunker_processChunkSequencePhased(OutputChunker *outputChunker, int64_t chunkOrdinal, char *sequenceName,
                                         Poa *poaHap1, Poa *poaHap2, stList *reads, stSet *readsBelongingToHap1,
                                         stSet *readsBelongingToHap2, stGenomeFragment *gF, Params *params) {
    // Create chunk name
    char *headerLinePrefix = stString_print("%s,%" PRIi64 ",", sequenceName, chunkOrdinal);

    outputChunker_processChunkSequencePhased2(outputChunker, headerLinePrefix,
                                              poaHap1, reads, readsBelongingToHap1, readsBelongingToHap2);

    outputChunker_processChunkSequencePhased2(outputChunker, headerLinePrefix,
                                              poaHap2, reads, readsBelongingToHap2, readsBelongingToHap1);

    // becasue we (may) have filtered reads now: readsBelongingToHapX has more reads than GF
    stSet *readIdsInGfHap1 = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    stSet *readIdsInGfHap2 = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    BamChunkRead *read = NULL;
    stSetIterator *itor = NULL;

    // Output the read partition hap1
    fprintf(outputChunker->outputReadPartitionFileHandle, "%s%" PRIi64 "\n", headerLinePrefix,
            stSet_size(readsBelongingToHap1) + 1);
    stGenomeFragment_printPartitionAsCSV(gF, outputChunker->outputReadPartitionFileHandle, params->phaseParams, 1,
                                         readIdsInGfHap1);
    itor = stSet_getIterator(readsBelongingToHap1);
    while ((read = stSet_getNext(itor)) != NULL) {
        if (stSet_search(readIdsInGfHap1, read->readName) == NULL) {
            fprintf(outputChunker->outputReadPartitionFileHandle, "%s,%f\n", read->readName, -1.0);
        }
    }
    stSet_destructIterator(itor);

    // Output the read partition hap2
    fprintf(outputChunker->outputReadPartitionFileHandle, "%s%" PRIi64 "\n", headerLinePrefix,
            stSet_size(readsBelongingToHap2) + 1);
    stGenomeFragment_printPartitionAsCSV(gF, outputChunker->outputReadPartitionFileHandle, params->phaseParams, 0,
                                         readIdsInGfHap2);
    itor = stSet_getIterator(readsBelongingToHap2);
    while ((read = stSet_getNext(itor)) != NULL) {
        if (stSet_search(readIdsInGfHap2, read->readName) == NULL) {
            fprintf(outputChunker->outputReadPartitionFileHandle, "%s,%f\n", read->readName, -1.0);
        }
    }
    stSet_destructIterator(itor);
    fflush(outputChunker->outputReadPartitionFileHandle);

    // Cleanup
    free(headerLinePrefix);
    stSet_destruct(readIdsInGfHap1);
    stSet_destruct(readIdsInGfHap2);
}

ChunkToStitch *outputChunker_readChunk(OutputChunker *outputChunker, bool phased) {
    /*
     * Read a chunk of output from the outputChunker.
     */
    ChunkToStitch *chunk = st_calloc(1, sizeof(ChunkToStitch));

    if (outputChunker->outputSequenceFile != NULL) {
        // primary "end of chunk" determinator
        if (!chunkToStitch_readSequenceChunk(outputChunker->outputSequenceFileHandle, chunk, phased)) {
            free(chunk);
            return NULL;
        }
    }
    if (phased) {
        // secondary "end of chunk" determinator.  used when we are not outputting sequence (for the very fast)
        if (!chunkToStitch_readReadPhasingChunk(outputChunker->outputReadPartitionFileHandle, chunk)) {
            if (outputChunker->outputSequenceFile != NULL) {
                st_errAbort("Expected chunk phasing info but found none! Expected chunk %"PRId64, chunk->chunkOrdinal);
            }
            free(chunk);
            return NULL;
        }
    }
    if (outputChunker->outputPoaFile != NULL) {
        chunkToStitch_readPoaChunk(outputChunker->outputPoaFileHandle, chunk, phased);
    }
    if (outputChunker->outputRepeatCountFile != NULL) {
        chunkToStitch_readRepeatCountChunk(outputChunker->outputRepeatCountFileHandle, chunk, phased);
    }

    return chunk;
}

void outputChunker_close(OutputChunker *outputChunker) {
    // Cleanup the sequence output file
    if (outputChunker->outputSequenceFileHandle != NULL) {
        fclose(outputChunker->outputSequenceFileHandle);
        outputChunker->outputSequenceFileHandle = NULL;
    }

    // Cleanup repeat count file
    if (outputChunker->outputRepeatCountFileHandle != NULL) {
        fclose(outputChunker->outputRepeatCountFileHandle);
        outputChunker->outputRepeatCountFileHandle = NULL;
    }

    // Cleanup poa file
    if (outputChunker->outputPoaFileHandle != NULL) {
        fclose(outputChunker->outputPoaFileHandle);
        outputChunker->outputPoaFileHandle = NULL;
    }

    // Cleanup read partition file
    if (outputChunker->outputReadPartitionFileHandle != NULL) {
        fclose(outputChunker->outputReadPartitionFileHandle);
        outputChunker->outputReadPartitionFileHandle = NULL;
    }
}

void outputChunker_closeAndDeleteFiles(OutputChunker *outputChunker) {
    /*
     * Closes the file streams and removes the output files (used for
     * chunker output to temporary files)
     */
    outputChunker_close(outputChunker); // Closes file streams

    if (!outputChunker->useMemoryBuffers) { // If not in memory need to delete underlying files
        // if in memory, buffers will be freed in destructor

        // Delete the sequence output file
        if (outputChunker->outputSequenceFile != NULL) {
            stFile_rmrf(outputChunker->outputSequenceFile);
        }

        // Delete repeat count file
        if (outputChunker->outputRepeatCountFile != NULL) {
            stFile_rmrf(outputChunker->outputRepeatCountFile);
        }

        // Delete the poa file
        if (outputChunker->outputPoaFile != NULL) {
            stFile_rmrf(outputChunker->outputPoaFile);
        }

        // Delete read partition file
        if (outputChunker->outputReadPartitionFile != NULL) {
            stFile_rmrf(outputChunker->outputReadPartitionFile);
        }
    }
}

void writeLines(FILE *fh, stList *lines) {
    for (int64_t i = 0; i < stList_length(lines); i++) {
        fprintf(fh, "%s\n", (char*)stList_get(lines, i));
    }
}

void outputChunker_writeChunkToFinalOutput(OutputChunker *outputChunker,
                                           char *seqName, char *seq, stList *poaLines, stList *repeatCountLines,
                                           bool startOfSequence) {
    /*
     * Writes the chunk to the final output files
     */

    // Write the sequence
    if (outputChunker->outputSequenceFile != NULL) {
        if (startOfSequence) {
            fastaWrite(seq, seqName, outputChunker->outputSequenceFileHandle);
        } else {
            fprintf(outputChunker->outputSequenceFileHandle, "%s\n", seq);
        }
    }

    // Write the POA
    if (outputChunker->outputPoaFile != NULL) {
        writeLines(outputChunker->outputPoaFileHandle, poaLines);
    }

    // Write the repeat counts
    if (outputChunker->outputRepeatCountFile != NULL) {
        writeLines(outputChunker->outputRepeatCountFileHandle, repeatCountLines);
    }
}

void outputChunker_destruct(OutputChunker *outputChunker) {
    // Close any open file descriptors
    outputChunker_close(outputChunker);

    // Cleanup the sequence output file
    if (outputChunker->outputSequenceFile != NULL) {
        free(outputChunker->outputSequenceFile);
    }

    // Cleanup repeat count file
    if (outputChunker->outputRepeatCountFile != NULL) {
        free(outputChunker->outputRepeatCountFile);
    }

    // Cleanup poa file
    if (outputChunker->outputPoaFile != NULL) {
        free(outputChunker->outputPoaFile);
    }

    // Cleanup read partition file
    if (outputChunker->outputReadPartitionFile != NULL) {
        free(outputChunker->outputReadPartitionFile);
    }

    // Cleanup residual
    free(outputChunker);
}

/*
 * OutputChunkers
 */

struct _outputChunkers {
    int64_t noOfOutputChunkers;
    stList *tempFileChunkers;
    int64_t tempFileChunkerCounter;
    int64_t chunkOrderNo;
    OutputChunker *outputChunkerHap1;
    OutputChunker *outputChunkerHap2;
    Params *params;
};

static char *printTempFileName(char *fileName, int64_t index) {
    if (fileName == NULL) {
        return NULL;
    }
    return stString_print("%s.%" PRIi64 ".temp", fileName, index);
}

static char *printFinalFileName(char *fileName, char *suffix) {
    if (fileName == NULL) {
        return NULL;
    }
    return stString_print("%s%s", fileName, suffix == NULL ? "" : suffix);
}

OutputChunkers *
outputChunkers_construct(int64_t noOfOutputChunkers, Params *params, char *outputSequenceFile, char *outputPoaFile,
                         char *outputReadPartitionFile, char *outputRepeatCountFile, char *hap1Suffix, char *hap2Suffix,
                         bool inMemoryBuffers) {

    char *outputReadPartitionFileForStitching = NULL;
    if (hap2Suffix != NULL) {
        if (stString_eq(hap1Suffix, hap2Suffix)) {
            st_errAbort("Hap1 and hap2 suffixes are identical, can not open distinct files for output\n");
        }
        // Make temporary read phasing file if not specified
        //TODO if inMemoryBuffers and !shouldOutputReadPartition, this results in shouldOutputReadPartition
        if (outputReadPartitionFile == NULL) {
            outputReadPartitionFileForStitching = "temp_read_phasing_file.csv";
            st_logInfo("> Making a temporary file to store read phasing in: %s\n", outputReadPartitionFileForStitching);
        } else {
            outputReadPartitionFileForStitching = outputReadPartitionFile;
        }
    } else {
        if (outputReadPartitionFile != NULL) {
            st_errAbort("Hap2 not specified but trying to output read partition\n");
        }
    }
    if (inMemoryBuffers) {
        st_logInfo("> Saving temporary data to in-memory buffers.\n");
    }

    OutputChunkers *outputChunkers = st_calloc(1, sizeof(OutputChunkers));
    outputChunkers->noOfOutputChunkers = noOfOutputChunkers;
    outputChunkers->params = params;
    // Make the temporary, parallel chunkers
    outputChunkers->tempFileChunkers = stList_construct3(0, (void (*)(void *)) outputChunker_destruct);
    for (int64_t i = 0; i < noOfOutputChunkers; i++) {
        stList_append(outputChunkers->tempFileChunkers, inMemoryBuffers ?
            outputChunker_constructInMemory(params, outputSequenceFile != NULL, outputPoaFile != NULL,
                    outputReadPartitionFileForStitching != NULL, outputRepeatCountFile != NULL) :
            outputChunker_construct(params, printTempFileName(outputSequenceFile, i), printTempFileName(outputPoaFile, i),
                printTempFileName(outputReadPartitionFileForStitching, i), printTempFileName(outputRepeatCountFile, i)));
    }
    // Make the final output chunkers
    outputChunkers->outputChunkerHap1 = outputChunker_construct(params,
                                                                printFinalFileName(outputSequenceFile, hap1Suffix),
                                                                printFinalFileName(outputPoaFile, hap1Suffix),
                                                                printFinalFileName(outputReadPartitionFile, hap1Suffix),
                                                                printFinalFileName(outputRepeatCountFile, hap1Suffix));

    if (hap2Suffix != NULL) {
        outputChunkers->outputChunkerHap2 = outputChunker_construct(params,
                                                                    printFinalFileName(outputSequenceFile, hap2Suffix),
                                                                    printFinalFileName(outputPoaFile, hap2Suffix),
                                                                    printFinalFileName(outputReadPartitionFile,
                                                                                       hap2Suffix),
                                                                    printFinalFileName(outputRepeatCountFile,
                                                                                       hap2Suffix));
    }

    return outputChunkers;
}

void outputChunkers_processChunkSequence(OutputChunkers *outputChunkers, int64_t chunker, int64_t chunkOrdinal,
                                    char *sequenceName, Poa *poa,
                                    stList *reads) {
    outputChunker_processChunkSequence(stList_get(outputChunkers->tempFileChunkers, chunker), chunkOrdinal,
                                       sequenceName, poa,
                                       reads);
}

void outputChunkers_processChunkSequencePhased(OutputChunkers *outputChunkers, int64_t chunker, int64_t chunkOrdinal,
                                               char *sequenceName, Poa *poaHap1, Poa *poaHap2, stList *reads,
                                               stSet *readsBelongingToHap1, stSet *readsBelongingToHap2,
                                               stGenomeFragment *gF, Params *params) {
    outputChunker_processChunkSequencePhased(stList_get(outputChunkers->tempFileChunkers, chunker), chunkOrdinal,
                                             sequenceName,
                                             poaHap1,
                                             poaHap2, reads, readsBelongingToHap1,
                                             readsBelongingToHap2, gF, params);
}

void outputChunkers_close(OutputChunkers *outputChunkers) {
    /*
     * Closes the file handles used.
     */
    for (int64_t i = 0; i < stList_length(outputChunkers->tempFileChunkers); i++) {
        outputChunker_close(stList_get(outputChunkers->tempFileChunkers, i));
    }
    outputChunker_close(outputChunkers->outputChunkerHap1);
    if (outputChunkers->outputChunkerHap2 != NULL) {
        outputChunker_close(outputChunkers->outputChunkerHap2);
    }
}

void outputChunkers_openForStitching(OutputChunkers *outputChunkers) {
    /*
     * Open the output chunkers for stitching
     */
    outputChunkers_close(outputChunkers);
    for (int64_t i = 0; i < stList_length(outputChunkers->tempFileChunkers); i++) {
        outputChunker_open(stList_get(outputChunkers->tempFileChunkers, i), "r");
    }
    outputChunker_open(outputChunkers->outputChunkerHap1, "w");
    if (outputChunkers->outputChunkerHap2 != NULL) {
        outputChunker_open(outputChunkers->outputChunkerHap2, "w");
    }
}

static ChunkToStitch *outputChunkers_readChunk(OutputChunkers *outputChunkers, bool phased) {
    for (int64_t i = 0; i < stList_length(outputChunkers->tempFileChunkers); i++) {
        ChunkToStitch *chunk = outputChunker_readChunk(stList_get(outputChunkers->tempFileChunkers,
                                                                  outputChunkers->tempFileChunkerCounter++ %
                                                                  stList_length(outputChunkers->tempFileChunkers)),
                                                       phased);
        if (chunk != NULL) {
            return chunk;
        }
    }
    return NULL;
}

static ChunkToStitch *outputChunkers_getNextChunkInSequence(OutputChunkers *outputChunkers,
                                                            stSortedSet *orderedChunks, bool phased) {
    ChunkToStitch *chunk;
    while (1) {
        chunk = outputChunkers_readChunk(outputChunkers, phased);
        if (chunk == NULL) {
            if (stSortedSet_size(orderedChunks) == 0) {
                return NULL;
            }
            chunk = stSortedSet_getFirst(orderedChunks);
            if (chunk->chunkOrdinal != outputChunkers->chunkOrderNo) {
                st_errAbort("Did not retrieve all the chunks from the temporary output");
            }
            break;
        }
        stSortedSet_insert(orderedChunks, chunk);
        chunk = stSortedSet_getFirst(orderedChunks);
        if (chunk->chunkOrdinal == outputChunkers->chunkOrderNo) {
            break;
        }
    }
    stSortedSet_remove(orderedChunks, chunk);
    outputChunkers->chunkOrderNo++;
    return chunk;
}

void outputChunkers_writeChunk(OutputChunkers *outputChunkers, ChunkToStitch *chunk) {
    /*
     * Writes the chunk to the final output files
     */
    outputChunker_writeChunkToFinalOutput(outputChunkers->outputChunkerHap1,
                                          chunk->seqName, chunk->seqHap1, chunk->poaHap1StringsLines,
                                          chunk->repeatCountLinesHap1, chunk->startOfSequence);
    if (outputChunkers->outputChunkerHap2 != NULL) {
        outputChunker_writeChunkToFinalOutput(outputChunkers->outputChunkerHap2,
                                              chunk->seqName, chunk->seqHap2, chunk->poaHap2StringsLines,
                                              chunk->repeatCountLinesHap2,
                                              chunk->startOfSequence);
    }
}


void writeReadPartition(stHash *readsInHap, FILE *fh) {
    /*
     * Write out the reads for a haplotype in the given file
     */
    fprintf(fh, "READ_NAME,PHRED_SCORE_OF_BEING_IN_PARTITION\n");
    stHashIterator *it = stHash_getIterator(readsInHap);
    char *readName;
    while ((readName = stHash_getNext(it)) != NULL) {
        double *prob = stHash_search(readsInHap, readName);
        fprintf(fh, "%s,%f\n", readName, *prob);
    }
    stHash_destructIterator(it);
}

void outputChunkers_stitchLinear(OutputChunkers *outputChunkers, bool phased) {
    /*
     * Stitch together the outputs using a single thread, but very minimal memory.
     */

    // Setup to write out the chunks
    outputChunkers_openForStitching(outputChunkers);

    // Create a cache to hold the chunks, ordered by their ordinal
    stSortedSet *orderedChunks = stSortedSet_construct3((int (*)(const void *, const void *)) chunkToStitch_cmp, NULL);

    // Get the first chunk
    ChunkToStitch *pChunk = outputChunkers_getNextChunkInSequence(outputChunkers, orderedChunks, phased), *chunk;

    if (pChunk == NULL) {
        // Nothing to do
        stSortedSet_destruct(orderedChunks);
        return;
    }

    // Length of the output sequences
    int64_t lengthOfSequenceOutputSoFarHap1 = 0;
    int64_t lengthOfSequenceOutputSoFarHap2 = 0;

    // Track the names of the reads in the two haplotypes, if phased
    stHash *hap1Reads = NULL, *hap2Reads = NULL;
    if (phased) {
        hap1Reads = getReadNames(pChunk->readsHap1Lines);
        hap2Reads = getReadNames(pChunk->readsHap2Lines);
    }

    // Indicate we're at the of beginning a sequences
    pChunk->startOfSequence = 1;

    // Get each successive chunk and stitch and phase progressively
    while ((chunk = outputChunkers_getNextChunkInSequence(outputChunkers, orderedChunks, phased)) != NULL) {
        assert(pChunk != NULL);

        // If phased, ensure the chunks phasing is consistent
        if (phased) {
            chunkToStitch_phaseAdjacentChunks(chunk, hap1Reads, hap2Reads);
        }

        // Set the flag determining if this is the start of a new sequence
        chunk->startOfSequence = !stString_eq(pChunk->seqName, chunk->seqName);

        if (chunk->startOfSequence) { // Reset the lengths of the new sequence output
            lengthOfSequenceOutputSoFarHap1 = 0;
            lengthOfSequenceOutputSoFarHap2 = 0;
        } else { // Trim the boundaries of the two chunks so that they don't overlap
            chunkToStitch_trimAdjacentChunks(pChunk, chunk, outputChunkers->params,
                                             &lengthOfSequenceOutputSoFarHap1, &lengthOfSequenceOutputSoFarHap2);
        }

        // Write out the pChunk
        outputChunkers_writeChunk(outputChunkers, pChunk);

        // Cleanup the pChunk
        chunkToStitch_destruct(pChunk);

        // Set the new pChunk
        pChunk = chunk;
    }

    // Write out the pChunk
    outputChunkers_writeChunk(outputChunkers, pChunk);

    // Cleanup the pChunk
    chunkToStitch_destruct(pChunk);

    // Write out the read name phasing, if needed
    if (phased) {
        writeReadPartition(hap1Reads, outputChunkers->outputChunkerHap1->outputReadPartitionFileHandle);
        writeReadPartition(hap2Reads, outputChunkers->outputChunkerHap2->outputReadPartitionFileHandle);
    }

    // Cleanup
    if (stSortedSet_size(orderedChunks) != 0) {
        st_errAbort("Got chunks left over after writing out chunks");
    }
    stSortedSet_destruct(orderedChunks);
    if (phased) {
        stHash_destruct(hap1Reads);
        stHash_destruct(hap2Reads);
    }
}

void updateStitchingChunk(ChunkToStitch *stitched, ChunkToStitch *pChunk, stList *hap1Seqs, stList *hap2Seqs,
                          bool phased, bool trackPoa, bool trackRepeatCounts) {

    // for very fast case where we don't write fasta out
    if (hap1Seqs == NULL && hap2Seqs == NULL) {
        return;
    }

    stList_append(hap1Seqs, stString_copy(pChunk->seqHap1));
    if (phased) {
        stList_append(hap2Seqs, stString_copy(pChunk->seqHap2));
    }
    if (trackPoa) {
        stList_appendAll(stitched->poaHap1StringsLines, pChunk->poaHap1StringsLines);
        stList_setDestructor(pChunk->poaHap1StringsLines, NULL);
        if (phased) {
            stList_appendAll(stitched->poaHap2StringsLines, pChunk->poaHap2StringsLines);
            stList_setDestructor(pChunk->poaHap2StringsLines, NULL);
        }
    }
    if (trackRepeatCounts) {
        stList_appendAll(stitched->repeatCountLinesHap1, pChunk->repeatCountLinesHap1);
        stList_setDestructor(pChunk->repeatCountLinesHap1, NULL);
        if (phased) {
            stList_appendAll(stitched->repeatCountLinesHap2, pChunk->repeatCountLinesHap2);
            stList_setDestructor(pChunk->repeatCountLinesHap2, NULL);
        }
    }
}

void convertReadPartitionToLines(stHash *readsInHap, stList *readPartitionLines) {
    /*
     * Format the output of the reads for a haplotype
     */
    stList_append(readPartitionLines, stString_print("READ_NAME,LOG_PROB_OF_BEING_IN_PARTITION\n"));
    stHashIterator *it = stHash_getIterator(readsInHap);
    char *readName;
    while ((readName = stHash_getNext(it)) != NULL) {
        double *prob = stHash_search(readsInHap, readName);
        stList_append(readPartitionLines, stString_print("%s,%f\n", readName, *prob));
    }
    stHash_destructIterator(it);
}

ChunkToStitch *mergeContigChunkz(ChunkToStitch **chunks, int64_t startIdx, int64_t endIdxExclusive, bool phased,
        Params *params) {
    // sanity check

    // Get the first chunk
    ChunkToStitch *pChunk = chunks[startIdx];
    ChunkToStitch *chunk = NULL;

    // Length of the output sequences
    int64_t lengthOfSequenceOutputSoFarHap1 = 0;
    int64_t lengthOfSequenceOutputSoFarHap2 = 0;

    // Our "stitched" chunk
    bool trackSequence = pChunk->seqHap1 != NULL;
    bool trackRepeatCounts = pChunk->repeatCountLinesHap1 != NULL;
    bool trackPoa = pChunk->poaHap1StringsLines != NULL;
    ChunkToStitch *stitched = chunkToStitch_construct(NULL, -1, phased, trackRepeatCounts, trackPoa);

    // Lists to keep track of haplotype strings
    stList *hap1Seqs = (trackSequence ? stList_construct3(0, free) : NULL);
    stList *hap2Seqs = (phased && trackSequence ? stList_construct3(0, free) : NULL);

    // Track the names of the reads in the two haplotypes, if phased
    stHash *hap1Reads, *hap2Reads;
    if (phased) {
        hap1Reads = getReadNames(pChunk->readsHap1Lines);
        hap2Reads = getReadNames(pChunk->readsHap2Lines);
    }

    // Get each successive chunk and stitch and phase progressively
    for (int64_t chunkIndex = startIdx + 1; chunkIndex < endIdxExclusive; chunkIndex++) {
        assert(pChunk != NULL);
        chunk = chunks[chunkIndex];

        // If phased, ensure the chunks phasing is consistent
        if (phased) {
            chunkToStitch_phaseAdjacentChunks(chunk, hap1Reads, hap2Reads);
        }

        // Trim the overlap between chunks
        chunkToStitch_trimAdjacentChunks(pChunk, chunk, params,
                                         &lengthOfSequenceOutputSoFarHap1, &lengthOfSequenceOutputSoFarHap2);

        // Save to stitched
        updateStitchingChunk(stitched, pChunk, hap1Seqs, hap2Seqs, phased, trackPoa, trackRepeatCounts);

        // cleanup
        chunkToStitch_destruct(pChunk);

        // Set the new pChunk
        pChunk = chunk;
    }

    // save the last chunk and the read phasing
    updateStitchingChunk(stitched, pChunk, hap1Seqs, hap2Seqs, phased, trackPoa, trackRepeatCounts);
    stitched->seqHap1 = trackSequence ? stString_join2("", hap1Seqs) : NULL;
    if (phased) {
        stitched->seqHap2 = trackSequence ? stString_join2("", hap2Seqs) : NULL;
        // Save back the reads in each haplotype to the stitched chunk
        convertReadPartitionToLines(hap1Reads, stitched->readsHap1Lines);
        convertReadPartitionToLines(hap2Reads, stitched->readsHap2Lines);
    }

    // cleanup
    chunkToStitch_destruct(pChunk);
    if (trackSequence) stList_destruct(hap1Seqs);
    if (phased) {
        stList_destruct(hap2Seqs);
        stHash_destruct(hap1Reads);
        stHash_destruct(hap2Reads);
    }

    // fin
    return stitched;
}

ChunkToStitch *mergeContigChunkzThreaded(ChunkToStitch **chunks, int64_t startIdx, int64_t endIdxExclusive, int64_t numThreads,
                                bool phased, Params *params, char *referenceSequenceName) {

    // special unthreaded case
    if (numThreads == 1) return mergeContigChunkz(chunks, startIdx, endIdxExclusive, phased, params);

    // divide into chunks
    int64_t totalChunks = endIdxExclusive - startIdx;
    int64_t chunksPerThread = (int64_t) ceil(1.0 * totalChunks / numThreads);
    while (startIdx + chunksPerThread * (numThreads - 1) >= endIdxExclusive) {numThreads--;}
    ChunkToStitch **outputChunks = st_calloc(numThreads, sizeof(ChunkToStitch*));

    // multithread loop
    st_logInfo("  Merging chunks for %s from (%"PRId64", %"PRId64"] with %"PRId64" chunks per thread on %"PRId64" threads \n",
               referenceSequenceName, startIdx, endIdxExclusive, chunksPerThread, numThreads);
    #pragma omp parallel for schedule(static,1)
    for (int64_t thread = 0; thread < numThreads; thread++) {
        int64_t threadedStartIdx = startIdx + chunksPerThread * thread;
        int64_t threadedEndIdxExclusive = threadedStartIdx + chunksPerThread;
        if (endIdxExclusive < threadedEndIdxExclusive) threadedEndIdxExclusive = endIdxExclusive;

        outputChunks[thread] = mergeContigChunkz(chunks, threadedStartIdx, threadedEndIdxExclusive, phased, params);
    }

    // finish
    ChunkToStitch *stitched = mergeContigChunkz(outputChunks, 0, numThreads, phased, params);
    free(outputChunks); //chunks are freed in mergeCC
    return stitched;
}


void outputChunkers_stitch(OutputChunkers *outputChunkers, bool phased, int64_t chunkCount) {
    outputChunkers_stitchAndTrackReadIds(outputChunkers, phased, chunkCount, NULL, NULL);
}
void outputChunkers_stitchAndTrackReadIds(OutputChunkers *outputChunkers, bool phased, int64_t chunkCount,
        stList *readIdsHap1, stList *readIdsHap2) {

    // prep for merge
    assert(chunkCount > 0);

    // Setup to write out the chunks
    outputChunkers_openForStitching(outputChunkers);

    // Create a cache to hold the chunks, ordered by their ordinal
    ChunkToStitch **chunks = st_calloc(chunkCount, sizeof(ChunkToStitch *));
    int64_t *foundChunksPerThread = st_calloc(outputChunkers->noOfOutputChunkers, sizeof(int64_t));

    /// get all chunks
    # ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1)
    # endif
    for (int64_t i = 0; i < outputChunkers->noOfOutputChunkers; i++) {
        foundChunksPerThread[i] = 0;
        ChunkToStitch *chunk = NULL;
        while ((chunk = outputChunker_readChunk(stList_get(outputChunkers->tempFileChunkers, i), phased)) != NULL) {
            if (chunks[chunk->chunkOrdinal] != NULL) {
                st_errAbort("Encountered chunk %"PRId64" twice while reading from temp files!");
            }
            chunks[chunk->chunkOrdinal] = chunk;
            foundChunksPerThread[i]++;
        }
    }

    // sanity check debugging
    int64_t foundChunks = 0;
    for (int64_t i = 0; i < outputChunkers->noOfOutputChunkers; i++) {
        foundChunks += foundChunksPerThread[i];
    }
    free(foundChunksPerThread);
    if (foundChunks != chunkCount) {
        int64_t i = 0;
        stList *missingChunks = stList_construct3(0, free);
        while (i < outputChunkers->noOfOutputChunkers && stList_length(missingChunks) < 10) {
            if (chunks[i] == NULL) stList_append(missingChunks, stString_print("%s", i));
        }
        if (stList_length(missingChunks) == 10 && i != outputChunkers->noOfOutputChunkers) stList_append(missingChunks, "..");
        st_errAbort("Missing %"PRId64" chunks: %s\n", chunkCount - foundChunks, stString_join2(", ", missingChunks));
    }

    // prep for merging
    int64_t contigStartIdx = 0;
    char *referenceSequenceName = stString_copy(chunks[0]->seqName);
    int64_t  lastReportedPercentage = 0;
    time_t mergeStartTime = time(NULL);
    st_logCritical("> Merging results from %"PRIu64" chunks.\n", chunkCount);

    // find which chunks belong to each contig, merge each contig threaded, write out
    for (int64_t chunkIdx = 1; chunkIdx <= chunkCount; chunkIdx++) {

        // we encountered the last chunk in the contig (end of list or new refSeqName)
        if (chunkIdx == chunkCount || !stString_eq(referenceSequenceName, chunks[chunkIdx]->seqName)) {

            // generate and save sequences
            ChunkToStitch *stitched = mergeContigChunkzThreaded(chunks, contigStartIdx, chunkIdx,
                    outputChunkers->noOfOutputChunkers, phased, outputChunkers->params, referenceSequenceName);
            stitched->seqName = stString_copy(referenceSequenceName);
            stitched->startOfSequence = true;
            outputChunkers_writeChunk(outputChunkers, stitched);

            // to write to bam, we need to add all these
            if (readIdsHap1 != NULL && readIdsHap2 != NULL) {
                stHash *chunkReadToProbHap1 = getReadNames(stitched->readsHap1Lines);
                stHash *chunkReadToProbHap2 = getReadNames(stitched->readsHap2Lines);
                stList *chunkReadsHap1 = stHash_getKeys(chunkReadToProbHap1);
                stList *chunkReadsHap2 = stHash_getKeys(chunkReadToProbHap2);
                stList_appendAll(readIdsHap1, chunkReadsHap1);
                stList_appendAll(readIdsHap2, chunkReadsHap2);
                stHash_setDestructKeys(chunkReadToProbHap1, NULL);
                stHash_setDestructKeys(chunkReadToProbHap2, NULL);
                stList_setDestructor(chunkReadsHap1, NULL);
                stList_setDestructor(chunkReadsHap2, NULL);
                stHash_destruct(chunkReadToProbHap1);
                stHash_destruct(chunkReadToProbHap2);
                stList_destruct(chunkReadsHap1);
                stList_destruct(chunkReadsHap2);
            }

            // log progress
            int64_t currentPercentage = (int64_t) (100 * chunkIdx / chunkCount);
            if (currentPercentage != lastReportedPercentage) {
                lastReportedPercentage = currentPercentage;
                int64_t timeTaken = (int64_t) (time(NULL) - mergeStartTime);
                int64_t secondsRemaining = (int64_t) floor(
                        1.0 * timeTaken / currentPercentage * (double) (100 - currentPercentage));
                char *timeDescriptor = (secondsRemaining == 0 && currentPercentage <= 50 ?
                                        stString_print("unknown") : getTimeDescriptorFromSeconds(secondsRemaining));
                st_logCritical("> Merging %2"PRId64"%% complete (%"PRId64"/%"PRId64").  Estimated time remaining: %s\n",
                               currentPercentage, chunkIdx, chunkCount, timeDescriptor);
                free(timeDescriptor);
            }

            // Clean up
            chunkToStitch_destruct(stitched);
            free(referenceSequenceName);

            // Reset for next reference sequence
            if (chunkIdx != chunkCount) {
                contigStartIdx = chunkIdx;
                referenceSequenceName = stString_copy(chunks[chunkIdx]->seqName);
            }
        }
        // nothing to do otherwise, just wait until end or new contig
    }

    // cleanup
    free(chunks); //chunks are freed as they're merged

}


void outputChunkers_destruct(OutputChunkers *outputChunkers) {
    // Close the file streams and delete the temporary files of the temp file chunkers
    for (int64_t i = 0; i < stList_length(outputChunkers->tempFileChunkers); i++) {
        outputChunker_closeAndDeleteFiles(stList_get(outputChunkers->tempFileChunkers, i));
    }
    time_t start = time(NULL);
    // Now cleanup the temp file chunkers
    stList_destruct(outputChunkers->tempFileChunkers);
    // Cleanup the final output chunkers
    outputChunker_destruct(outputChunkers->outputChunkerHap1);
    if (outputChunkers->outputChunkerHap2 != NULL) {
        outputChunker_destruct(outputChunkers->outputChunkerHap2);
    }
    free(outputChunkers);
    char *timeDes = getTimeDescriptorFromSeconds(time(NULL) - start);
    st_logInfo("    Closed remaining output chunking infrastructure in %s\n", timeDes);
    free(timeDes);

}