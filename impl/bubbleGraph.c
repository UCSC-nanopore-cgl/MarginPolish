/*
 * Copyright (C) 2019 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"

/*
 * Bubble graphs
 */

int64_t bubble_getReferenceAlleleIndex(Bubble *b) {
    for (int64_t i = 0; i < b->alleleNo; i++) {
        if (rleString_eq(b->refAllele, b->alleles[i])) {
            return i;
        }
    }
    return -1;
}

double rleString_calcLogProb(RleString *allele, PolishParams *p) {
    double lProb = 0.0;
    for (int64_t i = 0; i < allele->length; i++) {
        lProb += log(0.25) + log(0.01) + 2.3025 * p->repeatSubMatrix->baseLogProbs_AT[allele->repeatCounts[i]];
    }
    return lProb;
}

double bubble_getLogLikelihoodOfAllele(Bubble *b, int64_t allele, PolishParams *p) {
    double logLikelihood = 0.0;
    for (int64_t i = 0; i < b->readNo; i++) {
        logLikelihood += b->alleleReadSupports[allele * b->readNo + i];
    }
    return logLikelihood; // + rleString_calcLogProb(b->alleles[allele], p);
}

int64_t bubble_getIndexOfHighestLikelihoodAllele(Bubble *b, PolishParams *p) {
    int64_t maxAllele = 0;
    assert(b->alleleNo > 0);
    double maxAlleleLikelihood = bubble_getLogLikelihoodOfAllele(b, 0, p);
    for (int64_t i = 1; i < b->alleleNo; i++) {
        double alleleLikelihood = bubble_getLogLikelihoodOfAllele(b, i, p);
        if (alleleLikelihood > maxAlleleLikelihood) {
            maxAllele = i;
            maxAlleleLikelihood = alleleLikelihood;
        }
    }
    return maxAllele;
}

uint64_t *bubbleGraph_getConsensusPath(BubbleGraph *bg, PolishParams *polishParams) {
    uint64_t *consensusPath = st_calloc(bg->bubbleNo, sizeof(uint64_t));
    for (int64_t i = 0; i < bg->bubbleNo; i++) {
        Bubble *b = &(bg->bubbles[i]);
        consensusPath[i] = bubble_getIndexOfHighestLikelihoodAllele(b, polishParams);
    }
    return consensusPath;
}

RleString *bubbleGraph_getConsensusString(BubbleGraph *bg, uint64_t *consensusPath,
                                          int64_t **poaToConsensusMap, PolishParams *polishParams) {
    // Map to track alignment between the new consensus sequence and the current reference sequence
    *poaToConsensusMap = st_malloc(bg->refString->length * sizeof(int64_t));
    for (int64_t i = 0; i < bg->refString->length; i++) {
        (*poaToConsensusMap)[i] = -1;
    }

    // Substrings of the consensus string that when concatenated form the overall consensus string
    stList *consensusSubstrings = stList_construct3(0, free);
    char previousBase = '-';
    int64_t j = 0; // Index in the consensus substring
    int64_t k = 0; // Index in the reference string
    int64_t totalDiffs = 0; // Index to keep track of number of alleles changed for debug printing
    for (int64_t i = 0; i < bg->bubbleNo; i++) {
        Bubble *b = &(bg->bubbles[i]);

        // Add prefix after the last bubble (or start) but before the new bubble start
        if (k < b->refStart) {
            // Get substring
            RleString *refSubString = rleString_copySubstring(bg->refString, k, b->refStart - k);
            assert(refSubString->length > 0);
            stList_append(consensusSubstrings, rleString_expand(refSubString));

            // Update coordinate map between old and new reference

            // Skip an element in the consensus string if the same as the previous base
            // as will get squashed when run length encoded
            if (polishParams->useRunLengthEncoding && refSubString->rleString[0] == previousBase) {
                k++;
            }

            while (k < b->refStart) {
                (*poaToConsensusMap)[k++] = j++;
            }
            previousBase = refSubString->rleString[refSubString->length - 1];

            // Cleanup
            rleString_destruct(refSubString);
        }

        // Add the bubble string itself
        // Noting, if there are not sufficient numbers of sequences to call the consensus
        // use the current reference sequence
        RleString *consensusSubstring = b->alleles[consensusPath[i]];
        assert(consensusSubstring->length > 0);
        stList_append(consensusSubstrings, rleString_expand(consensusSubstring));

        if (st_getLogLevel() >= debug) {
            if (!rleString_eq(consensusSubstring, b->refAllele)) {
                st_logDebug(
                        "In bubbleGraph_getConsensus (diff %" PRIi64 " , from: %" PRIi64 " to: %" PRIi64 ", \nexisting string:\t",
                        totalDiffs++, k, k + b->refAllele->length);
                rleString_print(b->refAllele, stderr);
                st_logDebug("\nnew string:\t\t");
                rleString_print(consensusSubstring, stderr);
                st_logDebug("\n");

                for (int64_t l = 0; l < b->alleleNo; l++) {
                    st_logDebug("\tGot allele: \t");
                    rleString_print(b->alleles[l], stderr);
                    st_logDebug(" with log-likelihood: %f\n", bubble_getLogLikelihoodOfAllele(b, l, polishParams));
                }

                for (int64_t l = 0; l < b->readNo; l++) {
                    RleString *readSubstring = bamChunkReadSubstring_getRleString(b->reads[l]);
                    st_logDebug("\tGot read: \t");
                    rleString_print(readSubstring, stderr);
                    st_logDebug(", q-value: %f\n", b->reads[l]->qualValue);
                    rleString_destruct(readSubstring);
                }
            }
        }

        // Check if the same as the existing reference
        // in which case we can maintain the alignment
        if (rleString_eq(consensusSubstring, b->refAllele)) {
            if (polishParams->useRunLengthEncoding && consensusSubstring->rleString[0] == previousBase) {
                k++;
            }
            while (k < b->refStart + b->refAllele->length) {
                (*poaToConsensusMap)[k++] = j++;
            }
        } else {
            // Otherwise just update coordinates
            k += b->refAllele->length;
            j += consensusSubstring->length +
                 // Latter expression establishes if the first position will be compressed into the earlier one
                 (polishParams->useRunLengthEncoding && consensusSubstring->rleString[0] == previousBase ? -1 : 0);
        }
        previousBase = consensusSubstring->rleString[consensusSubstring->length - 1];
    }

    // Add the suffix of the reference after the last bubble
    if (k < bg->refString->length) {
        RleString *refSubString = rleString_copySubstring(bg->refString, k, bg->refString->length - k);
        stList_append(consensusSubstrings, rleString_expand(refSubString));

        if (polishParams->useRunLengthEncoding && refSubString->rleString[0] == previousBase) {
            k++;
        }
        while (k < bg->refString->length) {
            (*poaToConsensusMap)[k++] = j++;
        }

        rleString_destruct(refSubString);
    }

    // Build the new consensus string by concatenating the constituent pieces
    char *newExpandedConsensusString = stString_join2("", consensusSubstrings);
    RleString *newConsensusString = polishParams->useRunLengthEncoding ? rleString_construct(newExpandedConsensusString)
                                                                       : rleString_construct_no_rle(
                    newExpandedConsensusString);

    assert(newConsensusString->length == j);

    // Cleanup
    stList_destruct(consensusSubstrings);
    free(newExpandedConsensusString);

    return newConsensusString;
}

// New polish algorithm

double getTotalWeight(Poa *poa, PoaNode *node) {
    /*
     * Returns the total base weight of reads aligned to the given node.
     */
    double totalWeight = 0.0;
    for (int64_t i = 0; i < poa->alphabet->alphabetSize; i++) {
        totalWeight += node->baseWeights[i];
    }
    return totalWeight;
}

double getAvgCoverage(Poa *poa, int64_t from, int64_t to) {
    // Calculate average coverage, which is used to determine candidate variants
    double avgCoverage = 0.0;
    for (int64_t j = from; j < to; j++) {
        avgCoverage += getTotalWeight(poa, stList_get(poa->nodes, j));
    }
    return avgCoverage / (to - from);
}

char getNextCandidateBase(Poa *poa, PoaNode *node, int64_t *i, double candidateWeight) {
    /*
     * Iterates through candidate bases for a reference position returning those with sufficient weight.
     * Always returns the reference base
     */
    while (*i < poa->alphabet->alphabetSize) {
        char base = poa->alphabet->convertSymbolToChar(*i);
        if (node->baseWeights[(*i)++] > candidateWeight || toupper(node->base) == base) {
            return base;
        }
    }
    return '-';
}

int64_t getNextCandidateRepeatCount(Poa *poa, PoaNode *node, int64_t *i, double candidateWeight) {
    /*
     * Iterates through candidate repeat counts for a reference position returning those with sufficient weight.
     * Always returns the reference repeat count.
     */
    candidateWeight *= 2.0; // This is a hack to reduce the number of repeat counts investigated by making a repeat count need a larger change
    while (*i < poa->maxRepeatCount) {
        int64_t repeatCount = (*i)++;
        if (node->repeatCountWeights[repeatCount] > candidateWeight || node->repeatCount == repeatCount) {
            return repeatCount;
        }
    }
    return -1;
}

bool hasCandidateSubstitution(Poa *poa, PoaNode *node, double candidateWeight) {
    /*
     * Returns non-zero if the node has a candidate base that is different to the
     * current base.
     */
    int64_t i = 0;
    char base;
    while ((base = getNextCandidateBase(poa, node, &i, candidateWeight)) != '-') {
        if (base != node->base) {
            return 1;
        }
    }
    return 0;
}

bool hasCandidateRepeatCountChange(Poa *poa, PoaNode *node, double candidateWeight) {
    /*
     * Returns non-zero if the node has a candidate base repeat count that is different to the
     * current base's repeat count.
     */
    int64_t i = 0;
    int64_t repeatCount;
    while ((repeatCount = getNextCandidateRepeatCount(poa, node, &i, candidateWeight)) != -1) {
        if (repeatCount != node->repeatCount) {
            return 1;
        }
    }
    return 0;
}

RleString *getNextCandidateInsert(PoaNode *node, int64_t *i, double candidateWeight) {
    /*
     * Iterates through candidate insertions for a reference position returning those with sufficient weight.
     */
    while ((*i)++ < stList_length(node->inserts)) {
        PoaInsert *insert = stList_get(node->inserts, (*i) - 1);
        if (poaInsert_getWeight(insert) > candidateWeight) {
            return insert->insert;
        }
    }
    return NULL;
}

bool hasCandidateInsert(PoaNode *node, double candidateWeight) {
    /*
     * Returns non-zero if the node has a candidate insert.
     */
    int64_t i = 0;
    return getNextCandidateInsert(node, &i, candidateWeight) != NULL;
}

int64_t getNextCandidateDelete(PoaNode *node, int64_t *i, double candidateWeight) {
    /*
     * Iterates through candidate deletions for a reference position returning those with sufficient weight.
     */
    while ((*i)++ < stList_length(node->deletes)) {
        PoaDelete *delete = stList_get(node->deletes, (*i) - 1);
        if (poaDelete_getWeight(delete) > candidateWeight) {
            return delete->length;
        }
    }
    return -1;
}

bool maxCandidateDeleteLength(PoaNode *node, double candidateWeight) {
    /*
     * Returns maximum length of a candidate deletion starting after this position.
     */
    int64_t i = 0;
    int64_t deleteLength, maxDeleteLength = 0;
    while ((deleteLength = getNextCandidateDelete(node, &i, candidateWeight)) != -1) {
        if (deleteLength > maxDeleteLength) {
            maxDeleteLength = deleteLength;
        }
    }
    return maxDeleteLength;
}

static bool containsString(stList *strings, char *s) {
    for (int64_t i = 0; i < stList_length(strings); i++) {
        if (stString_eq(stList_get(strings, i), s)) {
            return 1;
        }
    }
    return 0;
}

stList *getCandidateConsensusSubstrings(Poa *poa, int64_t from, int64_t to,
                                        double *candidateWeights, double weightAdjustment,
                                        int64_t maximumStringNumber) {
    /*
     *  A candidate variant is an edit (either insert, delete or substitution) to the poa reference
     *  string with "high" weight. This function returns all possible combinations of candidate variants,
     *  each as a new consensus substring, for the interval of the reference string from "from" (inclusive)
     *  to "to" (exclusive). Returned list of strings always contains the reference string without edits (the no
     *  candidate variants string).
     */

    // Function is recursive

    // First get suffix substrings
    stList *suffixes;
    if (from + 1 < to) {
        suffixes = getCandidateConsensusSubstrings(poa, from + 1, to, candidateWeights, weightAdjustment,
                                                   maximumStringNumber);

        if (suffixes == NULL) { // If too many combinations, return null.
            return NULL;
        }
    } else {
        suffixes = stList_construct3(0, free);
        stList_append(suffixes, stString_copy("")); // Start with the empty string
    }

    // Now extend by adding on prefix variants
    stList *consensusSubstrings = stList_construct3(0, free);

    PoaNode *node = stList_get(poa->nodes, from);

    double candidateWeight = candidateWeights[from] * weightAdjustment;

    int64_t i = 0;
    char base;
    while ((base = getNextCandidateBase(poa, node, &i, candidateWeight)) !=
           '-') { // Enumerate the possible bases at the reference node.

        int64_t repeatCount, l = 1;
        while ((repeatCount = getNextCandidateRepeatCount(poa, node, &l, candidateWeight)) !=
               -1) { // Enumerate the possible repeat counts at the reference node.
            assert(repeatCount != 0);
            char *bases = expandChar(base, repeatCount);

            // Create the consensus substrings with no inserts or deletes starting at this node
            for (int64_t j = 0; j < stList_length(suffixes); j++) {
                stList_append(consensusSubstrings, stString_print("%s%s", bases, stList_get(suffixes, j)));
            }

            // Now add insert cases
            int64_t k = 0;
            RleString *insert;
            while ((insert = getNextCandidateInsert(node, &k, candidateWeight)) != NULL) {
                char *expandedInsert = rleString_expand(insert);
                assert(strlen(expandedInsert) > 0);
                for (int64_t j = 0; j < stList_length(suffixes); j++) {
                    stList_append(consensusSubstrings,
                                  stString_print("%s%s%s", bases, expandedInsert, stList_get(suffixes, j)));
                }
                free(expandedInsert);
            }

            // Add then deletes
            k = 0;
            int64_t deleteLength;
            while ((deleteLength = getNextCandidateDelete(node, &k, candidateWeight)) > 0) {
                for (int64_t j = 0; j < stList_length(suffixes); j++) {
                    char *suffixHaplotype = stList_get(suffixes, j);

                    // Build new deletion
                    char *s = stString_print("%s%s", bases,
                                             ((int64_t) strlen(suffixHaplotype) - deleteLength >= 0)
                                             ? &(suffixHaplotype[deleteLength]) : "");

                    // Add deletion if not already in the set of consensus strings
                    if (!containsString(consensusSubstrings, s)) {
                        stList_append(consensusSubstrings, s);
                    } else {
                        free(s);
                    }
                }
            }

            // Cleanup bases
            free(bases);
        }
    }

    // Clean up
    stList_destruct(suffixes);

    if (stList_length(consensusSubstrings) > maximumStringNumber) {
        // Clean up and return null (too many combinations)
        stList_destruct(consensusSubstrings);
        return NULL;
    }

    return consensusSubstrings;
}

BamChunkReadSubstring *
bamChunkRead_getSubstring(BamChunkRead *bamChunkRead, int64_t start, int64_t length, PolishParams *params) {
    assert(length >= 0);

    BamChunkReadSubstring *rs = st_calloc(1, sizeof(BamChunkReadSubstring));

    // Basic attributes
    rs->read = bamChunkRead;
    rs->start = start;
    rs->length = length;

    // Calculate the qual value
    if (bamChunkRead->qualities != NULL) {
        int64_t j = 0;
        for (int64_t i = 0; i < length; i++) {
            j += (int64_t) bamChunkRead->qualities[i + start];
        }
        rs->qualValue = (double) j / length; // Quals are phred, qual = -10 * log_10(p)
    } else {
        rs->qualValue = -1.0;
    }

    return rs;
}

RleString *bamChunkReadSubstring_getRleString(BamChunkReadSubstring *readSubstring) {
    return rleString_copySubstring(readSubstring->read->rleRead, readSubstring->start, readSubstring->length);
}

void bamChunkReadSubstring_destruct(BamChunkReadSubstring *rs) {
    free(rs);
}

int poaBaseObservation_cmp(const void *a, const void *b) {
    PoaBaseObservation *obs1 = (PoaBaseObservation *) a;
    PoaBaseObservation *obs2 = (PoaBaseObservation *) b;
    if (obs1->readNo != obs2->readNo) { // Sort first is ascending read number order
        return obs1->readNo < obs2->readNo ? -1 : 1;
    }
    if (obs1->weight != obs2->weight) { // Sort second in descending weight order
        return obs1->weight > obs2->weight ? -1 : 1;
    }
    return 0;
}

void sortBaseObservations(Poa *poa) {
    /*
     * Sort the POA base observations to make them appropriate for getReadSubstrings.
     */
    for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
        PoaNode *node = stList_get(poa->nodes, i);
        stList_sort(node->observations, poaBaseObservation_cmp);
    }
}

int64_t skipDupes(PoaNode *node, int64_t i, int64_t readNo) {
    while (i < stList_length(node->observations)) {
        PoaBaseObservation *obs = stList_get(node->observations, i);
        if (obs->readNo != readNo) {
            break;
        }
        i++;
    }
    return i;
}

int readSubstrings_cmpByQual(const void *a, const void *b) {
    /*
     * Compares read substrings by quality in descending order
     */
    BamChunkReadSubstring *rs1 = (BamChunkReadSubstring *) a;
    BamChunkReadSubstring *rs2 = (BamChunkReadSubstring *) b;

    return rs1->qualValue < rs2->qualValue ? 1 : (rs1->qualValue > rs2->qualValue ? -1 : 0);
}

stList *filterReadSubstrings(stList *readSubstrings, PolishParams *params) {
    // Sort the substrings by descending qvalue
    stList_sort(readSubstrings, readSubstrings_cmpByQual);

    while (stList_length(readSubstrings) > params->filterReadsWhileHaveAtLeastThisCoverage) {
        BamChunkReadSubstring *rs = stList_peek(readSubstrings);
        if (rs->qualValue >= params->minAvgBaseQuality ||
            rs->qualValue == -1) { //Filter by qvalue, but don't filter if some or all reads
            // don't have q-values
            break;
        }
        bamChunkReadSubstring_destruct(rs);
        stList_pop(readSubstrings);
    }

    return readSubstrings;
}

stList *getReadSubstrings2(stList *bamChunkReads, Poa *poa, int64_t from, int64_t to, PolishParams *params,
        bool shouldFilter) {
    /*
     * Get the substrings of reads aligned to the interval from (inclusive) to to
     * (exclusive) and their qual values. Adds them to readSubstrings and qualValues, respectively.
     */
    stList *readSubstrings = stList_construct3(0, (void (*)(void *)) bamChunkReadSubstring_destruct);

    // Deal with boundary cases
    if (from == 0) {
        if (to == stList_length(poa->nodes)) {
            // If from and to reference positions that bound the complete alignment just
            // copy the complete reads
            for (int64_t i = 0; i < stList_length(bamChunkReads); i++) {
                BamChunkRead *bamChunkRead = stList_get(bamChunkReads, i);
                stList_append(readSubstrings,
                              bamChunkRead_getSubstring(bamChunkRead, 0, bamChunkRead->rleRead->length, params));
            }
            return shouldFilter ? filterReadSubstrings(readSubstrings, params) : readSubstrings;
        }

        // Otherwise, include the read prefixes that end at to
        PoaNode *node = stList_get(poa->nodes, to);
        int64_t i = 0;
        while (i < stList_length(node->observations)) {
            PoaBaseObservation *obs = stList_get(node->observations, i);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, obs->readNo);
            // Trim the read substring, copy it and add to the substrings list
            stList_append(readSubstrings, bamChunkRead_getSubstring(bamChunkRead, 0, obs->offset, params));
            i = skipDupes(node, ++i, obs->readNo);
        }
        return shouldFilter ? filterReadSubstrings(readSubstrings, params) : readSubstrings;
    } else if (to == stList_length(poa->nodes)) {
        // Finally, include the read suffixs that start at from
        PoaNode *node = stList_get(poa->nodes, from);
        int64_t i = 0;
        while (i < stList_length(node->observations)) {
            PoaBaseObservation *obs = stList_get(node->observations, i);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, obs->readNo);
            // Trim the read substring, copy it and add to the substrings list
            stList_append(readSubstrings, bamChunkRead_getSubstring(bamChunkRead, obs->offset,
                                                                    bamChunkRead->rleRead->length - obs->offset,
                                                                    params));
            i = skipDupes(node, ++i, obs->readNo);
        }
        return shouldFilter ? filterReadSubstrings(readSubstrings, params) : readSubstrings;
    }

    PoaNode *fromNode = stList_get(poa->nodes, from);
    PoaNode *toNode = stList_get(poa->nodes, to);

    int64_t i = 0, j = 0;
    while (i < stList_length(fromNode->observations) && j < stList_length(toNode->observations)) {
        PoaBaseObservation *obsFrom = stList_get(fromNode->observations, i);
        PoaBaseObservation *obsTo = stList_get(toNode->observations, j);

        if (obsFrom->readNo == obsTo->readNo) {
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, obsFrom->readNo);
            if (obsTo->offset - obsFrom->offset > 0) { // If a non zero run of bases
                stList_append(readSubstrings,
                              bamChunkRead_getSubstring(bamChunkRead, obsFrom->offset, obsTo->offset - obsFrom->offset,
                                                        params));
            }
            i = skipDupes(fromNode, ++i, obsFrom->readNo);
            j = skipDupes(toNode, ++j, obsTo->readNo);
        } else if (obsFrom->readNo < obsTo->readNo) {
            i = skipDupes(fromNode, ++i, obsFrom->readNo);
        } else {
            assert(obsFrom->readNo > obsTo->readNo);
            j = skipDupes(toNode, ++j, obsTo->readNo);
        }
    }

    return shouldFilter ? filterReadSubstrings(readSubstrings, params) : readSubstrings;
}

stList *getReadSubstrings(stList *bamChunkReads, Poa *poa, int64_t from, int64_t to, PolishParams *params) {
    return getReadSubstrings2(bamChunkReads, poa, from, to, params, TRUE);
}

// Code to create anchors

double *getCandidateWeights(Poa *poa, PolishParams *params) {
    double *candidateWeights = st_calloc(stList_length(poa->nodes), sizeof(double));

    int64_t window = 100; // Size of window to average coverage over

    if (window >= stList_length(poa->nodes)) {
        double candidateWeight = getAvgCoverage(poa, 0, stList_length(poa->nodes)) * params->candidateVariantWeight;
        for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
            candidateWeights[i] = candidateWeight;
        }
        return candidateWeights;
    }

    double totalWeight = 0;
    for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
        totalWeight += getTotalWeight(poa, stList_get(poa->nodes, i));
        if (i >= window) {
            totalWeight -= getTotalWeight(poa, stList_get(poa->nodes, i - window));
            candidateWeights[i - window / 2] = totalWeight / window * params->candidateVariantWeight;
        }
    }

    // Fill in bounding bases
    for (int64_t i = 0; i < window / 2; i++) {
        candidateWeights[i] = candidateWeights[window / 2];
        candidateWeights[stList_length(poa->nodes) - 1 - i] = candidateWeights[stList_length(poa->nodes) - 1 -
                                                                               window / 2];
    }

    return candidateWeights;
}

bool *getCandidateVariantOverlapPositions(Poa *poa, double *candidateWeights) {
    /*
     * Return a boolean for each poaNode (as an array) indicating if the node is a candidate variant
     * site or is included in a candidate deletion.
     */

    bool *candidateVariantPositions = st_calloc(stList_length(poa->nodes), sizeof(bool));

    // Calculate positions that overlap candidate variants
    for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
        PoaNode *node = stList_get(poa->nodes, i);

        // Mark as variant if has a candidate substitution, repeat count change or an insert starts at this position
        if (hasCandidateSubstitution(poa, node, candidateWeights[i])
            || hasCandidateRepeatCountChange(poa, node, candidateWeights[i])
            || hasCandidateInsert(node, candidateWeights[i])) {
            candidateVariantPositions[i] = 1;
        }

        int64_t j = maxCandidateDeleteLength(node, candidateWeights[i]);
        if (j > 0) { // Mark as variant if precedes the start of a deletion
            candidateVariantPositions[i] = 1;
        }
        // Mark as variant if is included in candidate deletion
        while (j > 0) {
            assert(i + j < stList_length(poa->nodes));
            candidateVariantPositions[i + (j--)] = 1;
        }
    }

    return candidateVariantPositions;
}

bool *expand(bool *b, int64_t length, int64_t expansion) {
    /*
     * Returns a bool array in which a position is non-zero if a position
     * in b +/-i expansion is non-zero.
     */
    bool *b2 = st_calloc(length, sizeof(bool));
    for (int64_t i = 0; i < length; i++) {
        if (b[i]) {
            for (int64_t j = i - expansion; j < i + expansion; j++) {
                if (j >= 0 && j < length) {
                    b2[j] = 1;
                }
            }
        }
    }

    return b2;
}

void updateCandidateVariantPositionsByVcfEntries(bool *candidateVariantPositions, int64_t cvpLength, stList *vcfEntries) {
    stListIterator *itor = stList_getIterator(vcfEntries);
    VcfEntry *vcfEntry = stList_getNext(itor);

    int64_t tn = 0;
    int64_t fn = 0;
    int64_t fp = 0;
    int64_t tp = 0;
    for (int64_t i = 0; i < cvpLength; i++) {
        bool isCVP = candidateVariantPositions[i];
        bool isVCF = vcfEntry != NULL && vcfEntry->refPos == i;

        if (isCVP) {
            if (isVCF)  tp++;
            else        fp++;
        } else {
            if (isVCF)  fn++;
            else        tn++;
        }

        candidateVariantPositions[i] = isVCF;
        if (vcfEntry != NULL && vcfEntry->refPos <= i) {
            vcfEntry = stList_getNext(itor);
        }
    }
    char *logIdentifier = getLogIdentifier();
    st_logInfo(" %s Of %"PRId64" positions, got %"PRId64" TP, %"PRId64" FP, %"PRId64" FN, %"PRId64" TN, "
               "equating to a variation rate of %.5f, precision of %.5f, and recall of %.5f.\n",
            logIdentifier, cvpLength, tp, fp, fn, tn, 1.0*(tp+fn)/cvpLength, 1.0*tp/(tp+fp), 1.0*tp/(tp+fn));
    free(logIdentifier);

    stList_destructIterator(itor);
}


bool *getFilteredAnchorPositions(Poa *poa, double *candidateWeights, stList *vcfEntries, PolishParams *params) {
    /*
     * Create set of anchor positions, using positions not close to candidate variants
     */
    // Identity sites that overlap candidate variants, and expand to surrounding positions
    bool *candidateVariantPositions = getCandidateVariantOverlapPositions(poa, candidateWeights);

    if (vcfEntries != NULL) {
        updateCandidateVariantPositionsByVcfEntries(candidateVariantPositions, stList_length(poa->nodes), vcfEntries);
    }

    bool *expandedCandidateVariantPositions = expand(candidateVariantPositions, stList_length(poa->nodes),
                                                     params->columnAnchorTrim);

    // Anchors are those that are not close to expanded candidate variant positions
    bool *anchors = st_calloc(stList_length(poa->nodes), sizeof(bool));
    for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
        anchors[i] = !expandedCandidateVariantPositions[i];
    }

    // Cleanup
    free(candidateVariantPositions);
    free(expandedCandidateVariantPositions);

    // Log some stuff about the anchors
    if (st_getLogLevel() >= debug) {
        int64_t totalAnchorNo = 0;
        for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
            totalAnchorNo += anchors[i] ? 1 : 0;
        }
        st_logDebug(
                "Creating filtered anchor positions got: %" PRIi64 " anchors for ref seq of length: %" PRIi64 ", that's one every: %f bases\n",
                totalAnchorNo, stList_length(poa->nodes), stList_length(poa->nodes) / (double) totalAnchorNo);
    }

    return anchors;
}

uint64_t rleString_stringKey(const void *k) {
    return stHash_stringKey(((RleString *) k)->rleString);
}

int rleString_stringEqualKey(const void *key1, const void *key2) {
    return stString_eq(((RleString *) key1)->rleString, ((RleString *) key2)->rleString);
}

int rleString_expandedStringEqualKey(const void *key1, const void *key2) {
    if (!rleString_stringEqualKey(key1, key2)) {
        return 0;
    }
    RleString *r1 = (RleString *) key1;
    RleString *r2 = (RleString *) key2;
    if (r1->length != r2->length) {
        return 0;
    }
    for (int64_t i = 0; i < r1->length; i++) {
        if (r1->repeatCounts[i] != r2->repeatCounts[i]) {
            return 0;
        }
    }
    return 1;
}

stHash *groupRleStrings(stList *rleStrings) {
    /*
     * Input is a list of RleStrings. Returns a map whose keys are the compressed RLE strings and whose values are lists of the RleStrings with the given
     * compressed RLE string.
     */

    stHash *h = stHash_construct3(rleString_stringKey,
                                  rleString_expandedStringEqualKey, //rleString_stringEqualKey, //rleString_expandedStringEqualKey, //rleString_stringEqualKey, //
                                  NULL, (void (*)(void *)) stList_destruct);

    for (uint64_t i = 0; i < stList_length(rleStrings); i++) {
        RleString *rleString = stList_get(rleStrings, i);
        stList *l = stHash_search(h, rleString);
        if (l == NULL) {
            l = stList_construct();
            stHash_insert(h, rleString, l);
        }
        stList_append(l, rleString);
    }

    return h;
}

RleString *getConsensusRleString(stList *rleStrings) {
    /*
     * For a list of RleStrings all with the same RLE string return a
     * consensus RleString with consensus repeat counts.
     */
    assert(stList_length(rleStrings) > 0);
    RleString *r = stList_peek(rleStrings);
    uint8_t repeatCounts[r->length];

    for (int64_t j = 0; j < r->length; j++) {
        uint64_t k = 0;
        for (int64_t i = 0; i < stList_length(rleStrings); i++) {
            RleString *s = stList_get(rleStrings, i);
            assert(s->length == r->length);
            k += s->repeatCounts[j];
        }
        k = roundf(((float) k) / stList_length(rleStrings));
        repeatCounts[j] = k == 0 ? 1 : (k > 255 ? 255 : k);
    }

    return rleString_constructPreComputed(r->rleString, repeatCounts);
}

stList *getCandidateAllelesFromReadSubstrings(stList *readSubstrings, PolishParams *p) {
    // Get the rle strings for the bamChunkReadSubstrings
    stList *rleStrings = stList_construct3(0, (void (*)(void *)) rleString_destruct);
    for (int64_t i = 0; i < stList_length(readSubstrings); i++) {
        BamChunkReadSubstring *r = stList_get(readSubstrings, i);
        stList_append(rleStrings, bamChunkReadSubstring_getRleString(r));
    }

    // Group together the RleString by RLE string
    stHash *h = groupRleStrings(rleStrings);

    // For each RLE string get the consensus, expanded allele string
    stHashIterator *it = stHash_getIterator(h);
    RleString *rleString;
    stList *alleles = stList_construct3(0, (void (*)(void *)) free);
    while ((rleString = stHash_getNext(it)) != NULL) {
        stList *l = stHash_search(h, rleString);
        assert(l != NULL);
        //if(stList_length(l) > 1) {
        RleString *r = getConsensusRleString(l);
        stList_append(alleles, rleString_expand(r));
        rleString_destruct(r);
        //}
    }

    // Cleanup
    stHash_destructIterator(it);
    stHash_destruct(h);
    stList_destruct(rleStrings);

    return alleles;
}


void bubble_destruct(Bubble b) {
    // Cleanup the reads
    for (int64_t j = 0; j < b.readNo; j++) {
        free(b.reads[j]);
    }
    free(b.reads);
    // Cleanup the alleles
    for (int64_t j = 0; j < b.alleleNo; j++) {
        rleString_destruct(b.alleles[j]);
    }
    free(b.alleles);
    // Cleanup the allele supports
    free(b.alleleReadSupports);
    // Cleanup the reference allele
    rleString_destruct(b.refAllele);
}

void bubbleGraph_destruct(BubbleGraph *bg) {
    // Clean up the memory for each bubble
    for (int64_t i = 0; i < bg->bubbleNo; i++) {
        bubble_destruct(bg->bubbles[i]);
    }
    free(bg->bubbles);
    free(bg);
}


BubbleGraph *bubbleGraph_constructFromPoa(Poa *poa, stList *bamChunkReads, PolishParams *params) {
    return bubbleGraph_constructFromPoa2(poa, bamChunkReads, params, FALSE);
}

BubbleGraph *bubbleGraph_constructFromPoa2(Poa *poa, stList *bamChunkReads, PolishParams *params, bool phasing) {
    return bubbleGraph_constructFromPoaAndVCF(poa, bamChunkReads, NULL, params, phasing);
}

BubbleGraph *bubbleGraph_constructFromPoaAndVCF(Poa *poa, stList *bamChunkReads, stList *vcfEntries,
                                                PolishParams *params, bool phasing) {
    // Setup
    double *candidateWeights = getCandidateWeights(poa, params);

    // Log info about the alignment
    if (st_getLogLevel() >= info) {
        double avgCoverage = getAvgCoverage(poa, 0, stList_length(poa->nodes));
        double totalCandidateWeight = 0.0;
        for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
            totalCandidateWeight += candidateWeights[i];
        }
        st_logDebug("Got avg. coverage: %f for region of length: %" PRIi64 " and avg. candidate weight of: %f\n",
                    avgCoverage / PAIR_ALIGNMENT_PROB_1, stList_length(poa->nodes),
                    totalCandidateWeight / (PAIR_ALIGNMENT_PROB_1 * stList_length(poa->nodes)));
    }

    // Sort the base observations to make the getReadSubstrings function work
    sortBaseObservations(poa);

    // Identify anchor points, represented as a binary array, one bit for each POA node
    bool *anchors = getFilteredAnchorPositions(poa, candidateWeights, vcfEntries, params);

    // Make a list of bubbles
    stList *bubbles = stList_construct3(0, free);
    int64_t pAnchor = 0; // Previous anchor, starting from first position of POA, which is the prefix "N"
    for (int64_t i = 1; i < stList_length(poa->nodes); i++) {
        if (anchors[i]) { // If position i is an anchor
            assert(i > pAnchor);
            if (i - pAnchor != 1) { // In case anchors are not trivially adjacent there exists a potential bubble
                // with start coordinate on the reference sequence of pAnchor and length pAnchor-i

                // Get read substrings
                stList *readSubstrings = getReadSubstrings(bamChunkReads, poa, pAnchor+1, i, params);

                if(stList_length(readSubstrings) > 0) {
                    stList *alleles = NULL;
                    bool useReadAlleles = phasing ? params->useReadAllelesInPhasing : params->useReadAlleles;
                    if(useReadAlleles) {
                        alleles = getCandidateAllelesFromReadSubstrings(readSubstrings, params);
                    } else {
                        // Calculate the list of alleles
                        double weightAdjustment = 1.0;
                        do {
                            alleles = getCandidateConsensusSubstrings(poa, pAnchor+1, i,
                                                                      candidateWeights, weightAdjustment, params->maxConsensusStrings);
                            weightAdjustment *= 1.5; // Increase the candidate weight by 50%
                        } while(alleles == NULL);
                    }

                    // Get existing reference string
                    assert(i - 1 - pAnchor > 0);
                    RleString *existingRefSubstring = rleString_copySubstring(poa->refString, pAnchor, i - 1 - pAnchor);
                    assert(existingRefSubstring->length == i - pAnchor - 1);
                    char *expandedExistingRefSubstring = rleString_expand(existingRefSubstring);

                    // Check if the reference allele is in the set of alleles and add it if not
                    bool seenRefAllele = 0;
                    for (int64_t j = 0; j < stList_length(alleles); j++) {
                        if (stString_eq(expandedExistingRefSubstring, stList_get(alleles, j))) {
                            seenRefAllele = 1;
                            break;
                        }
                    }
                    if (!seenRefAllele) {
                        stList_append(alleles, stString_copy(expandedExistingRefSubstring));
                    }

                    // If it is not trivial because it contains more than one allele
                    if (stList_length(alleles) > 1) {

                        Bubble *b = st_malloc(sizeof(Bubble)); // Make a bubble and add to list of bubbles
                        stList_append(bubbles, b);

                        // Set the coordinates
                        b->refStart = pAnchor;
                        b->bubbleLength = i - 1 - pAnchor;

                        // The reference allele
                        b->refAllele = existingRefSubstring;

                        // Add read substrings
                        b->readNo = stList_length(readSubstrings);
                        b->reads = st_malloc(sizeof(BamChunkReadSubstring *) * b->readNo);
                        for (int64_t j = 0; j < b->readNo; j++) {
                            b->reads[j] = stList_pop(readSubstrings);
                        }

                        // Now copy the alleles list to the bubble's array of alleles
                        b->alleleNo = stList_length(alleles);
                        b->alleles = st_malloc(sizeof(RleString *) * b->alleleNo);
                        for (int64_t j = 0; j < b->alleleNo; j++) {
                            b->alleles[j] = params->useRunLengthEncoding ? rleString_construct(stList_get(alleles, j))
                                                                         : rleString_construct_no_rle(
                                            stList_get(alleles, j));
                        }

                        // Get allele supports
                        b->alleleReadSupports = st_calloc(b->readNo * b->alleleNo, sizeof(float));

                        stList *anchorPairs = stList_construct(); // Currently empty

                        SymbolString alleleSymbolStrings[b->alleleNo];
                        for (int64_t j = 0; j < b->alleleNo; j++) {
                            alleleSymbolStrings[j] = rleString_constructSymbolString(b->alleles[j], 0,
                                                                                     b->alleles[j]->length,
                                                                                     params->alphabet,
                                                                                     params->useRepeatCountsInAlignment,
                                                                                     poa->maxRepeatCount);
                        }

                        stHash *cachedScores = stHash_construct3(rleString_stringKey, rleString_expandedStringEqualKey,
                                                                 (void (*)(void *)) rleString_destruct, free);

                        for (int64_t k = 0; k < b->readNo; k++) {
                            RleString *readSubstring = bamChunkReadSubstring_getRleString(b->reads[k]);
                            SymbolString rS = rleString_constructSymbolString(readSubstring, 0, readSubstring->length,
                                                                              params->alphabet,
                                                                              params->useRepeatCountsInAlignment,
                                                                              poa->maxRepeatCount);
                            StateMachine *sM = b->reads[k]->read->forwardStrand
                                               ? params->stateMachineForForwardStrandRead
                                               : params->stateMachineForReverseStrandRead;

                            uint64_t *index = stHash_search(cachedScores, readSubstring);
                            if (index != NULL) {
                                for (int64_t j = 0; j < b->alleleNo; j++) {
                                    b->alleleReadSupports[j * b->readNo + k] = b->alleleReadSupports[j * b->readNo +
                                                                                                     *index];
                                }
                                rleString_destruct(readSubstring);
                            } else {
                                index = st_malloc(sizeof(uint64_t));
                                *index = k;
                                stHash_insert(cachedScores, readSubstring, index);
                                for (int64_t j = 0; j < b->alleleNo; j++) {
                                    b->alleleReadSupports[j * b->readNo + k] = computeForwardProbability(
                                            alleleSymbolStrings[j], rS, anchorPairs, params->p, sM, 0, 0);
                                }
                            }

                            symbolString_destruct(rS);
                        }

                        stHash_destruct(cachedScores);
                        for (int64_t j = 0; j < b->alleleNo; j++) {
                            symbolString_destruct(alleleSymbolStrings[j]);
                        }
                        stList_destruct(anchorPairs);
                    }
                        // Cleanup
                    else {
                        rleString_destruct(existingRefSubstring);
                    }

                    free(expandedExistingRefSubstring);
                    stList_destruct(alleles);
                    stList_destruct(readSubstrings);
                }
            }
            // Update previous anchor
            pAnchor = i;
        }
    }

    // Build the the graph

    BubbleGraph *bg = st_malloc(sizeof(BubbleGraph));
    bg->refString = poa->refString;

    // Copy the bubbles
    bg->bubbleNo = stList_length(bubbles);
    bg->bubbles = st_calloc(bg->bubbleNo, sizeof(Bubble)); // allocate bubbles
    for (int64_t i = 0; i < bg->bubbleNo; i++) {
        bg->bubbles[i] = *(Bubble *) stList_get(bubbles, i);
    }

    // Fill in the bubble allele offsets
    int64_t alleleOffset = 0;
    for (int64_t i = 0; i < bg->bubbleNo; i++) {
        bg->bubbles[i].alleleOffset = alleleOffset;
        alleleOffset += bg->bubbles[i].alleleNo;
    }
    bg->totalAlleles = alleleOffset;

    // Cleanup
    free(anchors);
    free(candidateWeights);
    stList_destruct(bubbles);

    return bg;
}


BubbleGraph *bubbleGraph_partitionFilteredReads(Poa *poa, stList *bamChunkReads, stGenomeFragment *gF,
                                                BubbleGraph *bg, BamChunk *bamChunk, uint64_t *reference_rleToNonRleCoordMap,
                                                stSet *hap1Reads, stSet *hap2Reads, PolishParams *params, FILE *out,
                                                char *logIdentifier) {
    // our eventual scores
    stHash *totalReadScore_hap1 = stHash_construct2(NULL, free);
    stHash *totalReadScore_hap2 = stHash_construct2(NULL, free);
    for (int64_t i = 0; i < stList_length(bamChunkReads); i++) {
        stHash_insert(totalReadScore_hap1, stList_get(bamChunkReads, i), st_calloc(1, sizeof(double)));
        stHash_insert(totalReadScore_hap2, stList_get(bamChunkReads, i), st_calloc(1, sizeof(double)));
    }

    // save to output
    if (out != NULL) fprintf(out, ",\n \"filtered\": [");
    bool firstBubble = TRUE;

    // loop over all primary bubbles
    for (uint64_t primaryBubbleIdx = 0; primaryBubbleIdx < gF->length; primaryBubbleIdx++) {

        // bubble and hap info
        Bubble *primaryBubble = &bg->bubbles[gF->refStart + primaryBubbleIdx];
        int64_t hap1AlleleNo = gF->haplotypeString1[primaryBubbleIdx];
        int64_t hap2AlleleNo = gF->haplotypeString2[primaryBubbleIdx];

        RleString *hap1 = primaryBubble->alleles[hap1AlleleNo];
        RleString *hap2 = primaryBubble->alleles[hap2AlleleNo];

        // we only care about hets
        if (hap1 == hap2) continue;

        // anchor positions
        uint64_t refStart = primaryBubble->refStart;

        // make allele list from primary haplotype alleles
        stList *alleles = stList_construct3(0, free);
        stList_append(alleles, rleString_expand(hap1));
        stList_append(alleles, rleString_expand(hap2));

        // get read substrings
        stList *readSubstrings = getReadSubstrings2(bamChunkReads, poa, refStart+1, refStart+primaryBubble->bubbleLength + 1,
                params, FALSE);

        // Get existing reference string
        //todo remove this?
        RleString *existingRefSubstring = rleString_copySubstring(poa->refString, refStart, primaryBubble->bubbleLength);
        assert(existingRefSubstring->length == primaryBubble->bubbleLength);
        char *expandedExistingRefSubstring = rleString_expand(existingRefSubstring);

        // Check if the reference allele is in the set of alleles and add it if not
        bool seenRefAllele = 0;
        for (int64_t j = 0; j < stList_length(alleles); j++) {
            if (stString_eq(expandedExistingRefSubstring, stList_get(alleles, j))) {
                seenRefAllele = 1;
                break;
            }
        }
        if (!seenRefAllele) {
            stList_append(alleles, stString_copy(expandedExistingRefSubstring));
        }

        // ensure we do have the alleles we expect
        assert(stList_length(alleles) == 2 || stList_length(alleles) == 3);


        Bubble *b = st_malloc(sizeof(Bubble)); // Make a bubble and add to list of bubbles

        // Set the coordinates
        b->refStart = (uint64_t) refStart;

        // The reference allele
        b->refAllele = existingRefSubstring;

        // Add read substrings
        b->readNo = (uint64_t) stList_length(readSubstrings);
        b->reads = st_malloc(sizeof(BamChunkReadSubstring *) * b->readNo);
        for (int64_t j = 0; j < b->readNo; j++) {
            b->reads[j] = stList_pop(readSubstrings);
        }

        // Now copy the alleles list to the bubble's array of alleles
        b->alleleNo = (uint64_t) stList_length(alleles);
        b->alleles = st_malloc(sizeof(RleString *) * b->alleleNo);
        for (int64_t j = 0; j < b->alleleNo; j++) {
            b->alleles[j] = params->useRunLengthEncoding ? rleString_construct(stList_get(alleles, j))
                                                         : rleString_construct_no_rle(stList_get(alleles, j));
        }

        // Get allele supports
        b->alleleReadSupports = st_calloc(b->readNo * b->alleleNo, sizeof(float));

        stList *anchorPairs = stList_construct(); // Currently empty

        SymbolString alleleSymbolStrings[b->alleleNo];
        for (int64_t j = 0; j < b->alleleNo; j++) {
            alleleSymbolStrings[j] = rleString_constructSymbolString(b->alleles[j], 0,
                                                                     b->alleles[j]->length,
                                                                     params->alphabet,
                                                                     params->useRepeatCountsInAlignment,
                                                                     poa->maxRepeatCount);
        }

        // get alignment likelihoods
        stHash *cachedScores = stHash_construct3(rleString_stringKey, rleString_expandedStringEqualKey,
                                                 (void (*)(void *)) rleString_destruct, free);
        for (int64_t k = 0; k < b->readNo; k++) {
            RleString *readSubstring = bamChunkReadSubstring_getRleString(b->reads[k]);
            SymbolString rS = rleString_constructSymbolString(readSubstring, 0, readSubstring->length,
                                                              params->alphabet,
                                                              params->useRepeatCountsInAlignment,
                                                              poa->maxRepeatCount);
            StateMachine *sM = b->reads[k]->read->forwardStrand
                               ? params->stateMachineForForwardStrandRead
                               : params->stateMachineForReverseStrandRead;

            uint64_t *index = stHash_search(cachedScores, readSubstring);
            if (index != NULL) {
                for (int64_t j = 0; j < b->alleleNo; j++) {
                    b->alleleReadSupports[j * b->readNo + k] = b->alleleReadSupports[j * b->readNo +
                                                                                     *index];
                }
                rleString_destruct(readSubstring);
            } else {
                index = st_malloc(sizeof(uint64_t));
                *index = (uint64_t) k;
                stHash_insert(cachedScores, readSubstring, index);
                for (int64_t j = 0; j < b->alleleNo; j++) {
                    b->alleleReadSupports[j * b->readNo + k] =
                            (float) computeForwardProbability(alleleSymbolStrings[j], rS, anchorPairs, params->p, sM, 0, 0);
                }
            }

            symbolString_destruct(rS);
        }



        // write to output
        if (out != NULL) {
            if (firstBubble) {
                fprintf(out, "\n  {\n");
                firstBubble = FALSE;
            } else {
                fprintf(out, ",\n  {\n");
            }
            int64_t trueRefStartPos = bamChunk->chunkBoundaryStart + reference_rleToNonRleCoordMap[b->refStart];
            fprintf(out, "   \"refPos\": %"PRId64",\n", trueRefStartPos);
            fprintf(out, "   \"reads\": [");
        }

        // rank reads for each bubble
        for (int64_t k = 0; k < b->readNo; k++) {
            BamChunkReadSubstring *bcrss = b->reads[k];
            BamChunkRead *bcr = bcrss->read;
            float supportHap1 = b->alleleReadSupports[0 * b->readNo + k];
            float supportHap2 = b->alleleReadSupports[1 * b->readNo + k];

            double *currRS = stHash_search(totalReadScore_hap1, bcr);
            *currRS += supportHap1 - stMath_logAddExact(supportHap1, supportHap2);
            currRS = stHash_search(totalReadScore_hap2, bcr);
            *currRS += supportHap2 - stMath_logAddExact(supportHap2, supportHap1);

            // write to output
            if (out != NULL) {
                if (k != 0) fprintf(out, ",");
                fprintf(out, "\n    {\n");
                fprintf(out, "     \"name\": \"%s\",\n", bcrss->read->readName);
                fprintf(out, "     \"qual\": %f,\n", bcrss->qualValue);
                fprintf(out, "     \"hapSupportH1\": %f,\n", supportHap1);
                fprintf(out, "     \"hapSupportH2\": %f\n", supportHap2);
                fprintf(out, "    }");
            }
        }

        // write to output
        if (out != NULL) {
            fprintf(out, "\n   ]");
            fprintf(out, "\n  }");
        }

        // cleanup
        stHash_destruct(cachedScores);
        for (int64_t j = 0; j < b->alleleNo; j++) {
            symbolString_destruct(alleleSymbolStrings[j]);
        }
        stList_destruct(anchorPairs);
        free(expandedExistingRefSubstring);
        stList_destruct(alleles);
        stList_destruct(readSubstrings);
        bubble_destruct(*b);
        free(b);
    }

    // write to output
    fprintf(out, "\n ]");

    // get scores and save to appropriate sets
    int64_t totalNoScoreLength = 0;
    int64_t noScoreCount = 0;
    int64_t unclassifiedCount = 0;
    int64_t hap1Count = 0;
    int64_t hap2Count = 0;
    for (int i = 0; i < stList_length(bamChunkReads); i++) {
        BamChunkRead *bcr = stList_get(bamChunkReads, i);
        double *totalSupportH1 = stHash_search(totalReadScore_hap1, bcr);
        double *totalSupportH2 = stHash_search(totalReadScore_hap2, bcr);

        if (*totalSupportH1 > *totalSupportH2) {
            stSet_insert(hap1Reads, bcr);
            hap1Count++;
        } else if (*totalSupportH2 > *totalSupportH1)  {
            stSet_insert(hap2Reads, bcr);
            hap2Count++;
        } else {
            if (*totalSupportH1 == 0) {
                totalNoScoreLength += bcr->rleRead->nonRleLength;
                noScoreCount++;
            }
            unclassifiedCount++;
        }
    }

    // loggit
    int64_t length = stList_length(bamChunkReads);
    st_logInfo(" %s Of %"PRId64" filtered reads: %"PRId64" (%.2f) were hap1, %"PRId64" (%.2f) were hap2, %"PRId64" (%.2f) were unclassified with %"PRId64" (%.2f) having no score (avg len %"PRId64").\n",
               logIdentifier, length, hap1Count, 1.0*hap1Count/length, hap2Count, 1.0*hap2Count/length,
               unclassifiedCount, 1.0*unclassifiedCount/length, noScoreCount,
               1.0*noScoreCount/(unclassifiedCount == 0 ? 1 : unclassifiedCount),
               totalNoScoreLength / (noScoreCount == 0 ? 1 : noScoreCount));


    // other cleanup
    stHash_destruct(totalReadScore_hap1);
    stHash_destruct(totalReadScore_hap2);
}

stHash *bubbleGraph_getProfileSeqs(BubbleGraph *bg, stReference *ref) {
    // First calculate the length of all the profile sequences

    stHash *readEnds = stHash_construct2(NULL, free); // The last bubble the read is observed to be part of

    // For each bubble in the bubble graph
    for (uint64_t i = 0; i < bg->bubbleNo; i++) {
        Bubble *b = &(bg->bubbles[i]);

        // For each read aligned to this bubble
        for (uint64_t j = 0; j < b->readNo; j++) {
            BamChunkReadSubstring *s = b->reads[j];
            assert(s->read != NULL);

            // Look up the corresponding profile sequence
            uint64_t *k = stHash_search(readEnds, s->read);

            if (k == NULL) { // We are starting a new read, so make a new
                // length entry
                k = st_calloc(1, sizeof(uint64_t));
                stHash_insert(readEnds, s->read, k);
            }

            k[0] = i; // Update last time read observed aligned to a bubble
        }
    }

    stHash *readsToPSeqs = stHash_construct();

    // Now build the profile sequences

    // For each bubble in the bubble graph
    for (uint64_t i = 0; i < bg->bubbleNo; i++) {
        Bubble *b = &(bg->bubbles[i]);

        // For each read aligned to this bubble
        for (uint64_t j = 0; j < b->readNo; j++) {
            BamChunkReadSubstring *s = b->reads[j];
            assert(stHash_search(readEnds, s->read) != NULL);

            // Look up the corresponding profile sequence
            stProfileSeq *pSeq = stHash_search(readsToPSeqs, s->read);

            if (pSeq == NULL) { // We are starting a new read, so make a new
                // profile sequence

                // Calculate the length in bubbles of the profile sequence
                uint64_t *k = stHash_search(readEnds, s->read);
                assert(k != NULL);
                assert(i <= k[0]); // The first bubble the read is part of must precede or be equal to the last
                uint64_t pSeqLength = k[0] - i + 1;
                assert(i + pSeqLength <= bg->bubbleNo);

                pSeq = stProfileSeq_constructEmptyProfile(ref, s->read->readName, i, pSeqLength);
                stHash_insert(readsToPSeqs, s->read, pSeq);
            }

            // Sanity check the pSeq
            assert(b->alleleOffset >= pSeq->alleleOffset);
            assert(i < pSeq->refStart + pSeq->length);

            // For each allele in bubble add the prob that the read was generated by
            // the read

            // First calculate the total log probability of the read given the alleles, to normalize
            // the log probabilities
            // This acts as a normalizing constant
            double totalLogProb = LOG_ZERO;
            for (uint64_t k = 0; k < b->alleleNo; k++) {
                totalLogProb = stMath_logAddExact(totalLogProb, b->alleleReadSupports[b->readNo * k + j]);
            }

            // Normalize prob by totalLogProb
            uint64_t alleleOffset = b->alleleOffset - pSeq->alleleOffset;
            for (uint64_t k = 0; k < b->alleleNo; k++) {
                float logProb = b->alleleReadSupports[b->readNo * k + j];
                assert(logProb <= totalLogProb);
                int64_t l = roundf(PROFILE_PROB_SCALAR * (totalLogProb - logProb));
                assert(l >= 0);
                pSeq->profileProbs[alleleOffset + k] = l > 255 ? 255 : l;
            }
        }
    }

    // Cleanup
    stHash_destruct(readEnds);

    return readsToPSeqs;
}

stReference *bubbleGraph_getReference(BubbleGraph *bg, char *refName, Params *params) {
    stReference *ref = st_calloc(1, sizeof(stReference));

    ref->referenceName = stString_copy(refName);
    ref->length = bg->bubbleNo;
    ref->sites = st_calloc(bg->bubbleNo, sizeof(stSite));
    ref->totalAlleles = 0;

    //stList *anchorPairs = stList_construct(); // Currently empty, and no anchor pairs will be created
    for (uint64_t i = 0; i < bg->bubbleNo; i++) {
        Bubble *b = &bg->bubbles[i];
        stSite *s = &ref->sites[i];
        s->alleleNumber = b->alleleNo;
        s->alleleOffset = b->alleleOffset;
        ref->totalAlleles += b->alleleNo;
        s->allelePriorLogProbs = st_calloc(b->alleleNo, sizeof(uint16_t)); // These are all set equal
        s->substitutionLogProbs = st_calloc(b->alleleNo * b->alleleNo, sizeof(uint16_t));

        for (uint64_t j = 0; j < b->alleleNo; j++) {
            for (uint64_t k = 0; k < b->alleleNo; k++) {
                s->substitutionLogProbs[j * b->alleleNo + k] =
                        j == k ? 0 : roundf(
                                -log(params->polishParams->hetSubstitutionProbability) * PROFILE_PROB_SCALAR); //l;
            }
        }
    }

    return ref;
}

/*
 * Phasing of bubble graphs
 */

void bubbleGraph_logPhasedBubbleGraph(BubbleGraph *bg, stRPHmm *hmm, stList *path,
                                      stHash *readsToPSeqs, stList *profileSeqs, stGenomeFragment *gF) {
    /*
     * Sanity checks / logging for phased bubble graph
     */

    if (st_getLogLevel() == debug) {
        // Check read partition is complete
        assert(stSet_size(gF->reads1) + stSet_size(gF->reads2) == stList_length(profileSeqs));
        stSet *intersection = stSet_getIntersection(gF->reads1, gF->reads2);
        assert(stSet_size(intersection) == 0);
        stSet_destruct(intersection);

        stRPColumn *column = hmm->firstColumn;
        assert(column->length > 0);
        uint64_t colIndex = 0, colCo = 0;

        for (uint64_t i = 0; i < gF->length; i++) {
            assert(column != NULL);
            Bubble *b = &bg->bubbles[gF->refStart + i];

            stSite *s = &(hmm->ref->sites[gF->refStart + i]);
            assert(s->alleleNumber == b->alleleNo);

            assert(gF->haplotypeString1[i] < b->alleleNo);
            assert(gF->haplotypeString2[i] < b->alleleNo);
            //assert(column->depth >= b->readNo);

            RleString *hap1 = b->alleles[gF->haplotypeString1[i]];
            RleString *hap2 = b->alleles[gF->haplotypeString2[i]];

            if (gF->haplotypeString1[i] != gF->haplotypeString2[i] || !rleString_eq(b->refAllele, hap1)) {
                stRPCell *cell = stList_get(path, colIndex);

                double strandSkew = bubble_phasedStrandSkew(b, readsToPSeqs, gF);

                st_logDebug(
                        ">>Phasing Bubble Graph: (Het: %s) At site: %i (of %i) with %i potential alleles got %s (%i) (log-prob: %f) for hap1 with %i reads and %s (%i) (log-prob: %f) for hap2 with %i reads (total depth %i), and ancestral allele %s (%i), genotype prob: %f, strand-skew p-value: %f\n",
                        gF->haplotypeString1[i] != gF->haplotypeString2[i] ? "True" : "False", (int) i,
                        (int) gF->length, (int) b->alleleNo,
                        b->alleles[gF->haplotypeString1[i]]->rleString, (int) gF->haplotypeString1[i],
                        gF->haplotypeProbs1[i], popcount64(cell->partition),
                        b->alleles[gF->haplotypeString2[i]]->rleString, (int) gF->haplotypeString2[i],
                        gF->haplotypeProbs2[i], (int) (column->depth - popcount64(cell->partition)),
                        (int) column->depth,
                        b->alleles[gF->ancestorString[i]]->rleString, (int) gF->ancestorString[i], gF->genotypeProbs[i],
                        (float) strandSkew);

                double strandSkews[b->alleleNo];
                bubble_calculateStrandSkews(b, strandSkews);

                for (uint64_t j = 0; j < b->alleleNo; j++) {
                    st_logDebug("\t>>Allele %i (ref allele: %s)\t strand-skew: %+.5f \t", (int) j,
                                rleString_eq(b->refAllele, b->alleles[j]) ? "True " : "False",
                                (float) strandSkews[j]);
                    rleString_print(b->alleles[j], stderr);
                    char *expandedAllele = rleString_expand(b->alleles[j]);
                    st_logDebug("\t%s\t", expandedAllele);
                    free(expandedAllele);
                    for (uint64_t k = 0; k < b->alleleNo; k++) {
                        st_logDebug("%i \t", (int) s->substitutionLogProbs[j * b->alleleNo + k]);
                    }
                    st_logDebug("\n");
                }


                for (uint64_t k = 0; k < 2; k++) {
                    uint64_t l = 0;
                    float supports[b->alleleNo];
                    for (uint64_t j = 0; j < b->alleleNo; j++) {
                        supports[j] = 0.0;
                    }

                    for (uint64_t j = 0; j < b->readNo; j++) {
                        BamChunkReadSubstring *s = b->reads[j];
                        stProfileSeq *pSeq = stHash_search(readsToPSeqs, s->read);
                        assert(pSeq != NULL);
                        if (stSet_search(k == 0 ? gF->reads1 : gF->reads2, pSeq) != NULL) {
                            st_logDebug("\t\t>>Partition %i, read %3i:\t strand %i\t ", (int) k + 1, (int) l++,
                                        (int) s->read->forwardStrand);

                            for (uint64_t m = 0; m < b->alleleNo; m++) {
                                st_logDebug("%+8.5f\t", b->alleleReadSupports[m * b->readNo + j]);
                                supports[m] += b->alleleReadSupports[m * b->readNo + j];
                            }

                            RleString *readSubstring = bamChunkReadSubstring_getRleString(s);
                            //st_logDebug("%s\n", readSubstring->rleString);
                            rleString_print(readSubstring, stderr);
                            st_logDebug(" qv: %8.5f\n", (float) s->qualValue);
                            rleString_destruct(readSubstring);
                        }
                    }

                    st_logDebug("\t\tCombined allele partition supports:\n");
                    st_logDebug("\t\t\t");
                    for (uint64_t j = 0; j < b->alleleNo; j++) {
                        st_logDebug("%8.5f\t", supports[j]);
                    }
                    st_logDebug("\n");

                }
            }

            if (++colCo >= column->length) {
                colCo = 0;
                colIndex++;
                column = colIndex < stList_length(path) ? column->nColumn->nColumn : NULL;
                assert(column == NULL || column->length > 0);
            }
        }
        assert(colIndex == stList_length(path));

        st_logDebug(">>Fraction of bubbles skewed %f (of %i total)\n",
                    (float) bubbleGraph_skewedBubbles(bg, readsToPSeqs, gF), (int) bg->bubbleNo);
    }
}


/*
 * Phasing of bubble graphs
 */

void bubbleGraph_saveBubblePhasingInfo(BamChunk *bamChunk, BubbleGraph *bg, stHash *readsToPSeqs, stGenomeFragment *gF,
                                       uint64_t *reference_rleToNonRleCoordMap, FILE *out) {
    fprintf(out, " \"primary\": [");
    bool firstBubble = TRUE;
    for (uint64_t i = 0; i < gF->length; i++) {

        // bubble and hap info
        Bubble *b = &bg->bubbles[gF->refStart + i];
        RleString *hap1 = b->alleles[gF->haplotypeString1[i]];
        RleString *hap2 = b->alleles[gF->haplotypeString2[i]];

        // we only care about hets
        if (hap1 == hap2) continue;

        if (firstBubble) {
            fprintf(out, "\n  {\n");
            firstBubble = FALSE;
        } else {
            fprintf(out, ",\n  {\n");
        }

        // bubble info
        int64_t trueRefStartPos = bamChunk->chunkBoundaryStart + reference_rleToNonRleCoordMap[b->refStart];
        double strandSkew = bubble_phasedStrandSkew(b, readsToPSeqs, gF);
        fprintf(out, "   \"refPos\": %"PRId64",\n", trueRefStartPos);
        fprintf(out, "   \"strandSkew\": %f,\n", strandSkew);

        // read info
        int64_t hap1AlleleNo = gF->haplotypeString1[i];
        int64_t hap2AlleleNo = gF->haplotypeString2[i];
        fprintf(out, "   \"reads\": [");
        for (uint64_t j = 0; j < b->readNo; j++) {
            if (j != 0) fprintf(out, ",");
            fprintf(out, "\n    {\n");

            BamChunkReadSubstring *bcrSubstring = b->reads[j];
            stProfileSeq *pSeq = stHash_search(readsToPSeqs, bcrSubstring->read);

            // read info
            char* readName = bcrSubstring->read->readName;
            double qualVal = bcrSubstring->qualValue;
            fprintf(out, "     \"name\": \"%s\",\n", readName);
            fprintf(out, "     \"qual\": %f,\n", qualVal);

            bool inHap1 = FALSE;
            bool inHap2 = FALSE;
            double readHap1Support = b->alleleReadSupports[hap1AlleleNo * b->readNo + j];
            double readHap2Support = b->alleleReadSupports[hap2AlleleNo * b->readNo + j];

            if (stSet_search(gF->reads1, pSeq) != NULL) {
                inHap1 = TRUE;
            }

            if (stSet_search(gF->reads2, pSeq) != NULL) {
                inHap2 = TRUE;
            }

            if (inHap1 == inHap2) {
                fprintf(out, "     \"hap\": %d,\n", 0);
            } else {
                fprintf(out, "     \"hap\": %d,\n", inHap1 ? 1 : 2);
            }
            fprintf(out, "     \"hapSupportH1\": %f,\n", readHap1Support);
            fprintf(out, "     \"hapSupportH2\": %f\n", readHap2Support);
            fprintf(out, "    }");
        }
        fprintf(out, "\n   ]");
        fprintf(out, "\n  }");
    }
    fprintf(out, "\n ]");
}

stSet *filterReadsByCoverageDepth2(stList *profileSeqs, Params *params) {
    stList *filteredProfileSeqs = stList_construct();
    stList *discardedProfileSeqs = stList_construct();
    filterReadsByCoverageDepth(profileSeqs, params->phaseParams, filteredProfileSeqs, discardedProfileSeqs);
    stSet *discardedReadsSet = stList_getSet(discardedProfileSeqs);
    stList_setDestructor(filteredProfileSeqs, NULL);
    stList_setDestructor(discardedProfileSeqs, NULL);
    stList_destruct(filteredProfileSeqs);
    stList_destruct(discardedProfileSeqs);

    return discardedReadsSet;
}

stGenomeFragment *bubbleGraph_phaseBubbleGraph(BubbleGraph *bg, stReference *ref, stList *reads, Params *params,
                                               stHash **readsToPSeqs) {
    /*
     * Runs phasing algorithm to split the reads embedded in the bubble graph into two partitions.
     *
     * Splits the forward and reverse strands to phase separately. After phasing them separately
     * joins them into one hmm.
     */

    // for logging
    char *logIdentifier = getLogIdentifier();

    // Generate profile sequences
    assert(ref->length == bg->bubbleNo);
    *readsToPSeqs = bubbleGraph_getProfileSeqs(bg, ref);
    stList *profileSeqs = stHash_getValues(*readsToPSeqs);

    assert(stList_length(reads) >= stList_length(profileSeqs));
    if (stList_length(reads) != stList_length(profileSeqs)) {
        st_logInfo(" %s In converting from reads to profile sequences have %" PRIi64 " reads and %" PRIi64 " profile sequences\n",
                logIdentifier, stList_length(reads), stList_length(profileSeqs));
    }

    // Remove excess coverage reads
    // Filter reads so that the maximum coverage depth does not exceed params->maxCoverageDepth
    st_logInfo(" %s Filtering reads by coverage depth\n", logIdentifier);
    stSet *discardedReadsSet = filterReadsByCoverageDepth2(profileSeqs, params);

    // Partition reads based upon strand
    st_logInfo(" %s Partitioning reads by strand for phasing\n", logIdentifier);
    stList *forwardStrandProfileSeqs = stList_construct();
    stList *reverseStrandProfileSeqs = stList_construct();
    for (int64_t i = 0; i < stList_length(reads); i++) {
        BamChunkRead *r = stList_get(reads, i);
        stProfileSeq *pSeq = stHash_search(*readsToPSeqs, r);
        if (pSeq != NULL &&
            stSet_search(discardedReadsSet, pSeq) == NULL) { // Has a pSeq and is not one of the filtered reads
            if (r->forwardStrand) {
                stList_append(forwardStrandProfileSeqs, pSeq);
            } else {
                stList_append(reverseStrandProfileSeqs, pSeq);
            }
        }
    }
    st_logInfo(" %s Got %" PRIi64 " forward strand reads for phasing and %" PRIi64 " negative strand reads for phasing\n",
               logIdentifier, stList_length(forwardStrandProfileSeqs), stList_length(reverseStrandProfileSeqs));

    // Deal with the case that the alignment is empty
    if (stList_length(profileSeqs) == 0) {
        stGenomeFragment *gf = stGenomeFragment_constructEmpty(ref, 0, 0, stSet_construct(), stSet_construct());
        stList_destruct(profileSeqs);
        stList_destruct(forwardStrandProfileSeqs);
        stList_destruct(reverseStrandProfileSeqs);
        stSet_destruct(discardedReadsSet);
        free(logIdentifier);
        return gf;
    }

    // Run phasing for each strand partition
    stRPHmmParameters *phaseParamsCopy = stRPHmmParameters_copy(params->phaseParams);
    phaseParamsCopy->includeAncestorSubProb = 0; // Switch off using ancestor substitution probabilities in calculating the hmm probs

    st_logInfo(" %s Phasing forward strand reads\n", logIdentifier);
    stList *tilingPathForward = getRPHmms(forwardStrandProfileSeqs, params->phaseParams);
    stList_setDestructor(tilingPathForward, NULL);

    st_logInfo(" %s Phasing reverse strand reads\n", logIdentifier);
    stList *tilingPathReverse = getRPHmms(reverseStrandProfileSeqs, params->phaseParams);
    stList_setDestructor(tilingPathReverse, NULL);

    // Join the hmms
    st_logInfo(" %s Joining forward and reverse strand phasing\n", logIdentifier);
    stRPHmm *hmm = fuseTilingPath(mergeTwoTilingPaths(tilingPathForward, tilingPathReverse));

    // Run the forward-backward algorithm
    phaseParamsCopy->includeAncestorSubProb = 1; // Now switch on using ancestor substitution probabilities in calculating the final, root hmm probs
    stRPHmm_forwardBackward(hmm);

    st_logInfo(" %s Forward probability of the hmm: %f, backward prob: %f\n", logIdentifier, (float) hmm->forwardLogProb,
               (float) hmm->backwardLogProb);

    // Now compute a high probability path through the hmm
    stList *path = stRPHmm_forwardTraceBack(hmm);

    assert(hmm->refStart >= 0);
    assert(hmm->refStart + hmm->refLength <= bg->bubbleNo);

    // Compute the genome fragment
    stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);

    // Refine the genome fragment by re-partitioning the reads iteratively
    stGenomeFragment_refineGenomeFragment(gF, hmm, path, params->phaseParams->roundsOfIterativeRefinement);

    // Sanity checks
    assert(gF->refStart >= 0);
    assert(gF->refStart + gF->length <= bg->bubbleNo);
    assert(gF->length == hmm->refLength);

    // For reads that exceeded the coverage depth, add them back to the haplotype they fit best
    stSetIterator *it = stSet_getIterator(discardedReadsSet);
    stProfileSeq *pSeq = NULL;
    while ((pSeq = stSet_getNext(it)) != NULL) {
        double i = getLogProbOfReadGivenHaplotype(gF->haplotypeString1, gF->refStart, gF->length, pSeq, gF->reference);
        double j = getLogProbOfReadGivenHaplotype(gF->haplotypeString2, gF->refStart, gF->length, pSeq, gF->reference);
        stSet_insert(i < j ? gF->reads2 : gF->reads1, pSeq);
    }
    stSet_destructIterator(it);

    // Set any homozygous alts back to being homozygous reference
    // This is really a hack because sometimes the phasing algorithm picks a non-reference allele for a homozygous
    // position
    /*for (uint64_t i = 0; i < gF->length; i++) {
        Bubble *b = &bg->bubbles[gF->refStart + i];

        if (gF->haplotypeString1[i] == gF->haplotypeString2[i]) {
            //|| binomialPValue(gF->readsSupportingHaplotype1[i] + gF->readsSupportingHaplotype2[i], gF->readsSupportingHaplotype1[i]) < 0.05) { // gF->readsSupportingHaplotype1[i] < 5 || gF->readsSupportingHaplotype2[i] < 5) { // is homozygous
            int64_t refAlleleIndex = bubble_getReferenceAlleleIndex(b);
            if (refAlleleIndex != -1) { // is homozygous alt
                gF->haplotypeString1[i] = refAlleleIndex; // set to reference allele
                gF->haplotypeString2[i] = refAlleleIndex;
            }
        }
    }*/

    // Check / log the result
    bubbleGraph_logPhasedBubbleGraph(bg, hmm, path, *readsToPSeqs, profileSeqs, gF);

    // Set destructors for later cleanup of profile sequences
    stSet_setDestructor(gF->reads1, (void (*)(void *)) stProfileSeq_destruct);
    stSet_setDestructor(gF->reads2, (void (*)(void *)) stProfileSeq_destruct);
    assert(stList_length(profileSeqs) == stSet_size(gF->reads1) + stSet_size(gF->reads2));
    assert(stSet_sizeOfIntersection(gF->reads1, gF->reads2) == 0);

    // Cleanup
    stRPHmmParameters_destruct(phaseParamsCopy);
    stSet_destruct(discardedReadsSet);
    stList_destruct(forwardStrandProfileSeqs);
    stList_destruct(reverseStrandProfileSeqs);
    stList_destruct(profileSeqs);
    stRPHmm_destruct(hmm, true);
    stList_destruct(path);
    free(logIdentifier);

    return gF;
}

Poa *bubbleGraph_getNewPoa(BubbleGraph *bg, uint64_t *consensusPath, Poa *poa, stList *reads, Params *params) {

    // Get new consensus string
    int64_t *poaToConsensusMap;
    RleString *newConsensusString = bubbleGraph_getConsensusString(bg, consensusPath, &poaToConsensusMap,
                                                                   params->polishParams);

    // Get anchor alignments
    stList *anchorAlignments = poa_getAnchorAlignments(poa, poaToConsensusMap, stList_length(reads),
                                                       params->polishParams);

    // Generated updated poa
    Poa *poa2 = poa_realign(reads, anchorAlignments, newConsensusString, params->polishParams);

    // Cleanup
    free(poaToConsensusMap);
    rleString_destruct(newConsensusString);
    stList_destruct(anchorAlignments);

    return poa2;
}

/*
 * Stuff to manage allele-strand-skew
 */

void bubble_calculateStrandSkews(Bubble *b, double *skews) {
    // Calculate the strand specific read supports
    double forwardStrandSupports[b->alleleNo];
    double reverseStrandSupports[b->alleleNo];
    uint64_t totalForward = 0, totalBackward = 0;
    for (int64_t j = 0; j < b->alleleNo; j++) {
        forwardStrandSupports[j] = 0.0;
        reverseStrandSupports[j] = 0.0;
    }
    for (int64_t i = 0; i < b->readNo; i++) {
        BamChunkReadSubstring *r = b->reads[i];
        double *d;
        if (r->read->forwardStrand) {
            totalForward++;
            d = forwardStrandSupports;
        } else {
            totalBackward++;
            d = reverseStrandSupports;
        }
        for (int64_t j = 0; j < b->alleleNo; j++) {
            d[j] += b->alleleReadSupports[j * b->readNo + i];
        }
    }

    // Calculate the average allele skew
    for (int64_t j = 0; j < b->alleleNo; j++) {
        skews[j] = (forwardStrandSupports[j] / totalForward - reverseStrandSupports[j] / totalBackward) /
                   (fabs(forwardStrandSupports[j] + reverseStrandSupports[j]) / (totalForward + totalBackward));
    }
}

uint128_t bionomialCoefficient(int64_t n, int64_t k) {
    uint128_t ans = 1;
    k = k > n - k ? n - k : k;
    int64_t j = 1;
    for (; j <= k; j++, n--) {
        if (n % j == 0) {
            ans *= n / j;
        } else if (ans % j == 0) {
            ans = ans / j * n;
        } else {
            ans = (ans * n) / j;
        }
    }
    return ans;
}

double binomialPValue(int64_t n, int64_t k) {
    uint128_t j = 0.0;
    k = k < n / 2 ? n - k : k;
    for (int64_t i = k; i <= n; i++) {
        j += bionomialCoefficient(n, i);
    }
    return j / pow(2.0, n);
}

double bubble_phasedStrandSkew(Bubble *b, stHash *readsToPSeqs, stGenomeFragment *gf) {
    int64_t reads = 0, positives = 0;
    for (int64_t i = 0; i < b->readNo; i++) {
        stProfileSeq *pSeq = stHash_search(readsToPSeqs, b->reads[i]->read);
        assert(pSeq != NULL);
        if (stSet_search(gf->reads1, pSeq) != NULL) {
            reads++;
            if (b->reads[i]->read->forwardStrand) {
                positives++;
            }
        } else if (stSet_search(gf->reads2, pSeq) != NULL) {
            reads++;
            if (!b->reads[i]->read->forwardStrand) {
                positives++;
            }
        }
    }
    return binomialPValue(reads, positives);
}

double bubbleGraph_skewedBubbles(BubbleGraph *bg, stHash *readsToPSeqs, stGenomeFragment *gf) {
    int64_t skewedBubbles = 0;
    for (int64_t i = 0; i < bg->bubbleNo; i++) {
        skewedBubbles += bubble_phasedStrandSkew(&bg->bubbles[i], readsToPSeqs, gf) < 0.05 ? 1 : 0;
    }
    return ((float) skewedBubbles) / bg->bubbleNo;
}

