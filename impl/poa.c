/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"
#include <sonLibListPrivate.h>
#include <htsIntegration.h>

char *getLogIdentifier() {
    char *logIdentifier;
# ifdef _OPENMP
    logIdentifier = stString_print(" T%02d", omp_get_thread_num());
# else
    logIdentifier = stString_copy("");
# endif
    return logIdentifier;
}

PoaBaseObservation *poaBaseObservation_construct(int64_t readNo, int64_t offset, double weight) {
    PoaBaseObservation *poaBaseObservation = st_calloc(1, sizeof(PoaBaseObservation));

    poaBaseObservation->readNo = readNo;
    poaBaseObservation->offset = offset;
    poaBaseObservation->weight = weight;

    return poaBaseObservation;
}

void poaBaseObservation_destruct(PoaBaseObservation *poaBaseObservation) {
    free(poaBaseObservation);
}

PoaInsert *poaInsert_construct(RleString *insert, double weight, bool strand) {
    PoaInsert *poaInsert = st_calloc(1, sizeof(PoaInsert));
    poaInsert->observations = stList_construct3(0, (void (*)(void *)) poaBaseObservation_destruct);

    poaInsert->insert = insert;
    if (strand) {
        poaInsert->weightForwardStrand = weight;
    } else {
        poaInsert->weightReverseStrand = weight;
    }

    return poaInsert;
}

void poaInsert_destruct(PoaInsert *poaInsert) {
    stList_destruct(poaInsert->observations);
    rleString_destruct(poaInsert->insert);
    free(poaInsert);
}

double poaInsert_getWeight(PoaInsert *insert) {
    return insert->weightForwardStrand + insert->weightReverseStrand;
    //return isBalanced(insert->weightForwardStrand, insert->weightReverseStrand, 100) ? insert->weightForwardStrand + insert->weightReverseStrand : 0.0;
}

PoaDelete *poaDelete_construct(int64_t length, double weight, bool strand) {
    PoaDelete *poaDelete = st_calloc(1, sizeof(PoaDelete));
    poaDelete->observations = stList_construct3(0, (void (*)(void *)) poaBaseObservation_destruct);

    poaDelete->length = length;
    if (strand) {
        poaDelete->weightForwardStrand = weight;
    } else {
        poaDelete->weightReverseStrand = weight;
    }

    return poaDelete;
}

void poaDelete_destruct(PoaDelete *poaDelete) {
    stList_destruct(poaDelete->observations);
    free(poaDelete);
}

double poaDelete_getWeight(PoaDelete *delete) {
    return delete->weightForwardStrand + delete->weightReverseStrand;
    //return isBalanced(delete->weightForwardStrand, delete->weightReverseStrand, 100) ? delete->weightForwardStrand + delete->weightReverseStrand : 0.0;
}

PoaNode *poaNode_construct(Poa *poa, char base, uint64_t repeatCount) {
    PoaNode *poaNode = st_calloc(1, sizeof(PoaNode));

    // for when people put crazy characters in their reference
    if (poa->alphabet->convertSymbolToChar(poa->alphabet->convertCharToSymbol(base)) == 'N') {
        base = 'N';
    }

    poaNode->inserts = stList_construct3(0, (void (*)(void *)) poaInsert_destruct);
    poaNode->deletes = stList_construct3(0, (void (*)(void *)) poaDelete_destruct);
    poaNode->base = base;
    poaNode->repeatCount = repeatCount;
    poaNode->baseWeights = st_calloc(poa->alphabet->alphabetSize, sizeof(double)); // Encoded using Symbol enum
    poaNode->repeatCountWeights = st_calloc(poa->maxRepeatCount, sizeof(double));
    poaNode->observations = stList_construct3(0, (void (*)(void *)) poaBaseObservation_destruct);

    return poaNode;
}

void poaNode_destruct(PoaNode *poaNode) {
    stList_destruct(poaNode->inserts);
    stList_destruct(poaNode->deletes);
    stList_destruct(poaNode->observations);
    free(poaNode->baseWeights);
    free(poaNode->repeatCountWeights);
    free(poaNode);
}

Poa *poa_getReferenceGraph(RleString *reference, Alphabet *alphabet, uint64_t maxRepeatCount) {
    Poa *poa = st_calloc(1, sizeof(Poa));

    poa->alphabet = alphabet;
    poa->maxRepeatCount = maxRepeatCount;
    poa->nodes = stList_construct3(0, (void (*)(void *)) poaNode_destruct);
    poa->refString = rleString_copy(reference);

    stList_append(poa->nodes, poaNode_construct(poa, 'N', 1)); // Add empty prefix node
    for (int64_t i = 0; i < reference->length; i++) {
        stList_append(poa->nodes,
                poaNode_construct(poa, (char) toupper(reference->rleString[i]), reference->repeatCounts[i]));
    }

    return poa;
}

void poa_destruct(Poa *poa) {
    rleString_destruct(poa->refString);
    stList_destruct(poa->nodes);
    free(poa);
}

int cmpAlignedPairsByCoordinates(const void *a, const void *b) {
    /*
     * Compares aligned pairs, represented as stIntTuples of the form (weight, x, y) first by
     * x coordinate and then by y coordinate, ignoring the weight.
     */
    stIntTuple *one = (stIntTuple *) a, *two = (stIntTuple *) b;

    int i = stIntTuple_get(one, 1) < stIntTuple_get(two, 1) ? -1 : stIntTuple_get(one, 1) > stIntTuple_get(two, 1) ? 1
                                                                                                                   : 0;
    if (i == 0) {
        i = stIntTuple_get(one, 2) < stIntTuple_get(two, 2) ? -1 : stIntTuple_get(one, 2) > stIntTuple_get(two, 2) ? 1
                                                                                                                   : 0;
    }

    return i;
}

int cmpAlignedPairsByInvertedCoordinates(const void *a, const void *b) {
    /*
     * As cmpAlignedPairsByCoordinates, but compares by y coordinate and then x coordinate.
     */
    stIntTuple *one = (stIntTuple *) a, *two = (stIntTuple *) b;

    int i = stIntTuple_get(one, 2) < stIntTuple_get(two, 2) ? -1 : stIntTuple_get(one, 2) > stIntTuple_get(two, 2) ? 1
                                                                                                                   : 0;
    if (i == 0) {
        i = stIntTuple_get(one, 1) < stIntTuple_get(two, 1) ? -1 : stIntTuple_get(one, 1) > stIntTuple_get(two, 1) ? 1
                                                                                                                   : 0;
    }

    return i;
}

bool isMatch(stSortedSet *matchesSet, int64_t x, int64_t y) {
    stIntTuple pair[4];
    pair[0] = 3;
    pair[2] = x;
    pair[3] = y;
    return stSortedSet_search(matchesSet, &pair) != NULL;
}

static void
addToInserts(PoaNode *node, RleString *insert, double weight, bool strand, PoaBaseObservation *observation) {
    /*
     * Add given insert to node.
     */

    PoaInsert *poaInsert = NULL;
    // Check if the complete insert is already in the poa graph:
    for (int64_t m = 0; m < stList_length(node->inserts); m++) {
        PoaInsert *tmp = stList_get(node->inserts, m);
        if (rleString_eq(tmp->insert, insert)) {
            poaInsert = tmp;
            break;
        }
    }
    // otherwise create and save it
    if (poaInsert == NULL) {
        poaInsert = poaInsert_construct(rleString_copy(insert), 0, FALSE);
        stList_append(node->inserts, poaInsert);
    }

    // update with (stranded) weight and observation
    if (strand) {
        poaInsert->weightForwardStrand += weight;
    } else {
        poaInsert->weightReverseStrand += weight;
    }
    stList_append(poaInsert->observations, observation);
}

static void addToDeletes(PoaNode *node, int64_t length, double weight, bool strand, PoaBaseObservation *observation) {
    /*
     * Add given deletion to node.
     */

    PoaDelete *poaDelete = NULL;
    // Check if the delete is already in the poa graph:
    for (int64_t m = 0; m < stList_length(node->deletes); m++) {
        PoaDelete *tmp = stList_get(node->deletes, m);
        if (tmp->length == length) {
            poaDelete = tmp;
            break;
        }
    }
    // otherwise create and save it
    if (poaDelete == NULL) {
        poaDelete = poaDelete_construct(length, 0, FALSE);
        stList_append(node->deletes, poaDelete);
    }

    // update with (stranded) weight and observation
    if (strand) {
        poaDelete->weightForwardStrand += weight;
    } else {
        poaDelete->weightReverseStrand += weight;
    }
    stList_append(poaDelete->observations, observation);
}

static bool matchesReferenceSubstring(RleString *refString, int64_t refStart, RleString *str, int64_t length,
                                      bool compareRepeatCounts) {
    /*
     * Returns true if the given string str matches the given reference substring starting
     * from the given reference position, refStart. Optionally also compares the repeat lengths for equality also.
     */
    for (int64_t l = 0; l < length; l++) {
        if (refString->rleString[refStart + l] != str->rleString[l] ||
            (compareRepeatCounts && refString->repeatCounts[refStart + l] != str->repeatCounts[l])) {
            return 0;
        }
    }
    return 1;
}

static bool hasInternalRepeat(RleString *str, int64_t repeatLength, bool compareRepeatCounts) {
    /*
     * Establishes if str has an internal repeat of length repeatLength.
     * e.g. if ATATAT, internal repeat AT (length 2) is an internal repeat, but ATA is not.
     */
    if (str->length % repeatLength != 0) { // If not divisible by repeatLength then can not be repeat
        return 0;
    }
    for (int64_t i = repeatLength; i < str->length; i += repeatLength) {
        for (int64_t j = 0; j < repeatLength; j++) {
            if (str->rleString[j] != str->rleString[j + i] ||
                (compareRepeatCounts && str->repeatCounts[j] != str->repeatCounts[j + i])) {
                return 0;
            }
        }
    }
    return 1;
}

int64_t getShift(RleString *refString, int64_t refStart, RleString *str, bool compareRepeatCounts) {
    // Walk back over reference sequence and see if indel can be shifted

    // Establish minimal internal repeat length
    // if ATATAT, minimal internal repeat is AT,
    // similarly if AAAAAAA then minimal internal repeat is A
    int64_t minRepeatLength = 0;
    while (minRepeatLength++ < str->length) {
        if (hasInternalRepeat(str, minRepeatLength, compareRepeatCounts)) {
            break;
        }
    }

    // Now walk back by multiples of minimal internal repeat length
    for (int64_t k = refStart - minRepeatLength; k >= 0; k -= minRepeatLength) {
        if (!matchesReferenceSubstring(refString, k, str, minRepeatLength, compareRepeatCounts)) {
            break;
        }
        refStart = k;
    }

    // Deal with special case that insert is of length 1 and compareRepeatsCounts=True, in which case further shift maybe possible
    // because repeat counts need not be equal for single base insert to be shifted back one base
    if (str->length == 1 && compareRepeatCounts && refStart > 0 &&
        refString->rleString[refStart - 1] == str->rleString[0]) {
        refStart--;
    }

    return refStart;
}

int64_t getMaxCommonSuffixLength(RleString *str1, int64_t length1, RleString *str2, bool compareRepeatCounts) {
    /*
     * Returns the length of the maximum suffix of the reference string ending at refStart (inclusive)
     * that is the same as a suffix of str.
     */
    int64_t i = 0;
    while (length1 - i - 1 >= 0 && (int64_t) str2->length - i - 1 >= 0) {
        if (str1->rleString[length1 - 1 - i] != str2->rleString[str2->length - 1 - i] ||
            (compareRepeatCounts && str1->repeatCounts[length1 - 1 - i] != str2->repeatCounts[str2->length - 1 - i])) {
            break;
        }
        i++;
    }

    return i;
}

void poa_augment(Poa *poa, RleString *read, bool readStrand, int64_t readNo, stList *matches, stList *inserts,
                 stList *deletes,
                 PolishParams *polishParams) {
    // Add weights of matches to the POA graph

    // For each match in alignment subgraph identify its corresponding node in the POA graph
    // add the weight of the match to the POA node
    for (int64_t i = 0; i < stList_length(matches); i++) {
        stIntTuple *match = stList_get(matches, i);

        PoaNode *node = stList_get(poa->nodes, stIntTuple_get(match, 1) + 1); // Get corresponding POA node

        // Add base weight to POA node
        int64_t j = stIntTuple_get(match, 2), weight = stIntTuple_get(match, 0);
        assert(poa->alphabet->convertCharToSymbol(read->rleString[j]) < poa->alphabet->alphabetSize);
        node->baseWeights[poa->alphabet->convertCharToSymbol(read->rleString[j])] += weight;
        assert(read->repeatCounts[j] >= 0);
        int64_t rc = read->repeatCounts[j] < poa->maxRepeatCount ? read->repeatCounts[j] : poa->maxRepeatCount - 1;
        node->repeatCountWeights[rc] += weight;

        // PoaObservation
        stList_append(node->observations,
                      poaBaseObservation_construct(readNo, stIntTuple_get(match, 2), stIntTuple_get(match, 0)));
    }

    // Create a set of match coordinates

    stSortedSet *matchesSet = stSortedSet_construct3(cmpAlignedPairsByCoordinates, NULL);
    for (int64_t i = 0; i < stList_length(matches); i++) {
        stIntTuple *match = stList_get(matches, i);

        assert(stSortedSet_search(matchesSet, match) == NULL);
        stSortedSet_insert(matchesSet, match); // Add to matches
    }

    // Add inserts to the POA graph

    // Sort the inserts first by the reference coordinate and then by read coordinate
    stList_sort(inserts, cmpAlignedPairsByCoordinates);

    // Let a complete-insert be a sequence of inserts with the same reference coordinate i
    // and consecutive read coordinates j, j+1, ..., j+n, such that the i,j-1 is a match or equal to (-1,-1) (the beginning)
    // in the alignment subgraph and i+1,j+n+1 is similarly a match or equal to (N, M) (the end of the alignment).

    // Enumerate set of complete inserts
    for (int64_t i = 0; i < stList_length(inserts);) {

        stIntTuple *insertStart = stList_get(inserts, i); // Start of putative complete-insert

        int64_t j = i + 1;
        for (; j < stList_length(inserts); j++) {
            stIntTuple *insertEnd = stList_get(inserts, j); // End of putative complete-insert

            // If they don't have the same reference coordinate then not part of same complete-insert
            if (stIntTuple_get(insertStart, 1) != stIntTuple_get(insertEnd, 1)) {
                break;
            }

            // If they don't form a contiguous sequence of read coordinates then not part of same complete-insert
            if (stIntTuple_get(insertStart, 2) + j - i != stIntTuple_get(insertEnd, 2)) {
                break;
            }
        }

        // At this point i (inclusive) and j (exclusive) form the start and end of a maximal putative
        // complete insert sequence in inserts

        // Now enumerate complete-inserts in this interval

        for (int64_t k = i; k < j; k++) {

            // If k position is not flanked by a preceding match or the beginning then can not be a complete insert
            if (!isMatch(matchesSet, stIntTuple_get(insertStart, 1), stIntTuple_get(insertStart, 2) + k - i - 1) &&
                stIntTuple_get(insertStart, 2) + k - i - 1 > -1) {
                continue;
            }

            for (int64_t l = k; l < j; l++) {

                // If l position is not flanked by a proceeding match or the end then can not be a complete insert
                if (!isMatch(matchesSet, stIntTuple_get(insertStart, 1) + 1,
                             stIntTuple_get(insertStart, 2) + l - i + 1) &&
                    stIntTuple_get(insertStart, 2) + l - i + 1 < read->length) {
                    continue;
                }

                // At this point k (inclusive) and l (inclusive) represent a complete-insert

                // Calculate weight and label, including repeat counts
                assert(stIntTuple_get(insertStart, 2) + k - i == stIntTuple_get(stList_get(inserts, k), 2));
                assert(stIntTuple_get(stList_get(inserts, k), 2) >= 0);

                RleString *insert = rleString_copySubstring(read, stIntTuple_get(stList_get(inserts, k), 2), l + 1 - k);
                double insertWeight = UINT_MAX;
                for (int64_t m = k; m < l + 1; m++) {
                    stIntTuple *insertBase = stList_get(inserts, m);
                    insertWeight =
                            insertWeight < stIntTuple_get(insertBase, 0) ? insertWeight : stIntTuple_get(insertBase, 0);
                }

                // Get the leftmost node in the poa graph to which the insert will connect

                // First find the left point to which the insert will be connected
                assert(stIntTuple_get(insertStart, 1) >= -1);
                int64_t insertPosition = stIntTuple_get(insertStart, 1) + 1;

                // Now walk back over reference sequence and see if insert can be left-shifted
                insertPosition = getShift(poa->refString, insertPosition, insert,
                                          polishParams->poaConstructCompareRepeatCounts);
                assert(insertPosition >= 0 && insertPosition <= poa->refString->length);

                // Finally see if can be shifted by common suffix
                int64_t commonSuffixLength = getMaxCommonSuffixLength(poa->refString, insertPosition, insert,
                                                                      polishParams->poaConstructCompareRepeatCounts);
                if (commonSuffixLength > 0) {
                    rleString_rotateString(insert, commonSuffixLength, polishParams->useRunLengthEncoding);
                    insertPosition -= commonSuffixLength;
                }
                assert(insertPosition >= 0);

                // Add insert to graph at leftmost position
                addToInserts(stList_get(poa->nodes, insertPosition), insert, insertWeight, readStrand,
                             poaBaseObservation_construct(readNo, stIntTuple_get(stList_get(inserts, k), 2),
                                                          insertWeight));

                // Cleanup
                rleString_destruct(insert);
            }
        }

        // Increase i to start of next maximal complete-insert
        i = j;
    }

    // Add deletes to the POA graph

    // Sort the deletes first by the read coordinate and then by reference coordinate
    stList_sort(deletes, cmpAlignedPairsByInvertedCoordinates);

    // Analogous to a complete-insert, let a complete-delete be a sequence of deletes with the same read coordinate j
    // and consecutive reference coordinates i, i+1, ..., i+m, such that the i-1,j is a match or equal to (-1,-1) (the beginning)
    // in the alignment subgraph and i+m+1,j+1 is similarly a match or (N, M) (the alignment end).

    // Enumerate set of complete-deletes, adding them to the graph
    for (int64_t i = 0; i < stList_length(deletes);) {
        stIntTuple *deleteStart = stList_get(deletes, i); // Start of putative complete-delete

        int64_t j = i + 1;
        for (; j < stList_length(deletes); j++) {
            stIntTuple *deleteEnd = stList_get(deletes, j); // End of putative complete-delete

            // If they don't have the same read coordinate then not part of same complete-insert
            if (stIntTuple_get(deleteStart, 2) != stIntTuple_get(deleteEnd, 2)) {
                break;
            }

            // If they don't form a contiguous sequence of read coordinates then not part of same complete-insert
            if (stIntTuple_get(deleteStart, 1) + j - i != stIntTuple_get(deleteEnd, 1)) {
                break;
            }
        }

        // At this point i (inclusive) and j (exclusive) form the start and end of a putative maximal
        // complete-delete sequence in deletes

        // Now enumerate complete-deletes in this interval

        for (int64_t k = i; k < j; k++) {

            // If k position is not flanked by a preceding match or the alignment beginning then can not be a complete-delete
            if (!isMatch(matchesSet, stIntTuple_get(deleteStart, 1) + k - i - 1, stIntTuple_get(deleteStart, 2)) &&
                stIntTuple_get(deleteStart, 1) + k - i - 1 > -1) {
                continue;
            }

            for (int64_t l = k; l < j; l++) {

                // If l position is not flanked by a proceeding match or the alignment end then can not be a complete-delete
                if (!isMatch(matchesSet, stIntTuple_get(deleteStart, 1) + l - i + 1,
                             stIntTuple_get(deleteStart, 2) + 1) &&
                    stIntTuple_get(deleteStart, 1) + l - i + 1 < poa->refString->length) {
                    continue;
                }

                // At this point k (inclusive) and l (inclusive) represent a complete-delete

                // Delete length
                int64_t deleteLength = l - k + 1;

                // Calculate weight
                double deleteWeight = UINT_MAX;
                for (int64_t m = k; m < l + 1; m++) {
                    stIntTuple *delete = stList_get(deletes, m);
                    deleteWeight = deleteWeight < stIntTuple_get(delete, 0) ? deleteWeight : stIntTuple_get(delete, 0);
                }

                // Get the leftmost node in the poa graph to which the delete will connect

                // First find the left point to which the delete would be connected
                assert(stIntTuple_get(deleteStart, 1) + k - i >= 0);
                int64_t deletePosition = stIntTuple_get(deleteStart, 1) + k - i;

                // Get string being deleted
                RleString *delete = rleString_copySubstring(poa->refString, deletePosition, deleteLength);

                // Now walk back over reference sequence and see if insert can be left-shifted
                deletePosition = getShift(poa->refString, deletePosition, delete,
                                          polishParams->poaConstructCompareRepeatCounts);

                // Finally see if can be shifted by common suffix
                deletePosition -= getMaxCommonSuffixLength(poa->refString, deletePosition, delete,
                                                           polishParams->poaConstructCompareRepeatCounts);
                rleString_destruct(delete);

                // Add delete to graph at leftmost position
                addToDeletes(stList_get(poa->nodes, deletePosition), deleteLength, deleteWeight, readStrand,
                             poaBaseObservation_construct(readNo, stIntTuple_get(deleteStart, 2), deleteWeight));
            }
        }

        // Increase i to start of next maximal complete-delete
        i = j;
    }

    // Cleanup
    stSortedSet_destruct(matchesSet);
}

stList *poa_getAnchorAlignments(Poa *poa, const int64_t *poaToConsensusMap, int64_t noOfReads, PolishParams *pp) {

    // Allocate anchor alignments
    stList *anchorAlignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
    for (int64_t i = 0; i < noOfReads; i++) {
        stList_append(anchorAlignments, stList_construct3(0, (void (*)(void *)) stIntTuple_destruct));
    }

    // Walk through the weights of the POA to construct the anchor alignments
    for (int64_t i = 1; i < stList_length(poa->nodes); i++) {
        PoaNode *poaNode = stList_get(poa->nodes, i);
        // If consensus index is null, then just use the underlying reference sequence of the POA.
        int64_t consensusIndex = poaToConsensusMap == NULL ? i - 1 : poaToConsensusMap[i - 1];
        if (consensusIndex != -1) { // Poa reference position is aligned to the consensus
            for (int64_t j = 0; j < stList_length(poaNode->observations); j++) {
                PoaBaseObservation *obs = stList_get(poaNode->observations, j);
                double normalizedObsWeight = obs->weight / PAIR_ALIGNMENT_PROB_1;

                if (normalizedObsWeight > pp->minPosteriorProbForAlignmentAnchors[0]) { // High confidence anchor pair
                    stList *anchorPairs = stList_get(anchorAlignments, obs->readNo);

                    // Figure out the exact diagonal expansion
                    int64_t expansion = (int64_t) pp->minPosteriorProbForAlignmentAnchors[1];
                    for (int64_t k = 2; k < pp->minPosteriorProbForAlignmentAnchorsLength; k += 2) {
                        if (normalizedObsWeight >= pp->minPosteriorProbForAlignmentAnchors[k]) {
                            expansion = (int64_t) pp->minPosteriorProbForAlignmentAnchors[k + 1];
                        } else {
                            break;
                        }
                    }

                    // The following is masking an underlying bug that allows for multiple high confidence
                    // alignments per position
                    if (stList_length(anchorPairs) == 0) {
                        stList_append(anchorPairs, stIntTuple_construct3(consensusIndex, obs->offset, expansion));
                    } else {
                        stIntTuple *pPair = stList_peek(anchorPairs);

                        if (stIntTuple_get(pPair, 0) < consensusIndex && stIntTuple_get(pPair, 1) < obs->offset) {
                            stList_append(anchorPairs, stIntTuple_construct3(consensusIndex, obs->offset, expansion));
                        }
                        //else {
                        //	fprintf(stderr, "Ooops read: %i x1: %i y1: %i x2: %i y2: %i %f\n", (int)obs->readNo,
                        //			(int)stIntTuple_get(pPair, 0), (int)stIntTuple_get(pPair, 1),
                        //			(int)consensusIndex, (int)obs->offset, (float)obs->weight/PAIR_ALIGNMENT_PROB_1);
                        //}
                    }

                }
            }
        }
    }

    return anchorAlignments;
}

static void adjustAnchors(stList *anchorPairs, int64_t index, int64_t adjustment) {
    for (int64_t i = 0; i < stList_length(anchorPairs); i++) {
        stIntTuple *pair = stList_get(anchorPairs, i);
        ((int64_t *) pair)[index + 1] += adjustment;
    }
}

/*
 * Generates aligned pairs and indel probs, but first crops reference to only include sequence from first
 * to last anchor position.
 */
void getAlignedPairsWithIndelsCroppingReference(RleString *reference,
                                                RleString *read, bool readStrand, stList *anchorPairs,
                                                stList **matches, stList **inserts, stList **deletes,
                                                PolishParams *polishParams) {
    // Crop reference, to avoid long unaligned prefix and suffix
    // that generates a lot of delete pairs

    // Get cropping coordinates
    // TODO I think we may want to extend refStart and refEnd by the length of the read before and after the first and last aligned positions
    int64_t firstRefPosition, endRefPosition;
    if (stList_length(anchorPairs) > 0) {
        stIntTuple *fPair = stList_get(anchorPairs, 0);
        firstRefPosition = stIntTuple_get(fPair, 0) - stIntTuple_get(fPair, 1);
        firstRefPosition = firstRefPosition < 0 ? 0 : firstRefPosition;

        stIntTuple *lPair = stList_peek(anchorPairs);
        endRefPosition = 1 + stIntTuple_get(lPair, 0) + (read->length - stIntTuple_get(lPair, 1));
        endRefPosition = endRefPosition > reference->length ? reference->length : endRefPosition;
    } else {
        firstRefPosition = 0;
        endRefPosition = reference->length;
    }
    assert(firstRefPosition <= reference->length && firstRefPosition >= 0);
    assert(endRefPosition <= reference->length && endRefPosition >= 0);

    // Adjust anchor positions
    adjustAnchors(anchorPairs, 0, -firstRefPosition);

    // Get symbol strings
    uint64_t maxRL = polishParams->useRunLengthEncoding ? (uint64_t) polishParams->repeatSubMatrix->maximumRepeatLength
                                                        : 2;
    SymbolString sX = rleString_constructSymbolString(reference, firstRefPosition, endRefPosition - firstRefPosition,
                                                      polishParams->alphabet, polishParams->useRepeatCountsInAlignment,
                                                      maxRL);
    SymbolString sY = rleString_constructSymbolString(read, 0, read->length, polishParams->alphabet,
                                                      polishParams->useRepeatCountsInAlignment,
                                                      maxRL);

    // Get alignment
    getAlignedPairsWithIndelsUsingAnchors(readStrand ? polishParams->stateMachineForForwardStrandRead :
                                          polishParams->stateMachineForReverseStrandRead, sX, sY,
                                          anchorPairs, polishParams->p, matches, deletes, inserts, 0, 0);

    // Cleanup symbol strings
    symbolString_destruct(sX);
    symbolString_destruct(sY);

    // Adjust back anchors
    adjustAnchors(anchorPairs, 0, firstRefPosition);

    // Shift matches/inserts/deletes
    adjustAnchors(*matches, 1, firstRefPosition);
    adjustAnchors(*inserts, 1, firstRefPosition);
    adjustAnchors(*deletes, 1, firstRefPosition);
}

Poa *poa_realign(stList *bamChunkReads, stList *anchorAlignments, RleString *reference,
                 PolishParams *polishParams) {
    // Build a reference graph with zero weights
    uint64_t maximumRepeatLength = 2; // MRL is exclusive
    if (polishParams->useRunLengthEncoding) {
        if (polishParams->repeatSubMatrix != NULL) {
            maximumRepeatLength = polishParams->repeatSubMatrix->maximumRepeatLength;
        } else {
            maximumRepeatLength = MAXIMUM_REPEAT_LENGTH;
        }
    }
    Poa *poa = poa_getReferenceGraph(reference, polishParams->alphabet, maximumRepeatLength);

    // For each read
    for (int64_t i = 0; i < stList_length(bamChunkReads); i++) {
        BamChunkRead *chunkRead = stList_get(bamChunkReads, i);

        // Generate set of posterior probabilities for matches, deletes and inserts with respect to reference.
        stList *matches = NULL, *inserts = NULL, *deletes = NULL;

        if(anchorAlignments == NULL) {
            SymbolString sX = rleString_constructSymbolString(reference, 0, reference->length, polishParams->alphabet,
                                                              polishParams->useRepeatCountsInAlignment, poa->maxRepeatCount - 1);
            SymbolString sY = rleString_constructSymbolString(chunkRead->rleRead, 0, chunkRead->rleRead->length,
                                                              polishParams->alphabet, polishParams->useRepeatCountsInAlignment, poa->maxRepeatCount - 1);

            getAlignedPairsWithIndels(chunkRead->forwardStrand ? polishParams->stateMachineForForwardStrandRead :
                                      polishParams->stateMachineForReverseStrandRead,
                                      sX, sY, polishParams->p, &matches, &deletes, &inserts, 0, 0);

            symbolString_destruct(sX);
            symbolString_destruct(sY);
        }
        else {
            getAlignedPairsWithIndelsCroppingReference(reference, chunkRead->rleRead, chunkRead->forwardStrand, stList_get(anchorAlignments, i),
                                                       &matches, &inserts, &deletes, polishParams);
        }

        // Add weights, edges and nodes to the poa
        poa_augment(poa, chunkRead->rleRead, chunkRead->forwardStrand, i, matches, inserts, deletes, polishParams);

        // Cleanup
        stList_destruct(matches);
        stList_destruct(inserts);
        stList_destruct(deletes);
    }

    return poa;
}

Poa *poa_realignOnlyAnchorAlignments(stList *bamChunkReads, stList *anchorAlignments, RleString *reference,
        PolishParams *polishParams) {
    // Build a reference graph with zero weights
    uint64_t maximumRepeatLength = 2; // MRL is exclusive
    if (polishParams->useRunLengthEncoding) {
        if (polishParams->repeatSubMatrix != NULL) {
            maximumRepeatLength = polishParams->repeatSubMatrix->maximumRepeatLength;
        } else {
            maximumRepeatLength = MAXIMUM_REPEAT_LENGTH;
        }
    }
    Poa *poa = poa_getReferenceGraph(reference, polishParams->alphabet, maximumRepeatLength);

    // For each read
    for (int64_t i = 0; i < stList_length(bamChunkReads); i++) {
        BamChunkRead *chunkRead = stList_get(bamChunkReads, i);
        stList *anchorAlignment = stList_get(anchorAlignments, i);

        // Generate set of posterior probabilities for matches, deletes and inserts with respect to reference.
        stList *matches = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
        stList *inserts = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
        stList *deletes = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);


        stListIterator *alignmentItor = stList_getIterator(anchorAlignment);
        stIntTuple *currAlign = stList_getNext(alignmentItor);
        int64_t posRef = stIntTuple_get(currAlign, 0);
        int64_t posRead = stIntTuple_get(currAlign, 1);

        while (TRUE) {
            if (currAlign == NULL) break;
            int64_t currAlignPosRef = stIntTuple_get(currAlign, 0);
            int64_t currAlignPosRead = stIntTuple_get(currAlign, 1);
            // Read delete
            if (posRef < currAlignPosRef) {
                stList_append(deletes, stIntTuple_construct3(PAIR_ALIGNMENT_PROB_1, currAlignPosRef, currAlignPosRead-1));
                posRef++;
            }

            // Read insert
            else if (posRead < currAlignPosRead) {
                stList_append(inserts, stIntTuple_construct3(PAIR_ALIGNMENT_PROB_1, currAlignPosRef-1, currAlignPosRead));
                posRead++;
            }

            // match
            else if (posRef == currAlignPosRef && posRead == currAlignPosRead) {
                stList_append(matches, stIntTuple_construct3(PAIR_ALIGNMENT_PROB_1, currAlignPosRef, currAlignPosRead));
                posRef++;
                posRead++;
                currAlign = stList_getNext(alignmentItor);
            }

            // should never happen
            else {
                assert(FALSE);
            }
        }

        // Add weights, edges and nodes to the poa
        poa_augment(poa, chunkRead->rleRead, chunkRead->forwardStrand, i, matches, inserts, deletes, polishParams);

        // Cleanup
        stList_destruct(matches);
        stList_destruct(inserts);
        stList_destruct(deletes);
        stList_destructIterator(alignmentItor);
    }

    return poa;
}

/*
 * Functions to calculate weights of poa nodes
 */

double poa_getReferenceNodeTotalMatchWeight(Poa *poa) {
    double weight = 0.0;
    for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
        PoaNode *node = stList_get(poa->nodes, i);
        weight += node->baseWeights[poa->alphabet->convertCharToSymbol(node->base)];
    }
    return weight;
}

double poa_getReferenceNodeTotalDisagreementWeight(Poa *poa) {
    double weight = 0.0;
    for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
        PoaNode *node = stList_get(poa->nodes, i);
        int64_t refSymbol = poa->alphabet->convertCharToSymbol(node->base);
        for (int64_t j = 0; j < poa->alphabet->alphabetSize; j++) {
            if (j != refSymbol) {
                weight += node->baseWeights[j];
            }
        }
    }
    return weight;
}

double poa_getInsertTotalWeight(Poa *poa) {
    double weight = 0.0;
    for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
        PoaNode *node = stList_get(poa->nodes, i);
        for (int64_t j = 0; j < stList_length(node->inserts); j++) {
            PoaInsert *insert = stList_get(node->inserts, j);
            weight += poaInsert_getWeight(insert) * insert->insert->length;
        }
    }
    return weight;
}

double poa_getDeleteTotalWeight(Poa *poa) {
    double weight = 0.0;
    for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
        PoaNode *node = stList_get(poa->nodes, i);
        for (int64_t j = 0; j < stList_length(node->deletes); j++) {
            PoaDelete *delete = stList_get(node->deletes, j);
            weight += poaDelete_getWeight(delete) * delete->length;
        }
    }
    return weight;
}

double poa_getTotalErrorWeight(Poa *poa) {
    return poa_getDeleteTotalWeight(poa) + poa_getInsertTotalWeight(poa) +
           poa_getReferenceNodeTotalDisagreementWeight(poa);
}

double *poaNode_getStrandSpecificBaseWeights(PoaNode *node, stList *bamChunkReads,
                                             double *totalWeight, double *totalPositiveWeight,
                                             double *totalNegativeWeight, Alphabet *a, stSet *readsToInclude) {
    /*
     * Calculate strand specific base weights.
     */
    *totalWeight = 0.0;
    *totalPositiveWeight = 0.0;
    *totalNegativeWeight = 0.0;
    double *baseWeights = st_calloc(a->alphabetSize * 2, sizeof(double));
    for (int64_t i = 0; i < stList_length(node->observations); i++) {
        PoaBaseObservation *baseObs = stList_get(node->observations, i);
        BamChunkRead *read = stList_get(bamChunkReads, baseObs->readNo);
        if (readsToInclude == NULL || stSet_search(readsToInclude, read) != NULL) {
            *totalWeight += baseObs->weight;
            char base = read->rleRead->rleString[baseObs->offset];
            baseWeights[a->convertCharToSymbol(base) * 2 +
                        (read->forwardStrand ? POS_STRAND_IDX : NEG_STRAND_IDX)] += baseObs->weight;
            if (read->forwardStrand) {
                *totalPositiveWeight += baseObs->weight;
            } else {
                *totalNegativeWeight += baseObs->weight;
            }
        }
    }

    return baseWeights;
}

/*
 * Functions to print poa
 */

void poa_printRepeatCountsCSV(Poa *poa, FILE *fH, stList *bamChunkReads) {
    fprintf(fH, "REF_INDEX,REF_BASE");
    fprintf(fH, ",REPEAT_COUNT_OBSxN(READ_BASE,READ_STRAND,REPEAT_COUNT,WEIGHT)\n");

    // Print info for each base in reference in turn
    for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
        PoaNode *node = stList_get(poa->nodes, i);

        fprintf(fH, "%" PRIi64 ",%c", i, node->base);

        for (int64_t j = 0; j < stList_length(node->observations); j++) {
            PoaBaseObservation *obs = stList_get(node->observations, j);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, obs->readNo);
            int64_t repeatCount = bamChunkRead->rleRead->repeatCounts[obs->offset];
            char base = bamChunkRead->rleRead->rleString[obs->offset];
            fprintf(fH, ",%c%c%" PRIi64 ",%.3f", base, bamChunkRead->forwardStrand ? '+' : '-', repeatCount,
                    obs->weight / PAIR_ALIGNMENT_PROB_1);
        }

        fprintf(fH, "\n");
    }
}

void poa_printDOT(Poa *poa, FILE *fH, stList *bamChunkReads) {

    char *insertColor = "\"darkgreen\"";
    char *backboneColor = "\"blue\"";
    char *deleteColor = "\"purple\"";

    fprintf(fH, "digraph poa {\nrankdir=\"LR\";\n");

    for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
        PoaNode *node = stList_get(poa->nodes, i);
        double runLengths[50] = {0};
        double weight = 0.0;
        for (int64_t o = 0; o < stList_length(node->observations); o++) {
            PoaBaseObservation *obvs = stList_get(node->observations, o);
            // node weight
            weight += obvs->weight;
            // run length obvs
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, obvs->readNo);
            RleString *rleString = bamChunkRead->rleRead;
            // skip non-match rates
            if (rleString->rleString[obvs->offset] != node->base) continue;
            // save run length
            int64_t rl = rleString->repeatCounts[obvs->offset];
            if (rl > 50) rl = 50;
            runLengths[rl - 1] += obvs->weight;
        }
        weight /= PAIR_ALIGNMENT_PROB_1;

        // build run length label
        stList *labelStrings = stList_construct3(0, free);
        stList_append(labelStrings, stString_print("%"PRId64, i));
        for (int64_t r = 0; r < 50; r++) {
            if (runLengths[r] != 0) {
                stList_append(labelStrings, stString_print("%2"PRId64"%c %2"PRId64,
                                                           r + 1, node->base,
                                                           (int64_t) (runLengths[r] / PAIR_ALIGNMENT_PROB_1)));
            }
        }
        char *label = stString_join2("\\n", labelStrings);
        fprintf(fH, "B%"PRId64" [label=\"%s\", fontcolor=%s, color=%s, penwidth=%f];\n", i, label,
                backboneColor, backboneColor, log(1 + weight));
        free(label);
        stList_destruct(labelStrings);

        if (i != 0) {
            fprintf(fH, "B%"PRId64" -> B%"PRId64" [label=\"%.2f\", fontcolor=%s, color=%s, weight=%d, penwidth=%f];\n",
                    i - 1, i, weight, backboneColor, backboneColor, (int) ceil(weight), log(1 + weight));
        }

        // Inserts
        for (int64_t j = 0; j < stList_length(node->inserts); j++) {
            PoaInsert *insert = stList_get(node->inserts, j);
            double iWeight = (insert->weightReverseStrand + insert->weightForwardStrand) / PAIR_ALIGNMENT_PROB_1;

            fprintf(fH, "I%"PRId64"_%"PRId64" [label=\"%s\", fontcolor=%s, color=%s, penwidth=%f];\n",
                    i, j, insert->insert->rleString, insertColor, insertColor, log(1 + iWeight));
            fprintf(fH,
                    "B%"PRId64" -> I%"PRId64"_%"PRId64" [label=\"%.2f\", fontcolor=%s, color=%s, weight=%d, penwidth=%f];\n",
                    i, i, j, iWeight, insertColor, insertColor, (int) ceil(iWeight), log(1 + iWeight));
            fprintf(fH, "I%"PRId64"_%"PRId64" -> B%"PRId64" [color=%s, weight=%d, penwidth=%f];\n",
                    i, j, i + 1, insertColor, (int) ceil(iWeight), log(1 + iWeight));
        }

        // Deletes
        for (int64_t j = 0; j < stList_length(node->deletes); j++) {
            PoaDelete *delete = stList_get(node->deletes, j);
            double dWeight = (delete->weightReverseStrand + delete->weightForwardStrand) / PAIR_ALIGNMENT_PROB_1;

            fprintf(fH, "B%"PRId64" -> B%"PRId64" [label=\"%.2f\", fontcolor=%s, color=%s, weight=%d, penwidth=%f];\n",
                    i, i + 1 + delete->length, dWeight, deleteColor, deleteColor, (int) ceil(dWeight),
                    log(1 + dWeight));
        }


    }

    fprintf(fH, "}\n");

}

void printMLRepeatCounts(RepeatSubMatrix *repeatSubMatrix, FILE *fh, Symbol base, stList *observations,
                         stList *bamChunkReads) {
    int64_t minRepeatLength, maxRepeatLength;

    // Calculate range of repeat counts observed
    repeatSubMatrix_getMinAndMaxRepeatCountObservations(repeatSubMatrix, observations,
                                                        bamChunkReads, &minRepeatLength, &maxRepeatLength);

    if (minRepeatLength == repeatSubMatrix->maximumRepeatLength) { // Case we have no valid observations
        for (int64_t i = 1; i < repeatSubMatrix->maximumRepeatLength; i++) {
            fprintf(fh, ",0");
        }
        return;
    }

    assert(maxRepeatLength - minRepeatLength >= 0);
    double logProbabilities[maxRepeatLength - minRepeatLength + 1];

    // Get weights for each repeat count
    repeatSubMatrix_getRepeatCountProbs(repeatSubMatrix, base, observations,
                                        bamChunkReads, logProbabilities, minRepeatLength, maxRepeatLength);

    // Calculate the normalizing constant for the probabilities
    double totalProb = LOG_ZERO;
    for (int64_t i = minRepeatLength; i <= maxRepeatLength; i++) {
        totalProb = stMath_logAddExact(logProbabilities[i - minRepeatLength] * 2.302585093,
                                       totalProb); // Use constant to move from base 10 to base e
    }

    // Print the repeat counts
    for (int64_t i = 1; i < minRepeatLength; i++) {
        fprintf(fh, ",0");
    }
    for (int64_t i = minRepeatLength; i <= maxRepeatLength; i++) {
        fprintf(fh, ",%f", (float) exp(logProbabilities[i - minRepeatLength] * 2.302585093 - totalProb));
    }
    for (int64_t i = maxRepeatLength + 1; i < repeatSubMatrix->maximumRepeatLength; i++) {
        fprintf(fh, ",0");
    }
}

static double nFloat(double numerator, double denominator) {
    return denominator == 0.0 ? 0.0 : numerator / denominator;
}

void poa_printCSV(Poa *poa, FILE *fH,
                  stList *bamChunkReads, RepeatSubMatrix *repeatSubMatrix,
                  float indelSignificanceThreshold) {

    fprintf(fH, "REF_INDEX,REF_BASE,REPEAT_COUNT,TOTAL_WEIGHT,FRACTION_POS_STRAND");
    for (int64_t j = 0; j < poa->alphabet->alphabetSize; j++) {
        char c = poa->alphabet->convertSymbolToChar(j);
        fprintf(fH, ",FRACTION_BASE_%c_WEIGHT,FRACTION_BASE_%c_POS_STRAND", c, c);
    }
    for (int64_t j = 1; j < repeatSubMatrix->maximumRepeatLength; j++) {
        fprintf(fH, ",PROB_REPEAT_COUNT_%" PRIi64 "", j);
    }
    fprintf(fH,
            ",INSERTS"); //xN(INSERT_SEQ,TOTAL_WEIGHT,FRACTION_POS_STRAND)'

    fprintf(fH,
            ",DELETES\n"); //xN(DELETE_LENGTH,TOTAL_WEIGHT,FRACTION_POS_STRAND)

    // Print info for each base in reference in turn
    for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
        PoaNode *node = stList_get(poa->nodes, i);

        // Print base weights first

        // Calculate strand specific base weights
        double totalWeight, totalPositiveWeight, totalNegativeWeight;
        double *baseWeights = poaNode_getStrandSpecificBaseWeights(node, bamChunkReads,
                                                                   &totalWeight, &totalPositiveWeight,
                                                                   &totalNegativeWeight, poa->alphabet, NULL);

        fprintf(fH, "%" PRIi64 ",%c,%" PRIi64 ",%f,%f", i, node->base, node->repeatCount,
                nFloat(totalWeight, PAIR_ALIGNMENT_PROB_1),
                nFloat(totalPositiveWeight, totalPositiveWeight + totalNegativeWeight));

        for (int64_t j = 0; j < poa->alphabet->alphabetSize; j++) {
            double positiveStrandBaseWeight = baseWeights[j * 2 + 1];
            double negativeStrandBaseWeight = baseWeights[j * 2 + 0];
            double totalBaseWeight = positiveStrandBaseWeight + negativeStrandBaseWeight;

            fprintf(fH, ",%f,%f", nFloat(node->baseWeights[j], totalWeight),
                    nFloat(positiveStrandBaseWeight, totalBaseWeight));
        }

        free(baseWeights);

        // Print repeat counts
        printMLRepeatCounts(repeatSubMatrix, fH, poa->alphabet->convertCharToSymbol(node->base),
                            node->observations, bamChunkReads);

        // Inserts
        fprintf(fH, ",");
        for (int64_t j = 0; j < stList_length(node->inserts); j++) {
            PoaInsert *insert = stList_get(node->inserts, j);
            if (poaInsert_getWeight(insert) / PAIR_ALIGNMENT_PROB_1 >= indelSignificanceThreshold) {
                char *s = rleString_expand(insert->insert);
                fprintf(fH, "|%s|%f|%f",
                        s, nFloat(poaInsert_getWeight(insert), PAIR_ALIGNMENT_PROB_1),
                        nFloat(insert->weightForwardStrand, poaInsert_getWeight(insert)));
                free(s);
            }
        }

        // Deletes
        fprintf(fH, ",");
        for (int64_t j = 0; j < stList_length(node->deletes); j++) {
            PoaDelete *delete = stList_get(node->deletes, j);
            if (poaDelete_getWeight(delete) / PAIR_ALIGNMENT_PROB_1 >= indelSignificanceThreshold) {
                fprintf(fH, "|%" PRIi64 "|%f|%f",
                        delete->length, nFloat(poaDelete_getWeight(delete), PAIR_ALIGNMENT_PROB_1),
                        nFloat(delete->weightForwardStrand, poaDelete_getWeight(delete)));
            }
        }
        fprintf(fH, "\n");
    }
}

void poa_printPhasedCSV_indelPrint(stList *observations, FILE *fH,
                                   stList *bamChunkReads, stSet *readsInHap1, stSet *readsInHap2) {
    double totalPositiveWeightHap1 = 0.0, totalNegativeWeightHap1 = 0.0;
    double totalPositiveWeightHap2 = 0.0, totalNegativeWeightHap2 = 0.0;

    for (int64_t k = 0; k < stList_length(observations); k++) {
        PoaBaseObservation *obs = stList_get(observations, k);
        BamChunkRead *read = stList_get(bamChunkReads, obs->readNo);
        if (stSet_search(readsInHap1, read) != NULL) {
            if (read->forwardStrand) {
                totalPositiveWeightHap1 += obs->weight;
            } else {
                totalNegativeWeightHap1 += obs->weight;
            }
        } else if (stSet_search(readsInHap2, read) != NULL) {
            if (read->forwardStrand) {
                totalPositiveWeightHap2 += obs->weight;
            } else {
                totalNegativeWeightHap2 += obs->weight;
            }
        }
    }

    double totalWeight =
            totalPositiveWeightHap1 + totalNegativeWeightHap1 + totalPositiveWeightHap2 + totalNegativeWeightHap2;

    fprintf(fH, "|%f|%f|%f|%f|%f",
            nFloat(totalWeight, PAIR_ALIGNMENT_PROB_1),
            nFloat(totalPositiveWeightHap1 + totalNegativeWeightHap1, totalWeight),
            nFloat(totalPositiveWeightHap2 + totalNegativeWeightHap2, totalWeight),
            nFloat(totalPositiveWeightHap1, totalPositiveWeightHap1 + totalNegativeWeightHap1),
            nFloat(totalPositiveWeightHap2, totalPositiveWeightHap2 + totalNegativeWeightHap2));
}

void poa_printPhasedCSV(Poa *poa, FILE *fH,
                        stList *bamChunkReads, stSet *readsInHap1, stSet *readsInHap2,
                        RepeatSubMatrix *repeatSubMatrix, float indelSignificanceThreshold) {

    fprintf(fH,
            "REF_INDEX,REF_BASE,REPEAT_COUNT,TOTAL_WEIGHT,FRACTION_HAP1_WEIGHT,FRACTION_HAP2_WEIGHT,FRACTION_POS_STRAND_HAP1,FRACTION_POS_STRAND_HAP2");
    for (int64_t j = 0; j < poa->alphabet->alphabetSize; j++) {
        char c = poa->alphabet->convertSymbolToChar(j);
        fprintf(fH,
                ",FRACTION_BASE_%c_WEIGHT,FRACTION_BASE_%c_HAP1,FRACTION_BASE_%c_HAP2,FRACTION_BASE_%c_POS_STRAND_HAP1,FRACTION_BASE_%c_POS_STRAND_HAP2",
                c, c, c, c, c);
    }
    for (int64_t j = 1; j < repeatSubMatrix->maximumRepeatLength; j++) {
        fprintf(fH, ",PROB_HAP1_REPEAT_COUNT_%" PRIi64 "", j);
    }
    for (int64_t j = 1; j < repeatSubMatrix->maximumRepeatLength; j++) {
        fprintf(fH, ",PROB_HAP2_REPEAT_COUNT_%" PRIi64 "", j);
    }
    fprintf(fH,
            ",INSERTS" //xN(INSERT_SEQ,TOTAL_WEIGHT,FRACTION_HAP1_WEIGHT,FRACTION_HAP2_WEIGHT,FRACTION_POS_STRAND_HAP1,FRACTION_POS_STRAND_HAP2)'
            ",DELETES\n"); //xN(DELETE_LENGTH,TOTAL_WEIGHT,FRACTION_HAP1_WEIGHT,FRACTION_HAP2_WEIGHT,FRACTION_POS_STRAND_HAP1,FRACTION_POS_STRAND_HAP2)'

    // Print info for each base in reference in turn
    for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
        PoaNode *node = stList_get(poa->nodes, i);

        // Print base weights

        // Calculate strand specific base weights
        double totalWeight, totalPositiveWeight, totalNegativeWeight;
        double *baseWeights = poaNode_getStrandSpecificBaseWeights(node, bamChunkReads,
                &totalWeight, &totalPositiveWeight, &totalNegativeWeight, poa->alphabet, NULL);

        // Calculate hap1, strand specific base weights
        double totalWeightHap1, totalPositiveWeightHap1, totalNegativeWeightHap1;
        double *baseWeightsHap1 = poaNode_getStrandSpecificBaseWeights(node, bamChunkReads,
                &totalWeightHap1, &totalPositiveWeightHap1, &totalNegativeWeightHap1, poa->alphabet, readsInHap1);

        // Calculate hap1, strand specific base weights
        double totalWeightHap2, totalPositiveWeightHap2, totalNegativeWeightHap2;
        double *baseWeightsHap2 = poaNode_getStrandSpecificBaseWeights(node, bamChunkReads,
                &totalWeightHap2, &totalPositiveWeightHap2, &totalNegativeWeightHap2, poa->alphabet, readsInHap2);

        //fprintf(fH, "REF_INDEX,REF_BASE,TOTAL_WEIGHT,FRACTION_HAP1_WEIGHT,FRACTION_HAP2_WEIGHT,FRACTION_POS_STRAND_HAP1,FRACTION_POS_STRAND_HAP2");
        fprintf(fH, "%" PRIi64 ",%c,%" PRIi64 ",%f,%f,%f,%f,%f", i, node->base, node->repeatCount,
                nFloat(totalWeight, PAIR_ALIGNMENT_PROB_1),
                nFloat(totalWeightHap1, totalWeight), nFloat(totalWeightHap2, totalWeight),
                nFloat(totalPositiveWeightHap1, totalWeightHap1), nFloat(totalPositiveWeightHap2, totalWeightHap2));

        for (int64_t j = 0; j < poa->alphabet->alphabetSize; j++) {
            double positiveStrandBaseWeight = baseWeights[j * 2 + 1];
            double negativeStrandBaseWeight = baseWeights[j * 2 + 0];
            double totalBaseWeight = positiveStrandBaseWeight + negativeStrandBaseWeight;

            double positiveStrandBaseWeightHap1 = baseWeightsHap1[j * 2 + 1];
            double negativeStrandBaseWeightHap1 = baseWeightsHap1[j * 2 + 0];

            double positiveStrandBaseWeightHap2 = baseWeightsHap2[j * 2 + 1];
            double negativeStrandBaseWeightHap2 = baseWeightsHap2[j * 2 + 0];

            //fprintf(fH, ",NORM_BASE_%c_WEIGHT,FRACTION_BASE_%c_HAP1,FRACTION_BASE_%c_HAP2,FRACTION_BASE_%c_POS_STRAND_HAP1,FRACTION_BASE_%c_POS_STRAND_HAP2", c, c, c, c);
            fprintf(fH, ",%f,%f,%f,%f,%f", nFloat(totalBaseWeight, totalWeight),
                    nFloat(positiveStrandBaseWeightHap1 + negativeStrandBaseWeightHap1, totalBaseWeight),
                    nFloat(positiveStrandBaseWeightHap2 + negativeStrandBaseWeightHap2, totalBaseWeight),
                    nFloat(positiveStrandBaseWeightHap1, positiveStrandBaseWeightHap1 + negativeStrandBaseWeightHap1),
                    nFloat(positiveStrandBaseWeightHap2, positiveStrandBaseWeightHap2 + negativeStrandBaseWeightHap2));
        }

        free(baseWeights);
        free(baseWeightsHap1);
        free(baseWeightsHap2);

        // Split observations between haplotypes (assume reads not in hap1 are in hap2)
        stList *observationsHap1 = stList_construct();
        stList *observationsHap2 = stList_construct();
        for (int64_t j = 0; j < stList_length(node->observations); j++) {
            PoaBaseObservation *observation = stList_get(node->observations, j);
            BamChunkRead *read = stList_get(bamChunkReads, observation->readNo);
            stList_append(stSet_search(readsInHap1, read) != NULL ? observationsHap1 : observationsHap2, observation);
        }

        // Print repeat counts for hap1
        printMLRepeatCounts(repeatSubMatrix, fH, poa->alphabet->convertCharToSymbol(node->base),
                            observationsHap1, bamChunkReads);

        // Print repeat counts for hap2
        printMLRepeatCounts(repeatSubMatrix, fH, poa->alphabet->convertCharToSymbol(node->base),
                            observationsHap2, bamChunkReads);

        // Cleanup
        stList_destruct(observationsHap1);
        stList_destruct(observationsHap2);

        // Inserts
        fprintf(fH, ",");
        for (int64_t j = 0; j < stList_length(node->inserts); j++) {
            PoaInsert *insert = stList_get(node->inserts, j);
            if (poaInsert_getWeight(insert) / PAIR_ALIGNMENT_PROB_1 >= indelSignificanceThreshold) {

                //fprintf(fH, "\tINSERTSxN(INSERT_SEQ, TOTAL_WEIGHT, FRACTION_HAP1_WEIGHT, FRACTION_HAP2_WEIGHT, FRACTION_POS_STRAND_HAP1, FRACTION_POS_STRAND_HAP2)\n");
                char *s = rleString_expand(insert->insert);
                fprintf(fH, "|%s", s);
                poa_printPhasedCSV_indelPrint(insert->observations, fH,
                                              bamChunkReads, readsInHap1, readsInHap2);
                free(s);
            }
        }

        // Deletes
        fprintf(fH, ",");
        for (int64_t j = 0; j < stList_length(node->deletes); j++) {
            PoaDelete *delete = stList_get(node->deletes, j);
            if (poaDelete_getWeight(delete) / PAIR_ALIGNMENT_PROB_1 >= indelSignificanceThreshold) {

                //fprintf(fH, "\tDELETESxN(DELETE_LENGTH, TOTAL_WEIGHT, FRACTION_HAP1_WEIGHT, FRACTION_HAP2_WEIGHT, FRACTION_POS_STRAND_HAP1, FRACTION_POS_STRAND_HAP2)\n");
                fprintf(fH, "|%" PRIi64 "", delete->length);
                poa_printPhasedCSV_indelPrint(delete->observations, fH,
                                              bamChunkReads, readsInHap1, readsInHap2);
            }
        }
        fprintf(fH, "\n");
    }
}

void poa_print(Poa *poa, FILE *fH,
               stList *bamChunkReads,
               float indelSignificanceThreshold) {
    // Print info for each base in reference in turn
    for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
        PoaNode *node = stList_get(poa->nodes, i);

        // Calculate strand specific base weights
        double totalWeight, totalPositiveWeight, totalNegativeWeight;
        double *baseWeights = poaNode_getStrandSpecificBaseWeights(node, bamChunkReads,
                                                                   &totalWeight, &totalPositiveWeight,
                                                                   &totalNegativeWeight, poa->alphabet, NULL);

        fprintf(fH, "%" PRIi64 "\t%c total-weight:%f\ttotal-pos-weight:%f\ttotal-neg-weight:%f", i, node->base,
                totalWeight / PAIR_ALIGNMENT_PROB_1, totalPositiveWeight / PAIR_ALIGNMENT_PROB_1,
                totalNegativeWeight / PAIR_ALIGNMENT_PROB_1);

        for (int64_t j = 0; j < poa->alphabet->alphabetSize; j++) {
            double positiveStrandBaseWeight = baseWeights[j * 2 + POS_STRAND_IDX];
            double negativeStrandBaseWeight = baseWeights[j * 2 + NEG_STRAND_IDX];
            double totalBaseWeight = positiveStrandBaseWeight + negativeStrandBaseWeight;

            if (totalBaseWeight / totalWeight > 0.25) {
                fprintf(fH, "\t%c:%f (%f) +str:%f, -str:%f,", poa->alphabet->convertSymbolToChar(j),
                        (float) node->baseWeights[j] / PAIR_ALIGNMENT_PROB_1, node->baseWeights[j] / totalWeight,
                        positiveStrandBaseWeight / totalPositiveWeight, negativeStrandBaseWeight / totalNegativeWeight);
            }
        }
        fprintf(fH, "\tTotal-weight:%f\n", (float) totalWeight / PAIR_ALIGNMENT_PROB_1);

        free(baseWeights);

        // Inserts
        for (int64_t j = 0; j < stList_length(node->inserts); j++) {
            PoaInsert *insert = stList_get(node->inserts, j);
            if (poaInsert_getWeight(insert) / PAIR_ALIGNMENT_PROB_1 >= indelSignificanceThreshold) {
                char *s = rleString_expand(insert->insert);
                fprintf(fH, "Insert\tSeq:%s\tTotal weight:%f\tForward Strand Weight:%f\tReverse Strand Weight:%f\n", s,
                        (float) poaInsert_getWeight(insert) / PAIR_ALIGNMENT_PROB_1,
                        (float) insert->weightForwardStrand / PAIR_ALIGNMENT_PROB_1,
                        (float) insert->weightReverseStrand / PAIR_ALIGNMENT_PROB_1);
                free(s);
            }
        }

        // Deletes
        for (int64_t j = 0; j < stList_length(node->deletes); j++) {
            PoaDelete *delete = stList_get(node->deletes, j);
            if (poaDelete_getWeight(delete) / PAIR_ALIGNMENT_PROB_1 >= indelSignificanceThreshold) {
                fprintf(fH,
                        "Delete\tLength:%" PRIi64 "\tTotal weight:%f\tForward Strand Weight:%f\tReverse Strand Weight:%f\n",
                        delete->length,
                        (float) poaDelete_getWeight(delete) / PAIR_ALIGNMENT_PROB_1,
                        (float) delete->weightForwardStrand / PAIR_ALIGNMENT_PROB_1,
                        (float) delete->weightReverseStrand / PAIR_ALIGNMENT_PROB_1);
            }
        }
    }
}

void poa_printSummaryStats(Poa *poa, FILE *fH) {
    double totalReferenceMatchWeight = poa_getReferenceNodeTotalMatchWeight(poa) / PAIR_ALIGNMENT_PROB_1;
    double totalReferenceMismatchWeight = poa_getReferenceNodeTotalDisagreementWeight(poa) / PAIR_ALIGNMENT_PROB_1;
    double totalInsertWeight = poa_getInsertTotalWeight(poa) / PAIR_ALIGNMENT_PROB_1;
    double totalDeleteWeight = poa_getDeleteTotalWeight(poa) / PAIR_ALIGNMENT_PROB_1;

    fprintf(fH,
            "Totals, reference match weight: %f reference mismatch weight: %f insert weight: %f delete weight: %f indel weight: %f, sum error: %f\n",
            totalReferenceMatchWeight, totalReferenceMismatchWeight,
            totalInsertWeight, totalDeleteWeight, totalInsertWeight + totalDeleteWeight,
            totalInsertWeight + totalDeleteWeight + totalReferenceMismatchWeight);
}

uint64_t
getMaxWeight(const double *weights, uint64_t weightNo, uint64_t referenceIndex, double referenceWeightPenalty) {
    // Used to pick bases and repeat counts by consensus algorithm
    double maxWeight = 0;
    int64_t maxIndex = -1;
    for (int64_t j = 0; j < weightNo; j++) {
        if (j != referenceIndex && weights[j] >= maxWeight) {
            maxWeight = weights[j];
            maxIndex = j;
        }
    }
    assert(maxIndex != -1);

    return weights[referenceIndex] * referenceWeightPenalty >= maxWeight ? referenceIndex : maxIndex;
}

RleString *poa_getConsensus(Poa *poa, int64_t **poaToConsensusMap, PolishParams *pp) {
    // Cheesy profile HMM like algorithm
    // Calculates forward probabilities through model, then
    // traces back through max prob local path, greedily

    // Probabilities/weights we keep track of.

    // Total weight of outgoing transitions
    double *totalOutgoingWeights = st_calloc(stList_length(poa->nodes), sizeof(double));

    // Forward probabilities
    double *nodeForwardLogProbs = st_calloc(stList_length(poa->nodes) + 1, sizeof(double));
    // Initialize, only start state has log(1) = 0 prob
    for (int64_t i = 1; i < stList_length(poa->nodes) + 1; i++) {
        nodeForwardLogProbs[i] = LOG_ZERO;
    }

    // Forward probabilities of transitioning from a node to the its successor without
    // an indel
    double *matchTransitionForwardLogProbs = st_calloc(stList_length(poa->nodes), sizeof(double));

    // Calculate incoming deletions for each node

    stList *incomingDeletions = stList_construct3(0, (void (*)(void *)) stList_destruct);
    for (int64_t i = 0; i < stList_length(poa->nodes) + 1; i++) {
        stList_append(incomingDeletions, stList_construct());
    }
    for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
        PoaNode *node = stList_get(poa->nodes, i);
        for (int64_t j = 0; j < stList_length(node->deletes); j++) {
            PoaDelete *delete = stList_get(node->deletes, j);
            stList_append(stList_get(incomingDeletions, i + delete->length + 1), delete);
        }
    }

    // Walk through the graph left-to-right calculating forward probabilities

    for (int64_t i = 0; i < stList_length(poa->nodes); i++) {
        PoaNode *node = stList_get(poa->nodes, i);

        // Calculate total weight of indels connecting from this node

        double totalIndelWeight = 0.0;

        for (int64_t j = 0; j < stList_length(node->inserts); j++) {
            PoaInsert *insert = stList_get(node->inserts, j);
            totalIndelWeight += poaInsert_getWeight(insert);
        }
        for (int64_t j = 0; j < stList_length(node->deletes); j++) {
            PoaDelete *delete = stList_get(node->deletes, j);
            totalIndelWeight += poaDelete_getWeight(delete);
        }

        // Calculate the match transition weight of node
        // that is, we make an estimate of the weight/expectation
        // of transitioning from this node to the next node without an indel

        double matchTransitionWeight = 0.0;
        if (i == 0) {
            if (stList_length(poa->nodes) == 1) { // In case is zero length reference
                matchTransitionWeight = 1.0;
            } else {
                // Set the initiation probability according to the average base weight
                for (int64_t j = 1; j < stList_length(poa->nodes); j++) {
                    PoaNode *nNode = stList_get(poa->nodes, j);
                    for (int64_t k = 0; k < poa->alphabet->alphabetSize; k++) {
                        matchTransitionWeight += nNode->baseWeights[k];
                    }
                }
                matchTransitionWeight /= (double) stList_length(poa->nodes) - 1;
                matchTransitionWeight -= totalIndelWeight;
            }
        } else {
            for (int64_t j = 0; j < poa->alphabet->alphabetSize; j++) {
                matchTransitionWeight += node->baseWeights[j];
            }
            matchTransitionWeight -= totalIndelWeight;
        }

        // Hack to stop zero weights
        matchTransitionWeight = matchTransitionWeight <= 0.0 ? 0.0001 : matchTransitionWeight; // Make a small value

        // Calculate the total weight of outgoing transitions
        totalOutgoingWeights[i] = matchTransitionWeight + totalIndelWeight;

        // Update the probabilities of nodes that connect by to this node

        // Inserts
        for (int64_t j = 0; j < stList_length(node->inserts); j++) {
            PoaInsert *insert = stList_get(node->inserts, j);
            nodeForwardLogProbs[i + 1] = logAdd(nodeForwardLogProbs[i + 1],
                                                nodeForwardLogProbs[i] +
                                                log(poaInsert_getWeight(insert) / totalOutgoingWeights[i]));
        }

        // Deletes
        for (int64_t j = 0; j < stList_length(node->deletes); j++) {
            PoaDelete *delete = stList_get(node->deletes, j);
            nodeForwardLogProbs[i + delete->length + 1] = logAdd(nodeForwardLogProbs[i + delete->length + 1],
                                                                 nodeForwardLogProbs[i] +
                                                                 log(poaDelete_getWeight(delete) /
                                                                     totalOutgoingWeights[i]));
        }

        // Match
        matchTransitionForwardLogProbs[i] =
                nodeForwardLogProbs[i] + log(matchTransitionWeight / totalOutgoingWeights[i]);
        nodeForwardLogProbs[i + 1] = logAdd(nodeForwardLogProbs[i + 1], matchTransitionForwardLogProbs[i]);
    }

    // Now traceback picking consensus greedily

    // Allocate consensus map, setting the alignment of reference
    // string positions initially all to gaps.
    *poaToConsensusMap = st_malloc((stList_length(poa->nodes) - 1) * sizeof(int64_t));
    for (int64_t i = 0; i < stList_length(poa->nodes) - 1; i++) {
        (*poaToConsensusMap)[i] = -1;
    }

    stList *consensusStrings = stList_construct3(0, free);
    int64_t runningConsensusLength = 0;
    char previousBase = '-';

    for (int64_t i = stList_length(poa->nodes); i > 0;) {

        //  Add base if not at end
        if (i < stList_length(poa->nodes)) {
            PoaNode *node = stList_get(poa->nodes, i);

            // Picks a base, giving a discount to the reference base,
            // because the alignment is biased towards it
            int64_t maxBaseIndex = getMaxWeight(node->baseWeights, poa->alphabet->alphabetSize,
                                                poa->alphabet->convertCharToSymbol(node->base),
                                                pp->referenceBasePenalty);
            char base = poa->alphabet->convertSymbolToChar(maxBaseIndex);

            if (pp->useRunLengthEncoding) {

                // Similarly picks a repeat count
                int64_t maxWeightRepeatCount = getMaxWeight(node->repeatCountWeights, poa->maxRepeatCount,
                                                            node->repeatCount, pp->referenceBasePenalty);
                maxWeightRepeatCount = maxWeightRepeatCount == 0 ? 1
                                                                 : maxWeightRepeatCount; // Avoid making a repeat count of zero (could
                // happen if coverage was 0 through region).

                // Add base * repeat count to list of consensus strings
                stList_append(consensusStrings, expandChar(base, maxWeightRepeatCount));

                // Update poa to consensus map and increase RLE consensus length
                if (previousBase != base) {
                    (*poaToConsensusMap)[i - 1] = runningConsensusLength++;
                }
                previousBase = base;
            } else { // Otherwise repeat counts always one
                stList_append(consensusStrings, expandChar(base, 1));
                (*poaToConsensusMap)[i - 1] = runningConsensusLength++;
            }
        }

        // Get max insert
        double maxInsertProb = LOG_ZERO;
        double totalInsertProb = LOG_ZERO;
        PoaInsert *maxInsert = NULL;
        PoaNode *pNode = stList_get(poa->nodes, i - 1);
        for (int64_t j = 0; j < stList_length(pNode->inserts); j++) {
            PoaInsert *insert = stList_get(pNode->inserts, j);
            double p = log(poaInsert_getWeight(insert) / totalOutgoingWeights[i - 1]) + nodeForwardLogProbs[i - 1];
            if (p > maxInsertProb) {
                maxInsertProb = p;
                maxInsert = insert;
            }
            totalInsertProb = logAdd(totalInsertProb, p);
        }

        // Get max delete
        double maxDeleteProb = LOG_ZERO;
        double totalDeleteProb = LOG_ZERO;
        PoaDelete *maxDelete = NULL;
        stList *incidentDeletes = stList_get(incomingDeletions, i);
        for (int64_t j = 0; j < stList_length(incidentDeletes); j++) {
            PoaDelete *delete = stList_get(incidentDeletes, j);
            double p = log(poaDelete_getWeight(delete) / totalOutgoingWeights[i - delete->length - 1]) +
                       nodeForwardLogProbs[i - delete->length - 1];
            if (p > maxDeleteProb) {
                maxDeleteProb = p;
                maxDelete = delete;
            }
            totalDeleteProb = logAdd(totalDeleteProb, p);
        }

        //fprintf(stderr, "%i Maxinsert: %f Maxdelete %f Match %f\n", (int)i, (float)maxInsertProb, (float)maxDeleteProb, (float)matchTransitionForwardLogProbs[i-1]);

        if (matchTransitionForwardLogProbs[i - 1] >= totalDeleteProb &&
            matchTransitionForwardLogProbs[i - 1] >= totalInsertProb) {
            // Is likely a match, move back to previous reference base
            i--;
        } else if (totalInsertProb >= totalDeleteProb) {
            // Is likely an insert, append insert to consensus string
            // and move to a previous reference base
            stList_append(consensusStrings, rleString_expand(maxInsert->insert));
            if (pp->useRunLengthEncoding) {
                char base = maxInsert->insert->rleString[maxInsert->insert->length - 1];
                runningConsensusLength += (int64_t) maxInsert->insert->length + (base != previousBase ? 0 : -1);
                previousBase = maxInsert->insert->rleString[0];
            } else {
                assert(maxInsert->insert->nonRleLength == maxInsert->insert->length);
                runningConsensusLength += maxInsert->insert->nonRleLength;
            }
            i--;
        } else {
            // Is likely a delete, jump back to skip deleted bases
            i -= maxDelete->length + 1;
        }
    }

    // Concatenate backwards to make consensus string
    stList_reverse(consensusStrings);
    char *expandedConsensusString = stString_join2("", consensusStrings);
    RleString *consensusString = pp->useRunLengthEncoding ? rleString_construct(expandedConsensusString)
                                                          : rleString_construct_no_rle(expandedConsensusString);
    free(expandedConsensusString);
    assert(runningConsensusLength == consensusString->length);

    // Now reverse the poaToConsensusMap, because offsets are from  end of string but need them to be from beginning
    for (int64_t i = 0; i < stList_length(poa->nodes) - 1; i++) {
        if ((*poaToConsensusMap)[i] != -1) {
            (*poaToConsensusMap)[i] = consensusString->length - 1 - (*poaToConsensusMap)[i];
        }
    }

    // Cleanup
    stList_destruct(consensusStrings);
    stList_destruct(incomingDeletions);
    free(nodeForwardLogProbs);
    free(matchTransitionForwardLogProbs);
    free(totalOutgoingWeights);

    return consensusString;
}

char *addInsert(char *string, char *insertString, int64_t editStart) {
    int64_t insertLength = strlen(insertString);
    int64_t stringLength = strlen(string);

    // Allocate new string
    char *editedString = st_malloc(sizeof(char) * (stringLength + insertLength + 1));

    // Make new reference string
    memcpy(editedString, string, sizeof(char) * editStart);
    memcpy(&(editedString[editStart]), insertString, sizeof(char) * insertLength);
    memcpy(&(editedString[editStart + insertLength]), &(string[editStart]),
           sizeof(char) * (1 + stringLength - editStart));

    return editedString;
}

char *removeDelete(char *string, int64_t deleteLength, int64_t editStart) {
    int64_t stringLength = strlen(string);

    // Allocate new string
    assert(stringLength >= deleteLength);
    char *editedString = st_malloc(sizeof(char) * (stringLength - deleteLength + 1));

    // Make new reference string
    memcpy(editedString, string, sizeof(char) * editStart);
    memcpy(&(editedString[editStart]), &(string[editStart + deleteLength]),
           sizeof(char) * (1 + stringLength - editStart - deleteLength));

    return editedString;
}

stList *poa_getReadAlignmentsToConsensus(Poa *poa, stList *bamChunkReads, PolishParams *polishParams) {
    // Generate anchor alignments
    stList *anchorAlignments = poa_getAnchorAlignments(poa, NULL, stList_length(bamChunkReads), polishParams);

    // Alignments
    stList *alignments = stList_construct3(0, (void (*)(void *)) stList_destruct);

    // Make the MEA alignments
    SymbolString refSymbolString = rleString_constructSymbolString(poa->refString, 0, poa->refString->length,
                                                                   polishParams->alphabet,
                                                                   polishParams->useRepeatCountsInAlignment,
                                                                   poa->maxRepeatCount);
    for (int64_t i = 0; i < stList_length(bamChunkReads); i++) {
        BamChunkRead *read = stList_get(bamChunkReads, i);

        // Generate the posterior alignment probabilities
        stList *matches, *inserts, *deletes;
        getAlignedPairsWithIndelsCroppingReference(poa->refString, read->rleRead, read->forwardStrand,
                                                   stList_get(anchorAlignments, i), &matches, &inserts, &deletes,
                                                   polishParams);

        // Get the MEA alignment
        double alignmentScore;
        stList *alignment = getMaximalExpectedAccuracyPairwiseAlignment(matches, deletes, inserts,
                                                                        poa->refString->length, read->rleRead->length,
                                                                        &alignmentScore, polishParams->p);

        // Symbol strings
        SymbolString readSymbolString = rleString_constructSymbolString(read->rleRead, 0, read->rleRead->length,
                                                                        polishParams->alphabet,
                                                                        polishParams->useRepeatCountsInAlignment,
                                                                        poa->maxRepeatCount);

        // Left shift the alignment
        stList *leftShiftedAlignment = leftShiftAlignment(alignment, refSymbolString, readSymbolString);

        // Cleanup
        stList_destruct(inserts);
        stList_destruct(deletes);
        stList_destruct(matches);
        stList_destruct(alignment);
        symbolString_destruct(readSymbolString);

        stList_append(alignments, leftShiftedAlignment);
    }

    // Cleanup
    stList_destruct(anchorAlignments);
    symbolString_destruct(refSymbolString);

    return alignments;
}

/*
 * Functions for run-length decoding with POAs
 */

int64_t getRunLengthMode(Alphabet *alphabet, Symbol base, stList *observations, stList *bamChunkReads) {
    stHash *runLengths = stHash_construct();
    int64_t maxCount = 0;
    int64_t maxRL = 0;
    for (int64_t i = 0; i < stList_length(observations); i++) {
        PoaBaseObservation *obs = stList_get(observations, i);
        BamChunkRead *bcr = stList_get(bamChunkReads, obs->readNo);
        RleString *rleString = bcr->rleRead;
        if (alphabet->convertCharToSymbol(rleString->rleString[obs->offset]) != base) continue;
        int64_t obvsRL = rleString->repeatCounts[obs->offset];
        int64_t currRlCount = (int64_t) stHash_remove(runLengths, (void *) obvsRL) + 1;
        if (currRlCount > maxCount) {
            maxCount = currRlCount;
            maxRL = obvsRL;
        }
        stHash_insert(runLengths, (void *) obvsRL, (void *) currRlCount);
    }
    stHash_destruct(runLengths);

    return maxRL;
}

static int64_t expandRLEConsensus2(Poa *poa, PoaNode *node, stList *bamChunkReads, RepeatSubMatrix *repeatSubMatrix) {
//	assert(repeatSubMatrix != NULL);
    // for an experiment, or case with no repeat matrix
    if (repeatSubMatrix == NULL) {
        return getRunLengthMode(poa->alphabet, poa->alphabet->convertCharToSymbol(node->base), node->observations,
                                bamChunkReads);
    }

    // Repeat count
    double logProbability;
    return repeatSubMatrix_getMLRepeatCount(repeatSubMatrix, poa->alphabet->convertCharToSymbol(node->base),
                                            node->observations,
                                            bamChunkReads, &logProbability);
}

void poa_estimateRepeatCountsUsingBayesianModel(Poa *poa, stList *bamChunkReads, RepeatSubMatrix *repeatSubMatrix) {
    poa->refString->nonRleLength = 0;
    for (uint64_t i = 1; i < stList_length(poa->nodes); i++) {
        PoaNode *node = stList_get(poa->nodes, i);
        poa->refString->repeatCounts[i - 1] = expandRLEConsensus2(poa, node, bamChunkReads,
                                                                  repeatSubMatrix);
        if (poa->refString->repeatCounts[i - 1] == 0) { // Prevent zero length estimates
            poa->refString->repeatCounts[i - 1] = 1;
        }
        node->repeatCount = poa->refString->repeatCounts[i - 1];
        poa->refString->nonRleLength += poa->refString->repeatCounts[i - 1];
    }
}

void poa_estimatePhasedRepeatCountsUsingBayesianModel(Poa *poa, stList *bamChunkReads,
                                                      RepeatSubMatrix *repeatSubMatrix, stSet *readsBelongingToHap1,
                                                      stSet *readsBelongingToHap2, PolishParams *params) {
    poa->refString->nonRleLength = 0;
    for (uint64_t i = 1; i < stList_length(poa->nodes); i++) {
        PoaNode *node = stList_get(poa->nodes, i);

        // Repeat count
        double logProbability;
        poa->refString->repeatCounts[i - 1] = repeatSubMatrix_getPhasedMLRepeatCount(repeatSubMatrix,
                                                                                     poa->refString->repeatCounts[i -
                                                                                                                  1],
                                                                                     poa->alphabet->convertCharToSymbol(
                                                                                             node->base),
                                                                                     node->observations,
                                                                                     bamChunkReads, &logProbability,
                                                                                     readsBelongingToHap1,
                                                                                     readsBelongingToHap2, params);

        if (poa->refString->repeatCounts[i - 1] == 0) { // Prevent zero length estimates
            poa->refString->repeatCounts[i - 1] = 1;
        }

        node->repeatCount = poa->refString->repeatCounts[i - 1]; // Update the repeat count of the node

        poa->refString->nonRleLength += poa->refString->repeatCounts[i - 1]; // Update the length of non-rle refString
    }
}

/*
 * Code to restimate bases using the phasing
 */

void getPhasedBaseWeights(double *logProbabilitiesHap, stList *observations, Poa *poa, stList *bamChunkReads) {
    for (int64_t i = 0; i < poa->alphabet->alphabetSize; i++) { // Initialize memory
        logProbabilitiesHap[i] = 0.0;
    }
    for (int64_t i = 0; i < stList_length(observations); i++) {
        PoaBaseObservation *observation = stList_get(observations, i);
        BamChunkRead *read = stList_get(bamChunkReads, observation->readNo);
        int64_t base = poa->alphabet->convertCharToSymbol(read->rleRead->rleString[observation->offset]);
        assert(base >= 0);
        assert(base < poa->alphabet->alphabetSize);
        logProbabilitiesHap[base] += observation->weight;
    }
    for (int64_t i = 0; i < poa->alphabet->alphabetSize; i++) { // Normalize
        logProbabilitiesHap[i] /= PAIR_ALIGNMENT_PROB_1;
    }
}

double getBaseProbForHaplotype(int64_t base, double *logProbabilitiesHap1, const double *logProbabilitiesHap2,
                               double logProbMLHap2, double logProbSubstitution) {
    double logProbHap2Same = logProbabilitiesHap2[base];
    return logProbabilitiesHap1[base] +
           ((logProbHap2Same > logProbMLHap2 + logProbSubstitution) ? logProbHap2Same : logProbMLHap2 +
                                                                                        logProbSubstitution);
}

uint64_t poaNode_getPhasedMLBase(PoaNode *node, stList *bamChunkReads, Poa *poa, stSet *readsBelongingToHap1) {

    // Split observations between haplotypes
    stList *observationsHap1 = stList_construct();
    stList *observationsHap2 = stList_construct();
    for (int64_t i = 0; i < stList_length(node->observations); i++) {
        PoaBaseObservation *observation = stList_get(node->observations, i);
        BamChunkRead *read = stList_get(bamChunkReads, observation->readNo);
        stList_append(stSet_search(readsBelongingToHap1, read) != NULL ? observationsHap1 : observationsHap2,
                      observation);
    }

    // Get probs for hap 1
    double logProbabilitiesHap1[poa->alphabet->alphabetSize];
    getPhasedBaseWeights(logProbabilitiesHap1, observationsHap1, poa, bamChunkReads);

    // Get probs for hap 2
    double logProbabilitiesHap2[poa->alphabet->alphabetSize];
    getPhasedBaseWeights(logProbabilitiesHap2, observationsHap2, poa, bamChunkReads);

    // Get ML prob for haplotype 2
    double logProbMLHap2;
    int64_t mLBaseHap2 = getMax(logProbabilitiesHap2, poa->alphabet->alphabetSize, &logProbMLHap2);

    // Calculate ML probability of repeat length for haplotype1
    double mlLogProb = getBaseProbForHaplotype(0, logProbabilitiesHap1, logProbabilitiesHap2,
                                               logProbMLHap2, log(0.1));
    int64_t mlBase = 0;
    for (int64_t i = 1; i < poa->alphabet->alphabetSize; i++) {
        double p = getBaseProbForHaplotype(i, logProbabilitiesHap1, logProbabilitiesHap2,
                                           logProbMLHap2, log(0.1));
        if (p > mlLogProb) {
            mlLogProb = p;
            mlBase = i;
        }
    }

    if (mlBase != poa->alphabet->convertCharToSymbol(node->base)) {
        st_logDebug(
                "Got %c base, previous base: %c, other hap ML base: %c (observation # hap1: %i, observqtion # hap2: %i\n",
                poa->alphabet->convertSymbolToChar(mlBase), node->base, poa->alphabet->convertSymbolToChar(mLBaseHap2),
                (int) stList_length(observationsHap1), (int) stList_length(observationsHap2));
    }

    // Cleanup
    stList_destruct(observationsHap1);
    stList_destruct(observationsHap2);

    return mlBase;
}

void poa_estimatePhasedBasesUsingBayesianModel(Poa *poa, stList *bamChunkReads, stSet *readsBelongingToHap1,
                                               PolishParams *params) {
    poa->refString->nonRleLength = 0;
    for (uint64_t i = 1; i < stList_length(poa->nodes); i++) {
        PoaNode *node = stList_get(poa->nodes, i);
        node->base = poa->alphabet->convertSymbolToChar(poaNode_getPhasedMLBase(node, bamChunkReads,
                                                                                poa, readsBelongingToHap1));
        poa->refString->rleString[i - 1] = node->base;
    }
}


// Core polishing logic functions

RleString *poa_polish(Poa *poa, stList *bamChunkReads, PolishParams *params,
                      int64_t **poaToConsensusMap) {
    /*
     * "Polishes" the given POA reference string to create a new consensus reference string.
     * Algorithm starts by dividing the reference into anchor points - points where the majority
     * of bamChunkReads are confidently aligned. Then, between each pair of consecutive anchors it looks
     * for "candidate variants", variants (either substitutions, insertions or deletions) with
     * significant weight. For every combination of candidate variants between the two anchor points
     * a candidate string is constructed and all the read substrings mapping between the two anchor points
     * are aligned to it. The candidate string, including the current reference substring,
     * with the highest likelihood is then selected.
     */
    BubbleGraph *bg = bubbleGraph_constructFromPoa(poa, bamChunkReads, params);

    uint64_t *consensusPath = bubbleGraph_getConsensusPath(bg, params);
    RleString *newConsensusString = bubbleGraph_getConsensusString(bg, consensusPath, poaToConsensusMap, params);

    free(consensusPath);
    bubbleGraph_destruct(bg);

    return newConsensusString;
}

// Functions to iteratively polish a sequence
Poa *poa_realignIterative(Poa *poa, stList *bamChunkReads, PolishParams *polishParams, bool hmmNotRealign,
                          int64_t minIterations, int64_t maxIterations) {
    assert(maxIterations >= 0);
    assert(minIterations <= maxIterations);

    time_t startTime = time(NULL);
    char *logIdentifier = getLogIdentifier();

    double score = poa_getReferenceNodeTotalMatchWeight(poa) - poa_getTotalErrorWeight(poa);

    st_logInfo(" %s Starting realignment with score: %6.4f\n", logIdentifier, score / PAIR_ALIGNMENT_PROB_1);

    int64_t i = 0;
    while (i < maxIterations) {
        i++;

        time_t consensusFindingStartTime = time(NULL);

        int64_t *poaToConsensusMap;
        RleString *reference = hmmNotRealign ? poa_getConsensus(poa, &poaToConsensusMap, polishParams) :
                               poa_polish(poa, bamChunkReads, polishParams, &poaToConsensusMap);

        st_logInfo(" %s Took %3d seconds to do round %" PRIi64 " of consensus finding using algorithm %s\n",
                   logIdentifier, (int) (time(NULL) - consensusFindingStartTime), i,
                   hmmNotRealign ? "consensus" : "polish");

        // Stop in case consensus string is same as old reference (i.e. greedy convergence)
        if (rleString_eq(reference, poa->refString)) {
            rleString_destruct(reference);
            free(poaToConsensusMap);
            break;
        }

        // Get anchor alignments
        stList *anchorAlignments = poa_getAnchorAlignments(poa, poaToConsensusMap, stList_length(bamChunkReads),
                                                           polishParams);

        time_t realignStartTime = time(NULL);

        // Generated updated poa
        Poa *poa2 = poa_realign(bamChunkReads, anchorAlignments, reference, polishParams);

        // Get updated repeat counts
        if (polishParams->useRunLengthEncoding) {
            poa_estimateRepeatCountsUsingBayesianModel(poa2, bamChunkReads, polishParams->repeatSubMatrix);
        }

        // Cleanup
        rleString_destruct(reference);
        free(poaToConsensusMap);
        stList_destruct(anchorAlignments);

        double score2 = poa_getReferenceNodeTotalMatchWeight(poa2) - poa_getTotalErrorWeight(poa2);

        st_logInfo(" %s Took %3d seconds to do round %" PRIi64 " of realignment, Have score: %6.4f (%f score diff)\n",
                   logIdentifier, (int) (time(NULL) - realignStartTime), i, score2 / PAIR_ALIGNMENT_PROB_1,
                   (score2 - score) / PAIR_ALIGNMENT_PROB_1);

        // Stop if score decreases (greedy stopping)
        if (score2 <= score && i > minIterations) {
            poa_destruct(poa2);
            break;
        }

        poa_destruct(poa);
        poa = poa2;
        score = score2;
    }

    st_logInfo(
            " %s Took %3d seconds to realign iterative using algorithm: %s through %" PRIi64 " iterations, got final score : %6.4f\n",
            logIdentifier, (int) (time(NULL) - startTime), hmmNotRealign ? "consensus" : "polish", i,
            score / PAIR_ALIGNMENT_PROB_1);

    free(logIdentifier);
    return poa;
}


Poa *poa_realignAll(stList *bamChunkReads, stList *anchorAlignments, RleString *reference,
                    PolishParams *polishParams) {
    time_t startTime = time(NULL);
    Poa *poa = poa_realign(bamChunkReads, anchorAlignments, reference, polishParams);
    char *logIdentifier = getLogIdentifier();

    st_logInfo(" %s Took %3d seconds to generate initial POA\n", logIdentifier, (int) (time(NULL) - startTime));
    free(logIdentifier);

    if (polishParams->maxPoaConsensusIterations > 0) {
        poa = poa_realignIterative(poa, bamChunkReads, polishParams, 1,
                polishParams->minPoaConsensusIterations, polishParams->maxPoaConsensusIterations);
    }

    if (polishParams->maxRealignmentPolishIterations > 0) {
        poa = poa_realignIterative(poa, bamChunkReads, polishParams, 0,
                polishParams->minRealignmentPolishIterations, polishParams->maxRealignmentPolishIterations);
    }

    return poa;
}
