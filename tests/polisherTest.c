/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "margin.h"

static char *polishParamsFile = "../params/ont/r9.4/allParams.np.human.r94-g344.json";
#define TEST_POLISH_FILES_DIR "../tests/data/polishTestExamples/"

Params *getParams() {
    Params *params = params_readParams(polishParamsFile);
    return params;
}

static void test_poa_getReferenceGraph(CuTest *testCase) {
    /*
     * Test building a trivial poa graph containing just a reference string.
     */

    RleString *reference = rleString_construct("GATTACA");

    Alphabet *alphabet = alphabet_constructNucleotide();
    Params *p = getParams();
    Poa *poa = poa_getReferenceGraph(reference, alphabet, p->polishParams->repeatSubMatrix->maximumRepeatLength);

    CuAssertTrue(testCase, stList_length(poa->nodes) == reference->length + 1);
    for (int64_t i = 0; i < reference->length; i++) {
        PoaNode *node = stList_get(poa->nodes, i + 1);

        CuAssertTrue(testCase, node->base == reference->rleString[i]);
        CuAssertTrue(testCase, stList_length(node->inserts) == 0);
        CuAssertTrue(testCase, stList_length(node->deletes) == 0);
    }

    PoaNode *node = stList_get(poa->nodes, 0);
    CuAssertTrue(testCase, node->base == 'N');
    CuAssertTrue(testCase, stList_length(node->inserts) == 0);
    CuAssertTrue(testCase, stList_length(node->deletes) == 0);

    poa_destruct(poa);
    rleString_destruct(reference);
    params_destruct(p);
}

static char *makeShiftedString(char *str, char *insert, int64_t insertPoint) {
    char *suffix = stString_copy(&str[insertPoint]);
    char *prefix = stString_copy(str);
    prefix[insertPoint] = '\0';
    char *shiftedStr = stString_print("%s%s%s", prefix, insert, suffix);
    free(suffix);
    free(prefix);
    return shiftedStr;
}

static void test_getShift(CuTest *testCase) {
    /*
     * Test left shifting code.
     */
    for (int64_t test = 0; test < 10000; test++) {
        // Make random rle string
        int64_t length = st_randomInt(1, 20);
        char *str = getRandomACGTSequence(length);
        RleString *str_rle = rleString_construct(str);

        // Make random rle insert of length m
        int64_t m = st_randomInt(1, 4);
        char *insert = getRandomACGTSequence(m);
        RleString *insert_rle = rleString_construct(insert);

        // Run get shift, getting shift in rle space
        int64_t i = getShift(str_rle, str_rle->length, insert_rle, 1);

        int64_t k = 0; // Calculate shift in non-rle space
        for (int64_t j = 0; j < i; j++) {
            k += str_rle->repeatCounts[j];
        }

        //if(k < length) {
        //	fprintf(stderr, "Str: %s, str-length:%" PRIi64 " insert: %s, insert:%" PRIi64 "\n", str, length, insert, k);
        //}

        // Test resulting transplanted string is same as concatenated str+insert
        char *shiftedStr = makeShiftedString(str, insert, k);
        char *concatenatedStr = stString_print("%s%s", str, insert);

        CuAssertStrEquals(testCase, concatenatedStr, shiftedStr);

        // Cleanup
        free(shiftedStr);

        // Test no further left shift would work in RLE space
        // TODO: REVISE THIS, as these checks don't work, consider
        // GGTGTTGT, str-length:8 insert: TGT, where the TGT can be inserted as position 2
        for (int64_t j = 0; j < k; j++) {
            if (j == 0 ||
                str[j - 1] != str[j]) { // This prevents shifts between what are the same positions in RLE space
                shiftedStr = makeShiftedString(str, insert, j);
                //CuAssertTrue(testCase, !stString_eq(shiftedStr, concatenatedStr)); // TODO FIX THIS
                free(shiftedStr);
            }
        }

        // Cleanup
        free(concatenatedStr);
        free(str);
        free(insert);
        rleString_destruct(str_rle);
        rleString_destruct(insert_rle);
    }
}

static void checkInserts(CuTest *testCase, Poa *poa, int64_t nodeIndex,
                         int64_t insertNumber, const char **inserts, const double *insertWeights, bool divideWeights) {
    PoaNode *node = stList_get(poa->nodes, nodeIndex);

    CuAssertIntEquals(testCase, stList_length(node->inserts), insertNumber);

    for (int64_t i = 0; i < insertNumber; i++) {
        PoaInsert *poaInsert = stList_get(node->inserts, i);
        CuAssertStrEquals(testCase, inserts[i], poaInsert->insert->rleString);
        CuAssertDblEquals(testCase, insertWeights[i],
                          poaInsert_getWeight(poaInsert) / (divideWeights ? PAIR_ALIGNMENT_PROB_1 : 1.0), 0.001);
    }
}

static void checkDeletes(CuTest *testCase, Poa *poa, int64_t nodeIndex,
                         int64_t deleteNumber, const int64_t *deleteLengths, const double *deleteWeights,
                         bool divideWeights) {
    PoaNode *node = stList_get(poa->nodes, nodeIndex);

    CuAssertIntEquals(testCase, stList_length(node->deletes), deleteNumber);

    for (int64_t i = 0; i < deleteNumber; i++) {
        PoaDelete *poaDelete = stList_get(node->deletes, i);
        CuAssertIntEquals(testCase, deleteLengths[i], poaDelete->length);
        CuAssertDblEquals(testCase, deleteWeights[i],
                          poaDelete_getWeight(poaDelete) / (divideWeights ? PAIR_ALIGNMENT_PROB_1 : 1.0), 0.001);
    }
}

static void checkNode(CuTest *testCase, Poa *poa, int64_t nodeIndex, char base, const double *baseWeights,
                      int64_t insertNumber, const char **inserts, const double *insertWeights,
                      int64_t deleteNumber, const int64_t *deleteLengths, const double *deleteWeights) {

    PoaNode *node = stList_get(poa->nodes, nodeIndex);
    CuAssertTrue(testCase, node->base == base);

    // Matches
    for (int64_t i = 0; i < poa->alphabet->alphabetSize; i++) {
        CuAssertDblEquals(testCase, node->baseWeights[i], baseWeights[i], 0.0);
    }

    // Inserts
    checkInserts(testCase, poa, nodeIndex, insertNumber, inserts, insertWeights, 0);

    // Deletes
    checkDeletes(testCase, poa, nodeIndex, deleteNumber, deleteLengths, deleteWeights, 0);
}

static void test_poa_augment_example(CuTest *testCase) {
    /*
     * Test poa_augment gives works as expected on a small example.
     */

    RleString *reference = rleString_construct_no_rle("GATTACA");

    Params *p = getParams();
    Alphabet *alphabet = alphabet_constructNucleotide();
    Poa *poa = poa_getReferenceGraph(reference, alphabet, p->polishParams->repeatSubMatrix->maximumRepeatLength);

    RleString *read = rleString_construct_no_rle("GATACGGT");

    stList *matches = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    stList *inserts = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    stList *deletes = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    stList_append(matches, stIntTuple_construct3(100, 0, 0));
    stList_append(matches, stIntTuple_construct3(100, 1, 1));
    stList_append(matches, stIntTuple_construct3(50, 2, 2));
    stList_append(matches, stIntTuple_construct3(50, 3, 2));
    stList_append(matches, stIntTuple_construct3(100, 4, 3));
    stList_append(matches, stIntTuple_construct3(100, 5, 4));
    stList_append(matches, stIntTuple_construct3(50, 6, 5));
    stList_append(matches, stIntTuple_construct3(25, 6, 6));
    stList_append(matches, stIntTuple_construct3(25, 6, 7));

    stList_append(inserts, stIntTuple_construct3(50, 5, 5));
    stList_append(inserts, stIntTuple_construct3(25, 5, 6));
    stList_append(inserts, stIntTuple_construct3(50, 6, 6));
    stList_append(inserts, stIntTuple_construct3(75, 6, 7));

    stList_append(deletes, stIntTuple_construct3(50, 2, 1));
    stList_append(deletes, stIntTuple_construct3(50, 3, 2));

    poa_augment(poa, read, 1, 0, matches, inserts, deletes, p->polishParams);

    // Check POA graph is what we expect

    CuAssertTrue(testCase, stList_length(poa->nodes) == 8); // Length + prefix node

    checkNode(testCase, poa, 0, 'N', (const double[]) {0.0, 0.0, 0.0, 0.0, 0.0},
              0, (const char *[]) {""}, (const double[]) {0.0},
              0, (const int64_t[]) {0}, (const double[]) {0.0});

    checkNode(testCase, poa, 1, 'G', (const double[]) {0.0, 0.0, 100.0, 0.0, 0.0},
              0, (const char *[]) {""}, (const double[]) {0.0},
              0, (const int64_t[]) {0}, (const double[]) {0.0});

    checkNode(testCase, poa, 2, 'A', (const double[]) {100.0, 0.0, 0.0, 0.0, 0.0},
              0, (const char *[]) {""}, (const double[]) {0.0},
              1, (const int64_t[]) {1}, (const double[]) {100.0});

    checkNode(testCase, poa, 3, 'T', (const double[]) {0.0, 0.0, 0.0, 50.0, 0.0},
              0, (const char *[]) {""}, (const double[]) {0.0},
              0, (const int64_t[]) {0}, (const double[]) {0.0});

    checkNode(testCase, poa, 4, 'T', (const double[]) {0.0, 0.0, 0.0, 50.0, 0.0},
              0, (const char *[]) {""}, (const double[]) {0.0},
              0, (const int64_t[]) {0}, (const double[]) {0.0});

    checkNode(testCase, poa, 5, 'A', (const double[]) {100.0, 0.0, 0.0, 0.0, 0.0},
              0, (const char *[]) {""}, (const double[]) {0.0},
              0, (const int64_t[]) {0}, (const double[]) {0.0});

    checkNode(testCase, poa, 6, 'C', (const double[]) {0.0, 100.0, 0.0, 0.0, 0.0},
              2, (const char *[]) {"G", "GG"}, (const double[]) {50.0, 25.0},
              0, (const int64_t[]) {0}, (const double[]) {0.0});

    checkNode(testCase, poa, 7, 'A', (const double[]) {0.0, 0.0, 75.0, 25.0, 0.0},
              2, (const char *[]) {"GT", "T"}, (const double[]) {50.0, 75.0},
              0, (const int64_t[]) {0}, (const double[]) {0.0});

    // Cleanup
    poa_destruct(poa);
    stList_destruct(matches);
    stList_destruct(inserts);
    stList_destruct(deletes);
    rleString_destruct(read);
    rleString_destruct(reference);
    params_destruct(p);
}

static void test_poa_realign_tiny_example1(CuTest *testCase) {
    /*
     * Tests that poa_realign builds the expected poa graph for a small example of input sequences
     */

    RleString *reference = rleString_construct_no_rle("GATACAGCGGG");
    BamChunkRead *read = bamChunkRead_construct2(stString_print("read"), stString_print("GATTACAGCG"), NULL, 1, 0);

    stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
    stList_append(reads, read);
    bool readStrand = 1;

    Params *params = params_readParams(polishParamsFile);
    PolishParams *polishParams = params->polishParams;

    Poa *poa = poa_realign(reads, NULL, reference, polishParams);
    // Check we get the set of inserts and deletes we expect

    // Basically two alignments probable:
    // One:
    //GAT-ACAGCGGG
    //GATTACAGC--G

    // Two:
    //GAT-ACAGCGGG
    //GATTACAGC-G-

    //Read:		GATTACAGCG
    //Reference:	GATACAGCGGG
    //0	N total-weight:0.000000	total-pos-weight:0.000000	total-neg-weight:0.000000	Total-weight:0.000000
    //1	G total-weight:0.999978	total-pos-weight:0.999978	total-neg-weight:0.000000	G:0.999978 (1.000000) +str:1.000000, -str:nan,	Total-weight:0.999978
    //2	A total-weight:0.994884	total-pos-weight:0.994884	total-neg-weight:0.000000	A:0.994884 (1.000000) +str:1.000000, -str:nan,	Total-weight:0.994884
    //Insert	Seq:T	Total weight:0.990160	Forward Strand Weight:0.990160	Reverse Strand Weight:0.000000
    //3	T total-weight:0.999554	total-pos-weight:0.999554	total-neg-weight:0.000000	T:0.999554 (1.000000) +str:1.000000, -str:nan,	Total-weight:0.999554
    //4	A total-weight:0.994393	total-pos-weight:0.994393	total-neg-weight:0.000000	A:0.994393 (1.000000) +str:1.000000, -str:nan,	Total-weight:0.994393
    //5	C total-weight:0.999478	total-pos-weight:0.999478	total-neg-weight:0.000000	C:0.999478 (1.000000) +str:1.000000, -str:nan,	Total-weight:0.999478
    //6	A total-weight:0.999600	total-pos-weight:0.999600	total-neg-weight:0.000000	A:0.999600 (1.000000) +str:1.000000, -str:nan,	Total-weight:0.999600
    //7	G total-weight:0.990489	total-pos-weight:0.990489	total-neg-weight:0.000000	G:0.990489 (1.000000) +str:1.000000, -str:nan,	Total-weight:0.990489
    //8	C total-weight:0.986294	total-pos-weight:0.986294	total-neg-weight:0.000000	C:0.986294 (1.000000) +str:1.000000, -str:nan,	Total-weight:0.986294
    //Delete	Length:1	Total weight:1.045831	Forward Strand Weight:1.045831	Reverse Strand Weight:0.000000
    //Delete	Length:2	Total weight:0.929195	Forward Strand Weight:0.929195	Reverse Strand Weight:0.000000
    //9	G total-weight:0.464145	total-pos-weight:0.464145	total-neg-weight:0.000000	G:0.464145 (1.000000) +str:1.000000, -str:nan,	Total-weight:0.464145
    //10	G total-weight:0.059048	total-pos-weight:0.059048	total-neg-weight:0.000000	G:0.059048 (1.000000) +str:1.000000, -str:nan,	Total-weight:0.059048
    //11	G total-weight:0.476354	total-pos-weight:0.476354	total-neg-weight:0.000000	G:0.476354 (1.000000) +str:1.000000, -str:nan,	Total-weight:0.476354

    st_logInfo("Read:\t\t%s\n", read->rleRead->rleString);
    st_logInfo("Reference:\t%s\n", reference->rleString);
    if (st_getLogLevel() >= info) {
        poa_print(poa, stderr, reads, 0.0);
    }

    // Check inserts
    checkInserts(testCase, poa, 0, 0, (const char *[]) {""}, (const double[]) {1}, 1);
    checkInserts(testCase, poa, 1, 0, (const char *[]) {""}, (const double[]) {1}, 1);
    checkInserts(testCase, poa, 2, 1, (const char *[]) {"T"}, (const double[]) {0.990160}, 1);
    checkInserts(testCase, poa, 3, 0, (const char *[]) {""}, (const double[]) {1}, 1);
    checkInserts(testCase, poa, 4, 0, (const char *[]) {""}, (const double[]) {1}, 1);
    checkInserts(testCase, poa, 5, 0, (const char *[]) {""}, (const double[]) {1}, 1);
    checkInserts(testCase, poa, 6, 0, (const char *[]) {""}, (const double[]) {1}, 1);
    checkInserts(testCase, poa, 7, 0, (const char *[]) {""}, (const double[]) {1}, 1);
    checkInserts(testCase, poa, 8, 0, (const char *[]) {""}, (const double[]) {1}, 1);
    checkInserts(testCase, poa, 9, 0, (const char *[]) {""}, (const double[]) {1}, 1);
    checkInserts(testCase, poa, 10, 0, (const char *[]) {""}, (const double[]) {1}, 1);
    checkInserts(testCase, poa, 11, 0, (const char *[]) {""}, (const double[]) {1}, 1);

    // Check deletes

    // No deletes first three positions
    checkDeletes(testCase, poa, 0, 0, (const int64_t[]) {1}, (const double[]) {1}, 1);
    checkDeletes(testCase, poa, 1, 0, (const int64_t[]) {1}, (const double[]) {1}, 1);
    checkDeletes(testCase, poa, 2, 0, (const int64_t[]) {1}, (const double[]) {1}, 1);
    checkDeletes(testCase, poa, 3, 0, (const int64_t[]) {1}, (const double[]) {1}, 1);
    checkDeletes(testCase, poa, 4, 0, (const int64_t[]) {1}, (const double[]) {1}, 1);
    checkDeletes(testCase, poa, 5, 0, (const int64_t[]) {1}, (const double[]) {1}, 1);
    checkDeletes(testCase, poa, 6, 0, (const int64_t[]) {1}, (const double[]) {1}, 1);
    checkDeletes(testCase, poa, 7, 0, (const int64_t[]) {1}, (const double[]) {1}, 1);
    checkDeletes(testCase, poa, 8, 2, (const int64_t[]) {1, 2}, (const double[]) {1.045831, 0.929195}, 1);
    checkDeletes(testCase, poa, 9, 0, (const int64_t[]) {1}, (const double[]) {1}, 1);
    checkDeletes(testCase, poa, 10, 0, (const int64_t[]) {1}, (const double[]) {1}, 1);
    checkDeletes(testCase, poa, 11, 0, (const int64_t[]) {1}, (const double[]) {1}, 1);

    params_destruct(params);
    poa_destruct(poa);
    stList_destruct(reads);
    rleString_destruct(reference);
}

static void test_poa_realign(CuTest *testCase) {
    /*
     * Test poa_realign by generating lots of random examples
     */

    for (int64_t test = 0; test < 100; test++) {

        Params *params = params_readParams(polishParamsFile);
        PolishParams *polishParams = params->polishParams;

        //Make true reference
        char *trueReference = getRandomSequence(st_randomInt(1, 100));

        // Make starting reference
        char *reference = evolveSequence(trueReference);
        RleString *reference_rle = params->polishParams->useRunLengthEncoding ?
                                   rleString_construct(reference) : rleString_construct_no_rle(reference);

        // Reads
        int64_t readNumber = st_randomInt(0, 20);
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        for (int64_t i = 0; i < readNumber; i++) {
            stList_append(reads,
                          bamChunkRead_construct2(stString_print("read_%d", i), evolveSequence(trueReference), NULL,
                                                  TRUE,
                                                  params->polishParams->useRunLengthEncoding));
        }

        Poa *poa = poa_realign(reads, NULL, reference_rle, polishParams);

        // Generate the read alignments and check the matches
        // Currently don't check the insert and deletes

        double *baseWeights = st_calloc(poa->alphabet->alphabetSize * reference_rle->length, sizeof(double));
        double *repeatCountWeights = st_calloc(poa->maxRepeatCount * reference_rle->length, sizeof(double));

        for (int64_t i = 0; i < readNumber; i++) {
            BamChunkRead *bamChunkRead = stList_get(reads, i);
            RleString *read = bamChunkRead->rleRead;

            // Make symbol strings
            SymbolString sX = rleString_constructSymbolString(reference_rle, 0, reference_rle->length, poa->alphabet,
                                                              polishParams->useRepeatCountsInAlignment, MAXIMUM_REPEAT_LENGTH);
            SymbolString sY = rleString_constructSymbolString(read, 0, read->length, poa->alphabet,
                                                              polishParams->useRepeatCountsInAlignment, MAXIMUM_REPEAT_LENGTH);

            // Generate set of posterior probabilities for matches, deletes and inserts with respect to reference.
            stList *matches = NULL, *inserts = NULL, *deletes = NULL;
            getAlignedPairsWithIndels(bamChunkRead->forwardStrand ?
                                      polishParams->stateMachineForForwardStrandRead :
                                      polishParams->stateMachineForReverseStrandRead, sX, sY, polishParams->p, &matches,
                                      &deletes, &inserts, 0, 0);

            // Collate matches
            for (int64_t j = 0; j < stList_length(matches); j++) {
                stIntTuple *match = stList_get(matches, j);
                baseWeights[stIntTuple_get(match, 1) * poa->alphabet->alphabetSize + poa->alphabet->convertCharToSymbol(
                        read->rleString[stIntTuple_get(match, 2)])] += stIntTuple_get(match, 0);
                repeatCountWeights[stIntTuple_get(match, 1) * poa->maxRepeatCount +
                                   read->repeatCounts[stIntTuple_get(match, 2)]] += stIntTuple_get(match, 0);
            }

            // Cleanup
            stList_destruct(matches);
            stList_destruct(inserts);
            stList_destruct(deletes);
            symbolString_destruct(sX);
            symbolString_destruct(sY);
        }

        // Check match and repeat count weights tally
        for (int64_t i = 0; i < reference_rle->length; i++) {
            PoaNode *poaNode = stList_get(poa->nodes, i + 1);
            for (int64_t j = 0; j < poa->alphabet->alphabetSize; j++) {
                CuAssertDblEquals(testCase, poaNode->baseWeights[j], baseWeights[i * poa->alphabet->alphabetSize + j],
                                  0.0001);
            }
            for (int64_t j = 0; j < poa->maxRepeatCount; j++) {
                CuAssertDblEquals(testCase, poaNode->repeatCountWeights[j],
                                  repeatCountWeights[i * poa->maxRepeatCount + j], 0.0001);
            }
        }

        st_logInfo("True-reference:%s\n", trueReference);
        if (st_getLogLevel() >= info) {
            poa_print(poa, stderr, reads, 5);
        }

        //Cleanup
        free(baseWeights);
        free(repeatCountWeights);
        free(trueReference);
        free(reference);
        rleString_destruct(reference_rle);
        stList_destruct(reads);
        poa_destruct(poa);
        params_destruct(params);
    }
}

static void test_poa_realignIterative(CuTest *testCase) {
    /*
     * Test random small examples against poa_realignIterative
     */

    for (int64_t test = 0; test < 100; test++) {
        Params *params = params_readParams(polishParamsFile);
        PolishParams *polishParams = params->polishParams;

        //Make true reference
        char *trueReference = getRandomSequence(st_randomInt(1, 100));

        // Make starting reference
        char *reference = evolveSequence(trueReference);
        RleString *reference_rle = params->polishParams->useRunLengthEncoding ?
                                   rleString_construct(reference) : rleString_construct_no_rle(reference);

        // Reads
        int64_t readNumber = st_randomInt(0, 20);
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        for (int64_t i = 0; i < readNumber; i++) {
            stList_append(reads, bamChunkRead_construct2(stString_print("Read_%d", i), evolveSequence(trueReference),
                                                         NULL, st_random() > 0.5,
                                                         params->polishParams->useRunLengthEncoding));
        }

        Poa *poa = poa_realignAll(reads, NULL, reference_rle, params->polishParams);

        st_logInfo("True-reference:%s\n", trueReference);
        if (st_getLogLevel() >= info) {
            poa_print(poa, stderr, reads, 5);
        }

        //Cleanup
        free(trueReference);
        free(reference);
        stList_destruct(reads);
        poa_destruct(poa);
        params_destruct(params);
        rleString_destruct(reference_rle);
    }
}

int64_t calcSequenceMatches(char *seq1, char *seq2) {
    Params *params = params_readParams(polishParamsFile);
    PolishParams *polishParams = params->polishParams;

    // Load non-rle statemachine
    StateMachine *sM = stateMachine3_constructNucleotide(threeState);

    // Make symbol strings
    SymbolString sX = symbolString_construct(seq1, 0, strlen(seq1), params->polishParams->alphabet);
    SymbolString sY = symbolString_construct(seq2, 0, strlen(seq2), params->polishParams->alphabet);

    //Get identity
    stList *allAlignedPairs = getAlignedPairs(sM, sX, sY, polishParams->p, 0, 0);
    stList *alignedPairs = filterPairwiseAlignmentToMakePairsOrdered(allAlignedPairs, sX, sY, params->polishParams->p);

    int64_t matches = getNumberOfMatchingAlignedPairs(sX, sY, alignedPairs);


    // Cleanup
    params_destruct(params);
    symbolString_destruct(sX);
    symbolString_destruct(sY);
    stateMachine_destruct(sM);
    //stList_destruct(alignedPairs);

    return matches;
}

typedef struct _alignmentMetrics {
    int64_t totalConsensusMatches;
    int64_t totalReferenceMatches;
    int64_t totalConsensusLength;
    int64_t totalReferenceLength;
    int64_t totalTrueReferenceLength;
} AlignmentMetrics;


static void test_poa_realign_example(CuTest *testCase, char *trueReference, char *reference,
                                     stList *originalReads, AlignmentMetrics *rleAlignmentMetrics,
                                     AlignmentMetrics *nonRleAlignmentMetrics,
                                     bool rle) {
    stList *reads = stList_construct();
    for (int64_t i = 0; i < stList_length(originalReads); i++) {
        BamChunkRead *bcr = stList_get(originalReads, i);
        stList_append(reads, bamChunkRead_constructCopy(bcr));
    }
    RleString *rleReference = rle ? rleString_construct(reference) : rleString_construct_no_rle(reference);
    RleString *rleTrueReference = rle ? rleString_construct(trueReference) : rleString_construct_no_rle(trueReference);

    Params *params = params_readParams(polishParamsFile);
    PolishParams *polishParams = params->polishParams;
    params->polishParams->useRunLengthEncoding = rle;
    // Set parameters
    params->polishParams->maxPoaConsensusIterations = 100;
    params->polishParams->minPoaConsensusIterations = 3;
    params->polishParams->referenceBasePenalty = 0.6;
    params->polishParams->maxRealignmentPolishIterations = 3;
    params->polishParams->minRealignmentPolishIterations = 3;

    Poa *poa = poa_realign(reads, NULL, rleReference, params->polishParams);
    Poa *poaRefined = poa_realignAll(reads, NULL, rleReference, params->polishParams);

    Poa *poaTrue = poa_realign(reads, NULL, rleTrueReference, params->polishParams);

    // Run phasing
    //stList *anchorAlignments = poa_getAnchorAlignments(poaRefined, NULL, stList_length(reads), params->polishParams);
    //stList *reads1, *reads2;
    //phaseReads(poaRefined->refString, strlen(poaRefined->refString), reads, anchorAlignments, &reads1, &reads2, params);
    //Poa *poaReads1 = poa_realignIterative(reads1, NULL, poaRefined->refString, polishParams);
    //Poa *poaReads2 = poa_realignIterative(reads2, NULL, poaRefined->refString, polishParams);

    // Look at non-rle comparison
    poa_estimateRepeatCountsUsingBayesianModel(poaRefined, reads, params->polishParams->repeatSubMatrix);
    char *nonRLEConsensusString = rleString_expand(poaRefined->refString);

    // Calculate alignments between true reference and consensus and starting reference sequences
    int64_t consensusMatches = calcSequenceMatches(rleTrueReference->rleString, poaRefined->refString->rleString);
    int64_t referenceMatches = calcSequenceMatches(rleTrueReference->rleString, rleReference->rleString);
    int64_t nonRLEConsensusMatches = calcSequenceMatches(trueReference, nonRLEConsensusString);
    int64_t nonRLEReferenceMatches = calcSequenceMatches(trueReference, reference);
    //int64_t consensusMatchesReads1 = calcSequenceMatches(rleTrueReference->rleString, poaReads1->refString);
    //int64_t consensusMatchesReads2 = calcSequenceMatches(rleTrueReference->rleString, poaReads2->refString);

    // Update the running total alignment metrics
    if (rleAlignmentMetrics != NULL) {
        rleAlignmentMetrics->totalConsensusMatches += consensusMatches;
        rleAlignmentMetrics->totalReferenceMatches += referenceMatches;
        rleAlignmentMetrics->totalConsensusLength += poaRefined->refString->length;
        rleAlignmentMetrics->totalReferenceLength += rleReference->length;
        rleAlignmentMetrics->totalTrueReferenceLength += rleTrueReference->length;
    }

    // Update the running total alignment metrics
    if (nonRleAlignmentMetrics != NULL) {
        nonRleAlignmentMetrics->totalConsensusMatches += nonRLEConsensusMatches;
        nonRleAlignmentMetrics->totalReferenceMatches += nonRLEReferenceMatches;
        nonRleAlignmentMetrics->totalConsensusLength += strlen(nonRLEConsensusString);
        nonRleAlignmentMetrics->totalReferenceLength += strlen(reference);
        nonRleAlignmentMetrics->totalTrueReferenceLength += strlen(trueReference);
    }

    // Log some stuff
    if (st_getLogLevel() >= info) {
        st_logInfo("Reference:      \t\t%s\n", rleReference->rleString);
        st_logInfo("True-reference: \t\t%s\n", rleTrueReference->rleString);
        st_logInfo("Consensus:      \t\t%s\n", poaRefined->refString->rleString);
        //st_logInfo("Consensus Reads1:\t%s\n", poaReads1->refString);
        //st_logInfo("Consensus Reads2:\t%s\n", poaReads2->refString);
        st_logInfo("Reference stats:     \t");
        poa_printSummaryStats(poa, stderr);
        st_logInfo("Consensus stats:     \t");
        poa_printSummaryStats(poaRefined, stderr);
        //st_logInfo("Reads 1 stats\t");
        //poa_printSummaryStats(poaReads1, stderr);
        //st_logInfo("Reads 2 stats\t");
        //poa_printSummaryStats(poaReads2, stderr);
        st_logInfo("True-reference stats:\t");
        poa_printSummaryStats(poaTrue, stderr);
        st_logInfo("Consensus : true-ref identity: %f\n",
                   2.0 * consensusMatches / (rleTrueReference->length + poaRefined->refString->length));
        //st_logInfo("Reads 1 consensus : true-ref identity: %f\n", 2.0*consensusMatchesReads1/(rleTrueReference->length + strlen(poaReads1->refString)));
        //st_logInfo("Reads 2 consensus : true-ref identity: %f\n", 2.0*consensusMatchesReads2/(rleTrueReference->length + strlen(poaReads2->refString)));
        st_logInfo("Start-ref : true-ref identity: %f\n",
                   2.0 * referenceMatches / (rleTrueReference->length + rleReference->length));
        //st_logInfo("Total reads: %i, # reads partition1: %i, # reads partition2: %i\n", (int)stList_length(reads), (int)stList_length(reads1), (int)stList_length(reads2));
        // Non-RLE stats
        st_logInfo("Non-RLE Reference:     \t\t%s\n", reference);
        st_logInfo("Non-RLE True-reference:\t\t%s\n", trueReference);
        st_logInfo("Non-RLE Consensus:     \t\t%s\n", nonRLEConsensusString);
        st_logInfo("Non-RLE Consensus : true-ref identity: %f\n",
                   2.0 * nonRLEConsensusMatches / (strlen(trueReference) + strlen(nonRLEConsensusString)));
        st_logInfo("Non-RLE Start-ref : true-ref identity: %f\n\n",
                   2.0 * nonRLEReferenceMatches / (strlen(trueReference) + strlen(reference)));
    }

    if (st_getLogLevel() >= debug && !stString_eq(rleTrueReference->rleString, poaRefined->refString->rleString)) {
        //poa_print(poa, stderr, 5);
        //poa_print(poaRefined, stderr, reads, 2, 0);

        //poa_printTSV(poa, stderr, reads, 2, 0);

        ///poa_printRepeatCounts(poa, stderr, reads);
    }

    // Cleanup
    params_destruct(params);
    poa_destruct(poa);
    poa_destruct(poaRefined);
    poa_destruct(poaTrue);
//	free(trueReference);
//	free(reference);
    //rleString_destruct(rleTrueReference);
    //rleString_destruct(rleReference);
    free(nonRLEConsensusString);
    //poa_destruct(poaReads1);
    //poa_destruct(poaReads2);
    //stList_destruct(anchorAlignments);
    //stList_destruct(reads1);
    //stList_destruct(reads2);
}

static struct List *readSequences(char *fastaFile, struct List **headers) {
    struct List *seqs = constructEmptyList(0, free);
    struct List *seqLengths = constructEmptyList(0, free);
    *headers = constructEmptyList(0, free);

    FILE *fH = fopen(fastaFile, "r");
    fastaRead(fH, seqs, seqLengths, *headers);
    fclose(fH);

    destructList(seqLengths);

    return seqs;
}

static void test_poa_realign_examples(CuTest *testCase, const char **examples, int64_t exampleNo, bool rle) {
    AlignmentMetrics *alignmentMetrics = st_calloc(1, sizeof(AlignmentMetrics));
    AlignmentMetrics *rleAlignmentMetrics = st_calloc(1, sizeof(AlignmentMetrics));
    for (int64_t example = 0; example < exampleNo; example++) {
        const char *readFile = examples[example * 2];
        const char *trueRefFile = examples[example * 2 + 1];

        st_logInfo("\nDoing polish test with %s read files and %s true ref file\n", readFile, trueRefFile);

        // Parse sequences
        struct List *readHeaders;
        struct List *nucleotides = readSequences((char *) readFile, &readHeaders);
        assert(nucleotides->length > 1);
        struct List *trueReferenceHeaders;
        struct List *trueReferenceList = readSequences((char *) trueRefFile, &trueReferenceHeaders);
        assert(trueReferenceList->length == 1);

        // Parse strands (and create reads)
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        for (int64_t i = 1; i < readHeaders->length; i++) {
            char *header = readHeaders->list[i];
            char strand = header[strlen(header) - 1];
            CuAssertTrue(testCase, strand == 'F' || strand == 'R');
            stList_append(reads, bamChunkRead_construct2(stString_print("read_%d", i),
                                                         stString_copy(nucleotides->list[i]), NULL, strand == 'F',
                                                         rle));
        }

        //if(strlen(reads->list[0]) > strlen(trueReferenceList->list[0]) * 0.8 || reads->length < 30) {
        //	fprintf(stderr, "Got short input ref:\n%s\n%s\n", reads->list[0], trueReferenceList->list[0]);
        //	continue;
        //}

        // Run poa iterative realign
        test_poa_realign_example(testCase, trueReferenceList->list[0], nucleotides->list[0], reads,
                                 rleAlignmentMetrics, alignmentMetrics, rle);

        // Cleanup
        destructList(nucleotides);
        destructList(readHeaders);
        destructList(trueReferenceList);
        destructList(trueReferenceHeaders);
        free(reads);
    }

    // Alignment metrics for set
    st_logInfo("Total consensus identity: %f for %i bases, vs. %f starting identity\n",
               2.0 * alignmentMetrics->totalConsensusMatches /
               (alignmentMetrics->totalConsensusLength + alignmentMetrics->totalTrueReferenceLength),
               alignmentMetrics->totalConsensusLength,
               2.0 * alignmentMetrics->totalReferenceMatches /
               (alignmentMetrics->totalReferenceLength + alignmentMetrics->totalTrueReferenceLength));
    if (rle) {
        st_logInfo("RLE Space: Total consensus identity: %f for %i bases, vs. %f starting identity\n",
                   2.0 * rleAlignmentMetrics->totalConsensusMatches /
                   (rleAlignmentMetrics->totalConsensusLength + rleAlignmentMetrics->totalTrueReferenceLength),
                   rleAlignmentMetrics->totalConsensusLength,
                   2.0 * rleAlignmentMetrics->totalReferenceMatches /
                   (rleAlignmentMetrics->totalReferenceLength + rleAlignmentMetrics->totalTrueReferenceLength));
    }

    free(alignmentMetrics);
    free(rleAlignmentMetrics);
}

static void test_poa_realign_examples_large(CuTest *testCase, int64_t exampleNo, const char *path, bool rle) {
    // Build strings
    const char **examples = st_calloc(2 * exampleNo, sizeof(char *));
    for (int64_t i = 0; i < exampleNo; i++) {
        examples[2 * i] = stString_print("%s/%i.fasta", path, (int) i);
        examples[2 * i + 1] = stString_print("%s/%i.ref.fasta", path, (int) i);
    }

    test_poa_realign_examples(testCase, examples, exampleNo, rle);

    // Cleanup
    for (int64_t i = 0; i < exampleNo * 2; i++) {
        free((char *) examples[i]);
    }
    free(examples);
}

void test_poa_realign_ecoli_examples_rle(CuTest *testCase) {
    test_poa_realign_examples_large(testCase, 20,
                                    TEST_POLISH_FILES_DIR"20_random_100bp_windows_directional_ecoli_guppy", 1);
}

void test_poa_realign_ecoli_examples_no_rle(CuTest *testCase) {
    test_poa_realign_examples_large(testCase, 20,
                                    TEST_POLISH_FILES_DIR"20_random_100bp_windows_directional_ecoli_guppy", 0);
}

void test_poa_realign_ecoli_many_examples_rle(CuTest *testCase) {
    test_poa_realign_examples_large(testCase, 100,
                                    TEST_POLISH_FILES_DIR"500_random_100bp_windows_directional_ecoli_guppy", 1);
}

void test_poa_realign_ecoli_many_examples_no_rle(CuTest *testCase) {
    test_poa_realign_examples_large(testCase, 100,
                                    TEST_POLISH_FILES_DIR"500_random_100bp_windows_directional_ecoli_guppy", 0);
}

static void test_rleString_example(CuTest *testCase, const char *testStr,
                                   int64_t rleLength, int64_t nonRleLength,
                                   const char *testStrRLE, const int64_t *repeatCounts,
                                   const int64_t *nonRleToRleCoordinateMap) {
    RleString *rleString = rleString_construct((char *) testStr);
    uint64_t *nonRleToRleCoordinateMap2 = rleString_getNonRleToRleCoordinateMap(rleString);

    CuAssertIntEquals(testCase, rleLength, rleString->length);
    CuAssertStrEquals(testCase, testStrRLE, rleString->rleString);
    for (int64_t i = 0; i < rleLength; i++) {
        CuAssertIntEquals(testCase, repeatCounts[i], rleString->repeatCounts[i]);
    }

    CuAssertIntEquals(testCase, nonRleLength, rleString->nonRleLength);
    for (int64_t i = 0; i < nonRleLength; i++) {
        CuAssertIntEquals(testCase, nonRleToRleCoordinateMap[i], nonRleToRleCoordinateMap2[i]);
    }

    char *expandedRleString = rleString_expand(rleString);
    CuAssertStrEquals(testCase, testStr, expandedRleString);

    free(expandedRleString);
    rleString_destruct(rleString);
    free(nonRleToRleCoordinateMap2);
}

static void test_rleString_examples(CuTest *testCase) {
    test_rleString_example(testCase, "GATTACAGGGGTT", 8, 13, "GATACAGT", (const int64_t[]) {1, 1, 2, 1, 1, 1, 4, 2},
                           (const int64_t[]) {0, 1, 2, 2, 3, 4, 5, 6, 6, 6, 6, 7, 7});

    test_rleString_example(testCase, "TTTTT", 1, 5, "T", (const int64_t[]) {5},
                           (const int64_t[]) {0, 0, 0, 0, 0});

    test_rleString_example(testCase, "", 0, 0, "", (const int64_t[]) {1},
                           (const int64_t[]) {0});

    test_rleString_example(testCase, "TTTTTCC", 2, 7, "TC", (const int64_t[]) {5, 2},
                           (const int64_t[]) {0, 0, 0, 0, 0, 1, 1});
}

void test_rle_rotateString(CuTest *testCase) {
    // Specific example
    RleString *e = rleString_construct("GATAACA");
    rleString_rotateString(e, 2, 1);
    RleString *e_rotated = rleString_construct("CAGATAA");
    CuAssertTrue(testCase, rleString_eq(e, e_rotated));
    rleString_destruct(e);
    rleString_destruct(e_rotated);

    // Another example, but where the resulting rle strings are shorter because the rotation merges the first and last position
    // of the initial string
    e = rleString_construct("ATAA");
    rleString_rotateString(e, 1, 1);
    e_rotated = rleString_construct("AAAT");
    CuAssertTrue(testCase, rleString_eq(e, e_rotated));
    rleString_destruct(e);
    rleString_destruct(e_rotated);

    for (int64_t test = 0; test < 500; test++) {
        // Generate random string and make rotation of it
        char *s = getRandomSequence(st_randomInt(0, 20));

        RleString *t_rotated = rleString_construct(s);
        RleString *t = rleString_construct(s);
        int64_t i = st_randomInt(0, 20);
        rleString_rotateString(t_rotated, i, 0);

        // Check the result
        for (int64_t j = 0; j < t->length; j++) {
            CuAssertIntEquals(testCase, t->rleString[j], t_rotated->rleString[(j + i) % t->length]);
            CuAssertIntEquals(testCase, t->repeatCounts[j], t_rotated->repeatCounts[(j + i) % t->length]);
        }

        // Cleanup
        rleString_destruct(t);
        rleString_destruct(t_rotated);
        free(s);
    }
}

void checkStringsAndFree(CuTest *testCase, const char *expected, char *temp) {
    CuAssertStrEquals(testCase, expected, temp);
    free(temp);
}

void test_addInsert(CuTest *testCase) {
    checkStringsAndFree(testCase, "GATTACA", addInsert("GAACA", "TT", 2));
    checkStringsAndFree(testCase, "GATTACA", addInsert("", "GATTACA", 0));
    checkStringsAndFree(testCase, "GATTACA", addInsert("ATTACA", "G", 0));
    checkStringsAndFree(testCase, "GATTACA", addInsert("GATTAC", "A", 6));
    checkStringsAndFree(testCase, "GATTACA", addInsert("GATTACA", "", 6));
    checkStringsAndFree(testCase, "GATTACA", addInsert("GATTACA", "", 3));
}

void test_removeDelete(CuTest *testCase) {
    checkStringsAndFree(testCase, "GATTACA", removeDelete("GATTGGACA", 2, 4));
    checkStringsAndFree(testCase, "GATTACA", removeDelete("GATTACA", 0, 0));
    checkStringsAndFree(testCase, "GATTACA", removeDelete("GATTACATT", 2, 7));
    checkStringsAndFree(testCase, "GATTACA", removeDelete("AGATTACA", 1, 0));
}

void test_polishParams(CuTest *testCase) {
    Params *params = params_readParams(polishParamsFile);
    PolishParams *polishParams = params->polishParams;

    CuAssertTrue(testCase, polishParams->useRunLengthEncoding);
    CuAssertDblEquals(testCase, polishParams->referenceBasePenalty, 0.5, 0);
    CuAssertDblEquals(testCase, polishParams->minPosteriorProbForAlignmentAnchorsLength, 6, 0);

    CuAssertDblEquals(testCase, polishParams->minPosteriorProbForAlignmentAnchors[0], 0.9, 0);
    CuAssertDblEquals(testCase, polishParams->minPosteriorProbForAlignmentAnchors[1], 10, 0);
    CuAssertDblEquals(testCase, polishParams->minPosteriorProbForAlignmentAnchors[2], 0.95, 0);
    CuAssertDblEquals(testCase, polishParams->minPosteriorProbForAlignmentAnchors[3], 4, 0);
    CuAssertDblEquals(testCase, polishParams->minPosteriorProbForAlignmentAnchors[4], 0.99, 0);
    CuAssertDblEquals(testCase, polishParams->minPosteriorProbForAlignmentAnchors[5], 0, 0);

    CuAssertDblEquals(testCase, polishParams->p->threshold, 0.01, 0);
    CuAssertDblEquals(testCase, polishParams->p->minDiagsBetweenTraceBack, 5000, 0);
    CuAssertDblEquals(testCase, polishParams->p->traceBackDiagonals, 40, 0);
    CuAssertDblEquals(testCase, polishParams->p->diagonalExpansion, 4, 0);
    CuAssertDblEquals(testCase, polishParams->p->constraintDiagonalTrim, 0, 0);
    CuAssertDblEquals(testCase, polishParams->p->splitMatrixBiggerThanThis, 250000000, 0);
    CuAssertDblEquals(testCase, polishParams->p->gapGamma, 0.25, 0);
    CuAssertTrue(testCase, !polishParams->p->alignAmbiguityCharacters);

    /*CuAssertDblEquals(testCase,  repeatSubMatrix_getLogProb(polishParams->repeatSubMatrix, a, 0, 0, 0), -0.059686935, 0);
    CuAssertDblEquals(testCase,  repeatSubMatrix_getLogProb(polishParams->repeatSubMatrix, c, 0, 0, 0), -0.055418707, 0);
    CuAssertDblEquals(testCase,  repeatSubMatrix_getLogProb(polishParams->repeatSubMatrix, g, 0, 0, 0), -0.05438334, 0);
    CuAssertDblEquals(testCase,  repeatSubMatrix_getLogProb(polishParams->repeatSubMatrix, t, 0, 0, 0), -0.035762809, 0);
    CuAssertDblEquals(testCase,  repeatSubMatrix_getLogProb(polishParams->repeatSubMatrix, a, 1, 0, 0), -0.036856437, 0);
    CuAssertDblEquals(testCase,  repeatSubMatrix_getLogProb(polishParams->repeatSubMatrix, c, 1, 0, 0), -0.062816805, 0);
    CuAssertDblEquals(testCase,  repeatSubMatrix_getLogProb(polishParams->repeatSubMatrix, g, 1, 0, 0), -0.055853556, 0);
    CuAssertDblEquals(testCase,  repeatSubMatrix_getLogProb(polishParams->repeatSubMatrix, t, 1, 0, 0), -0.065273937, 0);*/

    params_destruct(params);
}

void test_removeOverlapExample(CuTest *testCase) {
    Params *params = params_readParams(polishParamsFile);
    PolishParams *polishParams = params->polishParams;

    //Make prefix
    char *prefixString = stString_copy("ACGTACGTACGTACGTACGTACGTGATTTCAACGTACGT");

    // Make sufix
    char *suffixString = stString_copy("ACGTACGTGATTTCAACGTACGTACGTACGTACGTACGT");

    int64_t approxOverlap = 25;

    // Run overlap remover
    int64_t prefixStringCropEnd, suffixStringCropStart;
    double overlapWeight = removeOverlap(prefixString, strlen(prefixString), suffixString, strlen(suffixString),
                                         approxOverlap, polishParams,
                                         &prefixStringCropEnd, &suffixStringCropStart);

    char *pSeqCropped = stString_getSubString(prefixString, 0, prefixStringCropEnd);
    char *seqCropped = stString_getSubString(suffixString, suffixStringCropStart,
                                             strlen(suffixString) - suffixStringCropStart);
    char *stitchedSeq = stString_print("%s%s", pSeqCropped, seqCropped);

    CuAssertTrue(testCase, stString_eqcase(stitchedSeq, "ACGTACGTACGTACGTACGTACGTGATTTCAACGTACGTACGTACGTACGTACGT"));

    // Cleanup
    params_destruct(params);
    free(prefixString);
    free(suffixString);
    free(stitchedSeq);
}



void test_removeOverlap_RandomExamples(CuTest *testCase) {
    Params *params = params_readParams(polishParamsFile);
    PolishParams *polishParams = params->polishParams;
    params->polishParams->useRunLengthEncoding = FALSE;

    for (int64_t test = 0; test < 100; test++) {
        // Truth
        char *truth = getRandomSequence(st_randomInt(200, 300));
        int64_t halfway = strlen(truth) / 2;
        int64_t boundarySize = 30;

        //Make prefix
        int64_t prefixStop = st_randomInt(halfway + boundarySize/2, halfway + boundarySize*3/2);
        char *prefixString = stString_getSubString(truth, 0, prefixStop);

        // Make sufix
        int64_t suffixStart = st_randomInt(halfway - boundarySize*3/2, halfway - boundarySize/2);
        char *suffixString = stString_getSubString(truth, suffixStart, strlen(truth) - suffixStart);

        // Run overlap remover
        int64_t prefixStringCropEnd, suffixStringCropStart;
        removeOverlap(prefixString, strlen(prefixString), suffixString, strlen(suffixString), boundarySize*2, polishParams,
                      &prefixStringCropEnd, &suffixStringCropStart);

        CuAssertTrue(testCase, prefixStringCropEnd >= 0);
        CuAssertTrue(testCase, prefixStringCropEnd <= strlen(prefixString));

        CuAssertTrue(testCase, suffixStringCropStart >= 0);
        CuAssertTrue(testCase, suffixStringCropStart <= strlen(suffixString));

        char *pSeqCropped = stString_getSubString(prefixString, 0, prefixStringCropEnd);
        char *seqCropped = stString_getSubString(suffixString, suffixStringCropStart,
                strlen(suffixString) - suffixStringCropStart);
        char *stitchedSeq = stString_print("%s%s", pSeqCropped, seqCropped);

        CuAssertTrue(testCase, stString_eqcase(stitchedSeq, truth));

        free(pSeqCropped);
        free(seqCropped);
        free(prefixString);
        free(suffixString);
        free(stitchedSeq);
        free(truth);
    }

    params_destruct(params);
}


int64_t polishingTest(char *bamFile, char *referenceFile, char *paramsFile, char *region, bool verbose, bool diploid) {

    // Run margin phase
    char *logString = verbose ? "--logLevel DEBUG" : "--logLevel INFO";
    char *regionStr = region == NULL ? stString_print("") : stString_print("--region %s", region);
    char *diploidString = diploid ? "--diploid" : "";
    char *command = stString_print("./marginPolish %s %s %s %s %s %s", bamFile, referenceFile, paramsFile, regionStr,
                                   logString, diploidString);
    st_logInfo("> Running command: %s\n", command);

    int64_t i = st_system(command);
    free(regionStr);
    free(command);
    return i;
}

void test_polish5kb_rle(CuTest *testCase) {
    char *referenceFile = "../tests/data/realData/hg19.chr3.9mb.fa";
    bool verbose = false;
    char *bamFile = "../tests/data/realData/NA12878.np.chr3.5kb.bam";
    char *region = "chr3:2150000-2155000";

    st_logInfo("\n\nTesting polishing on %s\n", bamFile);
    int64_t i = polishingTest(bamFile, referenceFile, polishParamsFile, region, verbose, FALSE);
    CuAssertTrue(testCase, i == 0);

    st_logInfo("\n\nTesting diploid polishing on %s\n", bamFile);
    i = polishingTest(bamFile, referenceFile, polishParamsFile, region, verbose, TRUE);
    CuAssertTrue(testCase, i == 0);
}

void checkLargeGapOutput(CuTest *testCase) {
    //read output file, find non-n sequence
    char *outputFile = "output.fa";
    FILE *fh = fopen(outputFile, "r");
    stHash *referenceSequences = fastaReadToMap(fh);  //valgrind says blocks from this allocation are "still reachable"
    fclose(fh);
    char *polished = stHash_search(referenceSequences, "chr3");
    stList *splitByNs = stString_splitByString(polished, "N");
    CuAssertTrue(testCase, stList_length(splitByNs) > 2);

    for (int64_t i = 0; i < stList_length(splitByNs); i++) {
        int64_t partLength = strlen(stList_get(splitByNs, i));
        CuAssertTrue(testCase, partLength == 0 || (partLength > 1900 && partLength < 2100));
    }

    stList_destruct(splitByNs);
    stHash_destruct(referenceSequences);
}

void test_largeGap(CuTest *testCase) {
    char *referenceFile = "../tests/data/realData/hg19.chr3.9mb.fa";
    bool verbose = false;
    char *bamFile = "../tests/data/largeGapTest/largeGapTest.bam";
    char *region = "chr3:10000-17000";

    st_logInfo("\n\nTesting polishing on %s\n", bamFile);
    int64_t i = polishingTest(bamFile, referenceFile, polishParamsFile, region, verbose, FALSE);
    CuAssertTrue(testCase, i == 0);

    checkLargeGapOutput(testCase);
}

void test_largeGap2(CuTest *testCase) {
    char *referenceFile = "../tests/data/realData/hg19.chr3.9mb.fa";
    bool verbose = false;
    char *bamFile = "../tests/data/largeGapTest/largeGapTest2.bam";
    char *region = "chr3:10000-17000";

    st_logInfo("\n\nTesting polishing on %s\n", bamFile);
    int64_t i = polishingTest(bamFile, referenceFile, polishParamsFile, region, verbose, FALSE);
    CuAssertTrue(testCase, i == 0);

    checkLargeGapOutput(testCase);
}

void test_binomialPValue(CuTest *testCase) {
    CuAssertDblEquals(testCase, 252.0, bionomialCoefficient(10, 5), 0.001);
    CuAssertDblEquals(testCase, 15504.0, bionomialCoefficient(20, 15), 0.001);
    CuAssertDblEquals(testCase, 80347448443237920.0, bionomialCoefficient(64, 22), 0.001);
    CuAssertDblEquals(testCase, 151473214816.0, bionomialCoefficient(64, 10), 0.001);
    CuAssertDblEquals(testCase, 1832624140942590534.0, bionomialCoefficient(64, 32), 0.001);
}

CuSuite *polisherTestSuite(void) {
    CuSuite *suite = CuSuiteNew();

//    SUITE_ADD_TEST(suite, test_poa_getReferenceGraph);
//    SUITE_ADD_TEST(suite, test_getShift);
//    SUITE_ADD_TEST(suite, test_rleString_examples);
//    SUITE_ADD_TEST(suite, test_rle_rotateString);
//    SUITE_ADD_TEST(suite, test_poa_augment_example);
//    SUITE_ADD_TEST(suite, test_poa_realign_tiny_example1);
//    SUITE_ADD_TEST(suite, test_poa_realign);
//    SUITE_ADD_TEST(suite, test_getShift);
//    SUITE_ADD_TEST(suite, test_rleString_examples);
//    SUITE_ADD_TEST(suite, test_addInsert);
//    SUITE_ADD_TEST(suite, test_removeDelete);
//    SUITE_ADD_TEST(suite, test_polishParams);
    SUITE_ADD_TEST(suite, test_removeOverlapExample);
    SUITE_ADD_TEST(suite, test_removeOverlap_RandomExamples);
//    SUITE_ADD_TEST(suite, test_binomialPValue);
//    SUITE_ADD_TEST(suite, test_poa_realignIterative);
//    SUITE_ADD_TEST(suite, test_poa_realign_ecoli_examples_rle);
//    SUITE_ADD_TEST(suite, test_poa_realign_ecoli_examples_no_rle);
//    SUITE_ADD_TEST(suite, test_poa_realign_ecoli_many_examples_rle);
//    SUITE_ADD_TEST(suite, test_poa_realign_ecoli_many_examples_no_rle);
//    SUITE_ADD_TEST(suite, test_polish5kb_rle);
//    SUITE_ADD_TEST(suite, test_largeGap);
//    SUITE_ADD_TEST(suite, test_largeGap2);

    return suite;
}
