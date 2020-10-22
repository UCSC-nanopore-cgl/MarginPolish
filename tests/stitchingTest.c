//
// Created by Benedict Paten on 3/26/20.
//

#include "CuTest.h"
#include "margin.h"

static char *paramsFile = "../params/ont/r9.4/allParams.np.human.r94-g344.json";
static char *outputSequenceFile = "./testStitchingSequenceFile.fa";
static char *outputPoaFile = "./testStitchingPoaFile.csv";
static char *outputRepeatCountFile = "./testStitchingRepeatCountFile.csv";
static char *noRleParamsFile = "../params/misc/allParams.no_rle.json";


void checkCSV(CuTest *testCase, char *csvFile, char *sequence) {
    /*
     * Checks the CSV file we output is organized as expected
     */
    stList *csvLines = stFile_getLinesFromFile(csvFile);
    for (int64_t i = 1; i < stList_length(csvLines); i++) {
        stList *tokens = stString_splitByString(stList_get(csvLines, i), ",");
        // Check the sequence of positions is contiguous
        int64_t j = strtol(stList_get(tokens, 0), NULL, 10);
        CuAssertIntEquals(testCase, i - 1, j);
        // Check we have the expected reference base
        if (i > 1) {
            char *refBase = stList_get(tokens, 1);
            CuAssertTrue(testCase, strlen(refBase) == 1);
            CuAssertTrue(testCase, toupper(refBase[0]) == toupper(sequence[i - 2]));
        }
    }
    CuAssertIntEquals(testCase, strlen(sequence) + 2, stList_length(csvLines));
    stList_destruct(csvLines);
}

char *getSequence(CuTest *testCase, char *outputSequenceFile, char *sequenceName) {
    FILE *fh = fopen(outputSequenceFile, "r");
    stHash *seqMap = fastaReadToMap(fh);
    fclose(fh);
    stList *seqNames = stHash_getKeys(seqMap);
    stList *seqs = stHash_getValues(seqMap);
    CuAssertTrue(testCase, stList_length(seqNames) == 1);
    CuAssertStrEquals(testCase, sequenceName, stList_peek(seqNames));
    char *sequence = stString_copy(stList_pop(seqs));
    // Cleanup
    stList_destruct(seqs);
    stList_destruct(seqNames);
    stHash_destruct(seqMap);

    return sequence;
}

void test_stitching(CuTest *testCase) {
    /*
     * Runs the stitcher with a set of chunks and checks we get back the original sequence
     */
    setPairwiseAlignerKmerSize(2);
    setMinOverlapAnchorPairs(1);

    for (int64_t test = 0; test < 10; test++) {
        // Get params
        Params *params = params_readParams(paramsFile);
        params->polishParams->useRunLengthEncoding = FALSE; // Turn off RLE for this test to work
        params->polishParams->chunkBoundary = 3;

        // Sequences
        char *sequence = "AAAAAAAAAATTTTTTTTTTCCCCCCCCCCGGGGGGGGGG";
        char *sequenceName = "seq1";
        char *chunk0 = "AAAA";
        char *chunk1 = "AAAAAAAAAAT";
        char *chunk2 = "AAATTT";
        char *chunk3 = "AAATTTTTTTTTTCCCCC";
        char *chunk4 = "TTTCCCCCCCCCCG";
        char *chunk5 = "CGGGGGGGGGG";
        char *chunk6 = ""; // Test a 0 length chunk :)
        stList *chunks = stList_construct();
        stList_append(chunks, chunk0);
        stList_append(chunks, chunk1);
        stList_append(chunks, chunk2);
        stList_append(chunks, chunk3);
        stList_append(chunks, chunk4);
        stList_append(chunks, chunk5);
        stList_append(chunks, chunk6);
        stList *randomizedChunks = stList_copy(chunks, NULL);
        stList_shuffle(randomizedChunks); // This randomizes the order of the chunks
        int64_t noOfOutputChunkers = st_randomInt(1, 5);

        // Get chunker
        OutputChunkers *outputChunkers = outputChunkers_construct(noOfOutputChunkers, params,
                                                                  outputSequenceFile, outputPoaFile, NULL,
                                                                  outputRepeatCountFile,
                                                                  NULL, NULL, 0);

        // Now process the chunks
        for (int64_t i = 0; i < stList_length(randomizedChunks); i++) {
            char *chunk = stList_get(randomizedChunks, i);
            int64_t chunkOrdinal = stList_find(chunks, chunk);
            // Build the inputs for the chunker
            RleString *rle_chunk = rleString_construct_no_rle(chunk);
            Poa *poa = poa_getReferenceGraph(rle_chunk, params->polishParams->alphabet,
                                             (uint64_t ) params->polishParams->repeatSubMatrix->maximumRepeatLength);
            stList *reads = stList_construct(); // No reads, currently
            // Output the chunk
            outputChunkers_processChunkSequence(outputChunkers, st_randomInt(0, noOfOutputChunkers), chunkOrdinal,
                                                sequenceName, poa, reads);
            // Cleanup
            poa_destruct(poa);
            stList_destruct(reads);
            rleString_destruct(rle_chunk);
        }

        // Do stitching
        outputChunkers_stitch(outputChunkers, 0, stList_length(randomizedChunks));

        // Destroy chunkers
        outputChunkers_destruct(outputChunkers);

        // Read output sequence and check it agrees with what we expect
        char *seqFromFile = getSequence(testCase, outputSequenceFile, sequenceName);
        CuAssertStrEquals(testCase, sequence, seqFromFile);

        // Check poa
        checkCSV(testCase, outputPoaFile, sequence);

        // Check repeat count file
        checkCSV(testCase, outputRepeatCountFile, sequence);

        // Cleanup
        stList_destruct(randomizedChunks);
        stFile_rmrf(outputSequenceFile);
        stFile_rmrf(outputPoaFile);
        stFile_rmrf(outputRepeatCountFile);
        params_destruct(params);
        stList_destruct(chunks);
    }
}


ChunkToStitch **getChunksToStitchFromStrings(char **strings, int len) {
    ChunkToStitch **chunks = st_calloc(len, sizeof(ChunkToStitch*));
    for (int i = 0; i < len ; i++) {
        chunks[i] = chunkToStitch_construct(stString_copy("TestContig"),i,FALSE, FALSE, FALSE);
        chunks[i]->seqHap1 = stString_copy(strings[i]);
    }
    return chunks;
}

void test_mergeContigChunks(CuTest *testCase) {
    Params *params = params_readParams(noRleParamsFile);
    params->polishParams->chunkBoundary = 16;
    char **chunks = st_calloc(4, sizeof(char *));
    chunks[0] = stString_copy("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCC");
    chunks[1] = stString_copy("AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGG");
    chunks[2] = stString_copy("CCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT");
    chunks[3] = stString_copy("GGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
    ChunkToStitch **chunksToStitch = getChunksToStitchFromStrings(chunks, 4);
    ChunkToStitch *result = mergeContigChunkz(chunksToStitch, 0, 4, FALSE, params);
    CuAssertStrEquals(testCase, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", result->seqHap1);
}


void test_mergeContigChunksThreaded(CuTest *testCase) {
    Params *params = params_readParams(noRleParamsFile);
    //TODO needs to be 16
    params->polishParams->chunkBoundary = 2;
    setPairwiseAlignerKmerSize(2);
    char *chunks[16];// = st_calloc(16, sizeof(char*));
    chunks[0] = stString_copy("AAAAAAAACC");
    chunks[1] = stString_copy("AACCCCCCCCGG");
    chunks[2] = stString_copy("CCGGGGGGGGTT");
    chunks[3] = stString_copy("GGTTTTTTTTAA");
    chunks[4] = stString_copy("TTAAAAAAAACC");
    chunks[5] = stString_copy("AACCCCCCCCGG");
    chunks[6] = stString_copy("CCGGGGGGGGTT");
    chunks[7] = stString_copy("GGTTTTTTTTAA");
    chunks[8] = stString_copy("TTAAAAAAAACC");
    chunks[9] = stString_copy("AACCCCCCCCGG");
    chunks[10] = stString_copy("CCGGGGGGGGTT");
    chunks[11] = stString_copy("GGTTTTTTTTAA");
    chunks[12] = stString_copy("TTAAAAAAAACC");
    chunks[13] = stString_copy("AACCCCCCCCGG");
    chunks[14] = stString_copy("CCGGGGGGGGTT");
    chunks[15] = stString_copy("GGTTTTTTTT");
    char *truth = "AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTTAAAAAAAACCCCCCCCGGGGGGGGTTTTTTTTAAAAAAAACCCCCCCCGGGGGGGGTTTTTTTTAAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT";
    char* contig;

    ChunkToStitch **chunksToStitch = getChunksToStitchFromStrings(chunks, 16);
    contig = mergeContigChunkzThreaded(chunksToStitch, 0, 16, 1, FALSE, params, "testContig")->seqHap1;
    CuAssertTrue(testCase, strcmp(contig, truth) == 0);
    chunksToStitch = getChunksToStitchFromStrings(chunks, 16);
    contig = mergeContigChunkzThreaded(chunksToStitch, 0, 16, 2, FALSE, params, "testContig")->seqHap1;
    CuAssertTrue(testCase, strcmp(contig, truth) == 0);
    chunksToStitch = getChunksToStitchFromStrings(chunks, 16);
    contig = mergeContigChunkzThreaded(chunksToStitch, 0, 16, 3, FALSE, params, "testContig")->seqHap1;
    CuAssertTrue(testCase, strcmp(contig, truth) == 0);
    chunksToStitch = getChunksToStitchFromStrings(chunks, 16);
    contig = mergeContigChunkzThreaded(chunksToStitch, 0, 16, 4, FALSE, params, "testContig")->seqHap1;
    CuAssertTrue(testCase, strcmp(contig, truth) == 0);
    chunksToStitch = getChunksToStitchFromStrings(chunks, 16);
    contig = mergeContigChunkzThreaded(chunksToStitch, 0, 16, 5, FALSE, params, "testContig")->seqHap1;
    CuAssertTrue(testCase, strcmp(contig, truth) == 0);
    chunksToStitch = getChunksToStitchFromStrings(chunks, 16);
    contig = mergeContigChunkzThreaded(chunksToStitch, 0, 16, 6, FALSE, params, "testContig")->seqHap1;
    CuAssertTrue(testCase, strcmp(contig, truth) == 0);
    chunksToStitch = getChunksToStitchFromStrings(chunks, 16);
    contig = mergeContigChunkzThreaded(chunksToStitch, 0, 16, 7, FALSE, params, "testContig")->seqHap1;
    CuAssertTrue(testCase, strcmp(contig, truth) == 0 );
    chunksToStitch = getChunksToStitchFromStrings(chunks, 16);
    contig = mergeContigChunkzThreaded(chunksToStitch, 0, 16, 8, FALSE, params, "testContig")->seqHap1;
    CuAssertTrue(testCase, strcmp(contig, truth) == 0 );
}


CuSuite *stitchingTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
//    SUITE_ADD_TEST(suite, test_stitching);
    SUITE_ADD_TEST(suite, test_mergeContigChunks);
//    SUITE_ADD_TEST(suite, test_mergeContigChunksThreaded);
    return suite;
}
