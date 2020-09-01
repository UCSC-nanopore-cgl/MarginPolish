/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <htsIntegration.h>
#include "CuTest.h"
#include "margin.h"

#define INPUT_BAM "../tests/data/chunkingTest/chunkingTest.bam"
#define INPUT_PARAMS "../params/misc/allParams.no_rle.json"

static PolishParams *getParameters(uint64_t chunkSize, uint64_t chunkBoundary, bool includeSoftClipping) {
    PolishParams *params = polishParams_constructEmpty();
    params->chunkSize = chunkSize;
    params->chunkBoundary = chunkBoundary;
    params->includeSoftClipping = includeSoftClipping;
    params->useRunLengthEncoding = 0;

    return params;
}

static void test_getRegionChunker(CuTest *testCase) {
    // test part of contig (less than chunk size)
    PolishParams *params = getParameters(0, 0, FALSE);
    BamChunker *chunker = bamChunker_construct2(INPUT_BAM, "contig_1:100000-110000", params, false);
    CuAssertTrue(testCase, chunker->chunkCount == 1);
    CuAssertTrue(testCase, stString_eq(((BamChunk *) stList_get(chunker->chunks, 0))->refSeqName, "contig_1"));
    CuAssertTrue(testCase, ((BamChunk *) stList_get(chunker->chunks, 0))->chunkOverlapStart == 100000);
    CuAssertTrue(testCase, ((BamChunk *) stList_get(chunker->chunks, 0))->chunkOverlapEnd == 100008);

    // test whole contig by region
    chunker = bamChunker_construct2(INPUT_BAM, "contig_1:0-3000000", params, false);
    CuAssertTrue(testCase, chunker->chunkCount == 1);
    CuAssertTrue(testCase, stString_eq(((BamChunk *) stList_get(chunker->chunks, 0))->refSeqName, "contig_1"));
    CuAssertTrue(testCase, ((BamChunk *) stList_get(chunker->chunks, 0))->chunkOverlapStart == 100000);
    CuAssertTrue(testCase, ((BamChunk *) stList_get(chunker->chunks, 0))->chunkOverlapEnd == 2100008);
    free(chunker->params);

    params = getParameters(100000, 0, FALSE);
    chunker = bamChunker_construct2(INPUT_BAM, "contig_1:100000-300000", params, false);
    CuAssertTrue(testCase, chunker->chunkCount == 2);
    CuAssertTrue(testCase, stString_eq(((BamChunk *) stList_get(chunker->chunks, 0))->refSeqName, "contig_1"));
    CuAssertTrue(testCase, ((BamChunk *) stList_get(chunker->chunks, 0))->chunkOverlapStart == 100000);
    CuAssertTrue(testCase, ((BamChunk *) stList_get(chunker->chunks, 0))->chunkOverlapEnd == 200000);
    CuAssertTrue(testCase, stString_eq(((BamChunk *) stList_get(chunker->chunks, 1))->refSeqName, "contig_1"));
    CuAssertTrue(testCase, ((BamChunk *) stList_get(chunker->chunks, 1))->chunkOverlapStart == 200000);
    CuAssertTrue(testCase, ((BamChunk *) stList_get(chunker->chunks, 1))->chunkOverlapEnd ==
                           210020); //end pos stops at last aligned pos
    free(chunker->params);

    bamChunker_destruct(chunker);
}

static void test_getChunksByChrom(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct(INPUT_BAM, getParameters(0, 0, FALSE));
    CuAssertTrue(testCase, chunker->chunkCount == 2);
    free(chunker->params);
    bamChunker_destruct(chunker);
}

static void test_getChunksBy100kb(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct(INPUT_BAM, getParameters(100000, 0, FALSE));

    // contig_1 alignments start at 100 000 and go to 2 100 008 (21 @ 100k size)
    // contig_2 alignments start at 100 000 and go to   100 032 ( 1 @ 100k size)
    CuAssertTrue(testCase, chunker->chunkCount == 22);
    free(chunker->params);
    bamChunker_destruct(chunker);
}

static void test_getQualityScores(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct(INPUT_BAM, getParameters(10000, 0, FALSE));

    // has 9 reads of 8 characters aligned to position 100 000 and every 4 bases after (last read aligned to 100 032)
    BamChunk *chunk = NULL;

    for (int64_t chunkIdx = 0; chunkIdx < chunker->chunkCount; chunkIdx++) {
        chunk = bamChunker_getChunk(chunker, chunkIdx);
        // margin testing is on contig_2
        if (!stString_eq(chunk->refSeqName, "contig_2")) continue;

        // see how many reads there are
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *alignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
        uint32_t readCount = convertToReadsAndAlignments(chunk, NULL, reads, alignments, chunker->params);
        CuAssertTrue(testCase, readCount == stList_length(reads));
        CuAssertTrue(testCase, readCount == 9);

        for (int64_t i = 0; i < 9; i++) {
            BamChunkRead *read = stList_get(reads, i);
            switch (i) {
                case 0:
                    CuAssertTrue(testCase, read->qualities != NULL);
                    for (uint8_t j = 0; j < 8; j++) {
                        CuAssertTrue(testCase, read->qualities[j] == 15 + j);
                    }
                    break;
                case 1:
                    CuAssertTrue(testCase, read->qualities != NULL);
                    for (uint8_t j = 0; j < 8; j++) {
                        CuAssertTrue(testCase, read->qualities[j] == 22 - j);
                    }
                    break;
                case 2:
                    CuAssertTrue(testCase, read->qualities != NULL);
                    for (uint8_t j = 0; j < 8; j++) {
                        CuAssertTrue(testCase, read->qualities[j] == 32 + j);
                    }
                    break;
                case 3:
                    CuAssertTrue(testCase, read->qualities != NULL);
                    for (uint8_t j = 0; j < 8; j++) {
                        CuAssertTrue(testCase, read->qualities[j] == 0);
                    }
                    break;
                case 4:
                    CuAssertTrue(testCase, read->qualities != NULL);
                    for (uint8_t j = 0; j < 8; j++) {
                        CuAssertTrue(testCase, read->qualities[j] == 9);
                    }
                    break;
                case 5:
                case 6:
                case 7:
                case 8:
                    CuAssertTrue(testCase, read->qualities == NULL);
                    break;
                default:
                    CuAssertTrue(testCase, FALSE);
            }
        }


        // increment
        stList_destruct(reads);
        stList_destruct(alignments);
    }
    free(chunker->params);
    bamChunker_destruct(chunker);

}

static void test_getChunksWithBoundary(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct(INPUT_BAM, getParameters(8, 4, FALSE));

    // has 9 reads of 8 characters aligned to position 100 000 and every 4 bases after (last read aligned to 100 032)
    int contig2ChunkCount = 0;
    BamChunk *chunk = NULL;
    for (int64_t chunkIdx = 0; chunkIdx < chunker->chunkCount; chunkIdx++) {
        chunk = bamChunker_getChunk(chunker, chunkIdx);
        // margin testing is on contig_2
        if (!stString_eq(chunk->refSeqName, "contig_2")) continue;

        // see how many reads there are
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *alignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
        uint32_t readCount = convertToReadsAndAlignments(chunk, NULL, reads, alignments, chunker->params);
        CuAssertTrue(testCase, readCount == stList_length(reads));

        // all these should cover any read whose alignment overlaps (inclusive) with chunkMarginStart and (exclusive)
        // with chunkMarginEnd.  exclusive as in, if an alignment starts at 0 and has length 100, it should be included
        // if chunkMarginEnd is 99 but not 100
        if (contig2ChunkCount == 0) {
            //0: cbs:99996 cs:100000 ce:100008 cbe:100012 -> reads 0, 4, 8 ; !12+
            CuAssertTrue(testCase, readCount == 3);
        } else if (contig2ChunkCount == 1) {
            //1: cbs:100004 cs:100008 ce:1000016 cbe:100020 -> reads 0, 4, 8, 12, 16 ; !20+
            CuAssertTrue(testCase, readCount == 5);
        } else if (contig2ChunkCount == 2) {
            //2: cbs:100012 cs:100016 ce:1000024 cbe:100028 -> reads !4- ; 8, 12, 16, 20, 24 ; !28+
            CuAssertTrue(testCase, readCount == 5);
        } else if (contig2ChunkCount == 3) {
            //3: cbs:100020 cs:100024 ce:1000032 cbe:100036 -> reads !12- ; 16, 20, 24, 28, 32
            CuAssertTrue(testCase, readCount == 5);
        } else if (contig2ChunkCount == 4) {
            //4: cbs:100028 cs:100032 ce:1000040 cbe:100044 -> reads !20- ; 24, 28, 32
            CuAssertTrue(testCase, readCount == 3);
        }

        // increment
        contig2ChunkCount++;
        stList_destruct(reads);
        stList_destruct(alignments);
    }
    CuAssertTrue(testCase, contig2ChunkCount == 5);
    free(chunker->params);
    bamChunker_destruct(chunker);
}


static void test_getChunksWithoutBoundary(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct(INPUT_BAM, getParameters(8, 0, FALSE));

    // has 9 reads of 8 characters aligned to position 100 000 and every 4 bases after (last read aligned to 100 032)
    int contig2ChunkCount = 0;
    BamChunk *chunk = NULL;
    for (int64_t chunkIdx = 0; chunkIdx < chunker->chunkCount; chunkIdx++) {
        chunk = bamChunker_getChunk(chunker, chunkIdx);
        // margin testing is on contig_2
        if (!stString_eq(chunk->refSeqName, "contig_2")) continue;

        // see how many reads there are
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *alignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
        uint32_t readCount = convertToReadsAndAlignments(chunk, NULL, reads, alignments, chunker->params);
        CuAssertTrue(testCase, readCount == stList_length(reads));

        // all these should cover any read whose alignment overlaps (inclusive) with chunkMarginStart and (exclusive)
        // with chunkMarginEnd.  exclusive as in, if an alignment starts at 0 and has length 100, it should be included
        // if chunkMarginEnd is 99 but not 100
        if (contig2ChunkCount == 0) {
            //0: cbs:100000 cs:100000 ce:100008 cbe:100008 -> reads 0, 4 ; !8+
            CuAssertTrue(testCase, readCount == 2);
        } else if (contig2ChunkCount == 1) {
            //1: cbs:100008 cs:100008 ce:1000016 cbe:1000016 -> reads !0 ; 4, 8, 12 ; !16+
            CuAssertTrue(testCase, readCount == 3);
        } else if (contig2ChunkCount == 2) {
            //2: cbs:100016 cs:100016 ce:1000024 cbe:1000024 -> reads !8- ; 12, 16, 20 ; !24+
            CuAssertTrue(testCase, readCount == 3);
        } else if (contig2ChunkCount == 3) {
            //3: cbs:100024 cs:100024 ce:1000032 cbe:1000032 -> reads !16- ; 20, 24, 28 ; !32
            CuAssertTrue(testCase, readCount == 3);
        } else if (contig2ChunkCount == 4) {
            //4: cbs:100032 cs:100032 ce:1000040 cbe:1000040 -> reads !24- ; 28, 32
            CuAssertTrue(testCase, readCount == 2);
        }

        // increment
        contig2ChunkCount++;
        stList_destruct(reads);
        stList_destruct(alignments);
    }
    CuAssertTrue(testCase, contig2ChunkCount == 5);
    free(chunker->params);
    bamChunker_destruct(chunker);
}

void assertClippingAlignmentMatchCount(CuTest *testCase, int64_t idx, stList *alignment) {
    switch (idx) {
        case 0: //8S8M
        case 1: //8M8S
        case 2: //4S8M4S
        case 4: //4S4M2D4M4S
        case 6: //4S1M1D6M1D1M4S
        case 7: //4H8S8M
        case 8: //8M8S4H
        case 9: //4H4S8M4S4H
            CuAssertTrue(testCase, stList_length(alignment) == 8);
            break;
        case 3: //4S2M4I2M4S
            CuAssertTrue(testCase, stList_length(alignment) == 4);
            break;
        case 5: //4S1M1I4M1I1M4S
            CuAssertTrue(testCase, stList_length(alignment) == 6);
            break;
        default:
            CuAssertTrue(testCase, FALSE);
    }
}

static void test_getReadsWithoutSoftClipping(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct(INPUT_BAM, getParameters(100000, 0, FALSE));

    // have reads aligned to a region wholly between 200 000 and 299 999 for soft clip testing
    BamChunk *chunk = NULL;
    bool foundChunk = FALSE;
    for (int64_t chunkIdx = 0; chunkIdx < chunker->chunkCount; chunkIdx++) {
        chunk = bamChunker_getChunk(chunker, chunkIdx);
        if (!stString_eq(chunk->refSeqName, "contig_1") || !(chunk->chunkOverlapStart == 200000 &&
                                                             chunk->chunkOverlapEnd == 300000))
            continue;
        CuAssertTrue(testCase, !foundChunk);
        foundChunk = TRUE;

        // analyze reads and alignments
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *alignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
        uint32_t readCount = convertToReadsAndAlignments(chunk, NULL, reads, alignments, chunker->params);
        CuAssertTrue(testCase, readCount == 10);
        for (int64_t i = 0; i < 10; i++) {
            // check the length of the cigar strings
            assertClippingAlignmentMatchCount(testCase, i, stList_get(alignments, i));
        }
        stList_destruct(reads);
        stList_destruct(alignments);
    }
    CuAssertTrue(testCase, foundChunk);
    free(chunker->params);
    bamChunker_destruct(chunker);
}

static void test_getReadsWithSoftClipping(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct(INPUT_BAM, getParameters(100000, 0, TRUE));

    // have reads aligned to a region wholly between 200 000 and 299 999 for soft clip testing
    BamChunk *chunk = NULL;
    bool foundChunk = FALSE;
    for (int64_t chunkIdx = 0; chunkIdx < chunker->chunkCount; chunkIdx++) {
        chunk = bamChunker_getChunk(chunker, chunkIdx);
        if (!stString_eq(chunk->refSeqName, "contig_1") || !(chunk->chunkOverlapStart == 200000 &&
                                                             chunk->chunkOverlapEnd == 300000))
            continue;
        CuAssertTrue(testCase, !foundChunk);
        foundChunk = TRUE;

        // analyze reads and alignments
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *alignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
        uint32_t readCount = convertToReadsAndAlignments(chunk, NULL, reads, alignments, chunker->params);
        CuAssertTrue(testCase, readCount == 10);
        for (int64_t i = 0; i < 10; i++) {
            // check the length of the cigar strings
            assertClippingAlignmentMatchCount(testCase, i, stList_get(alignments, i));
        }
        stList_destruct(reads);
    }
    CuAssertTrue(testCase, foundChunk);
    free(chunker->params);
    bamChunker_destruct(chunker);
}

void assertAlignmentMatching(CuTest *testCase, stList *list1, int64_t *onethValues, int64_t *zerothValues, int len) {
    CuAssertTrue(testCase, stList_length(list1) == len);
    for (int i = 0; i < len; i++) {
        stIntTuple *tup = (stIntTuple *) stList_get(list1, i);
        CuAssertTrue(testCase, stIntTuple_get(tup, 0) == zerothValues[i]);
        CuAssertTrue(testCase, stIntTuple_get(tup, 1) == onethValues[i]);
    }
}

static void test_readAlignmentsWithoutSoftclippingChunkStart(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct(INPUT_BAM, getParameters(1000, 0, FALSE));

    BamChunk *chunk = NULL;
    bool foundChunk = FALSE;
    for (int64_t chunkIdx = 0; chunkIdx < chunker->chunkCount; chunkIdx++) {
        chunk = bamChunker_getChunk(chunker, chunkIdx);
        if (!stString_eq(chunk->refSeqName, "contig_1") || !(chunk->chunkOverlapStart == 400000 &&
                                                             chunk->chunkOverlapEnd == 401000))
            continue;
        CuAssertTrue(testCase, !foundChunk);
        foundChunk = TRUE;

        // analyze reads and alignments
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *alignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
        uint32_t readCount = convertToReadsAndAlignments(chunk, NULL, reads, alignments, chunker->params);
        CuAssertTrue(testCase, readCount == 24);
        for (int64_t i = 0; i < 24; i++) {
            BamChunkRead *read = stList_get(reads, i);
            switch (i) {
                case 0: //399996      8M        ACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3}, (int64_t[]) {0, 1, 2, 3}, 4);
                    break;
                case 1: //399996      4M1D3M    ACGTCGT
                    CuAssertTrue(testCase, stString_eq("CGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2}, (int64_t[]) {1, 2, 3}, 3);
                    break;
                case 2: //399996      4M1I4M    ACGTAACGT
                    CuAssertTrue(testCase, stString_eq("AACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {1, 2, 3, 4}, (int64_t[]) {0, 1, 2, 3}, 4);
                    break;
                case 3: //399996      4S8M      AAAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3}, (int64_t[]) {0, 1, 2, 3}, 4);
                    break;
                case 4: //399996      4S4M1D3M  AAAAACGTCGT
                    CuAssertTrue(testCase, stString_eq("CGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2}, (int64_t[]) {1, 2, 3}, 3);
                    break;
                case 5: //399996      4S4M1I4M  AAAAACGTAACGT
                    CuAssertTrue(testCase, stString_eq("AACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {1, 2, 3, 4}, (int64_t[]) {0, 1, 2, 3}, 4);
                    break;
                case 6: //400000      8M        ACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7}, (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            8);
                    break;
                case 7: //400000      1D7M      CGTACGT
                    CuAssertTrue(testCase, stString_eq("CGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6}, (int64_t[]) {1, 2, 3, 4, 5, 6, 7}, 7);
                    break;
                case 8: //400000      1I8M      AACGTACGT
                    CuAssertTrue(testCase, stString_eq("AACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {1, 2, 3, 4, 5, 6, 7, 8}, (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            8);
                    break;
                case 9: //400000      4S8M      AAAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7}, (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            8);
                    break;
                case 10: //400000      4S1D7M    AAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("CGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6}, (int64_t[]) {1, 2, 3, 4, 5, 6, 7}, 7);
                    break;
                case 11: //400000      4S1I8M    AAAAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("AACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {1, 2, 3, 4, 5, 6, 7, 8}, (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            8);
                    break;
                case 12: //400002      8M        ACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7}, (int64_t[]) {2, 3, 4, 5, 6, 7, 8, 9},
                                            8);
                    break;
                case 13: //400002      1D7M      CGTACGT
                    CuAssertTrue(testCase, stString_eq("CGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6}, (int64_t[]) {3, 4, 5, 6, 7, 8, 9}, 7);
                    break;
                case 14: //400002      1I8M      AACGTACGT
                    CuAssertTrue(testCase, stString_eq("AACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {1, 2, 3, 4, 5, 6, 7, 8}, (int64_t[]) {2, 3, 4, 5, 6, 7, 8, 9},
                                            8);
                    break;
                case 15: //400002      4S8M      AAAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7}, (int64_t[]) {2, 3, 4, 5, 6, 7, 8, 9},
                                            8);
                    break;
                case 16: //400002      4S1D7M    AAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("CGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6}, (int64_t[]) {3, 4, 5, 6, 7, 8, 9}, 7);
                    break;
                case 17: //400002      4S1I8M    AAAAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("AACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {1, 2, 3, 4, 5, 6, 7, 8}, (int64_t[]) {2, 3, 4, 5, 6, 7, 8, 9},
                                            8);
                    break;
                case 18: //400008  8M      ACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            (int64_t[]) {8, 9, 10, 11, 12, 13, 14, 15}, 8);
                    break;
                case 19: //400008  1D7M    CGTACGT
                    CuAssertTrue(testCase, stString_eq("CGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6}, (int64_t[]) {9, 10, 11, 12, 13, 14, 15},
                                            7);
                    break;
                case 20: //400008  1I8M    AACGTACGT
                    CuAssertTrue(testCase, stString_eq("AACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {1, 2, 3, 4, 5, 6, 7, 8},
                                            (int64_t[]) {8, 9, 10, 11, 12, 13, 14, 15}, 8);
                    break;
                case 21: //400008  4S8M    AAAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            (int64_t[]) {8, 9, 10, 11, 12, 13, 14, 15}, 8);
                    break;
                case 22: //400008  4S1D7M  AAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("CGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6}, (int64_t[]) {9, 10, 11, 12, 13, 14, 15},
                                            7);
                    break;
                case 23: //400008  4S1I8M  AAAAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("AACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {1, 2, 3, 4, 5, 6, 7, 8},
                                            (int64_t[]) {8, 9, 10, 11, 12, 13, 14, 15}, 8);
                    break;
                default:
                    CuAssertTrue(testCase, FALSE);
            }
        }
        stList_destruct(reads);
        stList_destruct(alignments);
    }
    CuAssertTrue(testCase, foundChunk);
    free(chunker->params);
    bamChunker_destruct(chunker);
}


static void test_readAlignmentsWithSoftclippingChunkStart(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct(INPUT_BAM, getParameters(1000, 0, TRUE));

    // have reads aligned to a region wholly between 200 000 and 299 999 for soft clip testing
    BamChunk *chunk = NULL;
    bool foundChunk = FALSE;
    for (int64_t chunkIdx = 0; chunkIdx < chunker->chunkCount; chunkIdx++) {
        chunk = bamChunker_getChunk(chunker, chunkIdx);
        if (!stString_eq(chunk->refSeqName, "contig_1") || !(chunk->chunkOverlapStart == 400000 &&
                                                             chunk->chunkOverlapEnd == 401000))
            continue;
        CuAssertTrue(testCase, !foundChunk);
        foundChunk = TRUE;

        // analyze reads and alignments
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *alignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
        uint32_t readCount = convertToReadsAndAlignments(chunk, NULL, reads, alignments, chunker->params);
        CuAssertTrue(testCase, readCount == 24);
        for (int64_t i = 0; i < 24; i++) {
            BamChunkRead *read = stList_get(reads, i);
            switch (i) {
                case 0: //399996      8M        ACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3}, (int64_t[]) {0, 1, 2, 3}, 4);
                    break;
                case 1: //399996      4M1D3M    ACGTCGT
                    CuAssertTrue(testCase, stString_eq("CGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2}, (int64_t[]) {1, 2, 3}, 3);
                    break;
                case 2: //399996      4M1I4M    ACGTAACGT
                    CuAssertTrue(testCase, stString_eq("AACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {1, 2, 3, 4}, (int64_t[]) {0, 1, 2, 3}, 4);
                    break;
                case 3: //399996      4S8M      AAAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3}, (int64_t[]) {0, 1, 2, 3}, 4);
                    break;
                case 4: //399996      4S4M1D3M  AAAAACGTCGT
                    CuAssertTrue(testCase, stString_eq("CGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2}, (int64_t[]) {1, 2, 3}, 3);
                    break;
                case 5: //399996      4S4M1I4M  AAAAACGTAACGT
                    CuAssertTrue(testCase, stString_eq("AACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {1, 2, 3, 4}, (int64_t[]) {0, 1, 2, 3}, 4);
                    break;
                case 6: //400000      8M        ACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7}, (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            8);
                    break;
                case 7: //400000      1D7M      CGTACGT
                    CuAssertTrue(testCase, stString_eq("CGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6}, (int64_t[]) {1, 2, 3, 4, 5, 6, 7}, 7);
                    break;
                case 8: //400000      1I8M      AACGTACGT
                    CuAssertTrue(testCase, stString_eq("AACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {1, 2, 3, 4, 5, 6, 7, 8}, (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            8);
                    break;
                case 9: //400000      4S8M      AAAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7}, (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            8);
                    break;
                case 10: //400000      4S1D7M    AAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("CGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6}, (int64_t[]) {1, 2, 3, 4, 5, 6, 7}, 7);
                    break;
                case 11: //400000      4S1I8M    AAAAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("AACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {1, 2, 3, 4, 5, 6, 7, 8}, (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            8);
                    break;
                case 12: //400002      8M        ACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7}, (int64_t[]) {2, 3, 4, 5, 6, 7, 8, 9},
                                            8);
                    break;
                case 13: //400002      1D7M      CGTACGT
                    CuAssertTrue(testCase, stString_eq("CGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6}, (int64_t[]) {3, 4, 5, 6, 7, 8, 9}, 7);
                    break;
                case 14: //400002      1I8M      AACGTACGT
                    CuAssertTrue(testCase, stString_eq("AACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {1, 2, 3, 4, 5, 6, 7, 8}, (int64_t[]) {2, 3, 4, 5, 6, 7, 8, 9},
                                            8);
                    break;
                case 15: //400002      4S8M      AAAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("AAACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {2, 3, 4, 5, 6, 7, 8, 9}, (int64_t[]) {2, 3, 4, 5, 6, 7, 8, 9},
                                            8);
                    break;
                case 16: //400002      4S1D7M    AAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("AACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {2, 3, 4, 5, 6, 7, 8}, (int64_t[]) {3, 4, 5, 6, 7, 8, 9}, 7);
                    break;
                case 17: //400002      4S1I8M    AAAAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("AAAACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {3, 4, 5, 6, 7, 8, 9, 10}, (int64_t[]) {2, 3, 4, 5, 6, 7, 8, 9},
                                            8);
                    break;
                case 18: //400008  8M      ACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            (int64_t[]) {8, 9, 10, 11, 12, 13, 14, 15}, 8);
                    break;
                case 19: //400008  1D7M    CGTACGT
                    CuAssertTrue(testCase, stString_eq("CGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6}, (int64_t[]) {9, 10, 11, 12, 13, 14, 15},
                                            7);
                    break;
                case 20: //400008  1I8M    AACGTACGT
                    CuAssertTrue(testCase, stString_eq("AACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {1, 2, 3, 4, 5, 6, 7, 8},
                                            (int64_t[]) {8, 9, 10, 11, 12, 13, 14, 15}, 8);
                    break;
                case 21: //400008  4S8M    AAAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("AAAAACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {4, 5, 6, 7, 8, 9, 10, 11},
                                            (int64_t[]) {8, 9, 10, 11, 12, 13, 14, 15}, 8);
                    break;
                case 22: //400008  4S1D7M  AAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("AAAACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {4, 5, 6, 7, 8, 9, 10}, (int64_t[]) {9, 10, 11, 12, 13, 14, 15},
                                            7);
                    break;
                case 23: //400008  4S1I8M  AAAAAACGTACGT
                    CuAssertTrue(testCase, stString_eq("AAAAAACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {5, 6, 7, 8, 9, 10, 11, 12},
                                            (int64_t[]) {8, 9, 10, 11, 12, 13, 14, 15}, 8);
                    break;
                default:
                    CuAssertTrue(testCase, FALSE);
            }
        }
        stList_destruct(reads);
        stList_destruct(alignments);
    }
    CuAssertTrue(testCase, foundChunk);
    free(chunker->params);
    bamChunker_destruct(chunker);
}


static void test_readAlignmentsWithoutSoftclippingChunkEnd(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct(INPUT_BAM, getParameters(20, 0, FALSE));

    BamChunk *chunk = NULL;
    bool foundChunk = FALSE;
    for (int64_t chunkIdx = 0; chunkIdx < chunker->chunkCount; chunkIdx++) {
        chunk = bamChunker_getChunk(chunker, chunkIdx);
        if (!stString_eq(chunk->refSeqName, "contig_1") || !(chunk->chunkOverlapStart == 410000 &&
                                                             chunk->chunkOverlapEnd == 410020))
            continue;
        CuAssertTrue(testCase, !foundChunk);
        foundChunk = TRUE;

        // analyze reads and alignments
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *alignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
        uint32_t readCount = convertToReadsAndAlignments(chunk, NULL, reads, alignments, chunker->params);
        CuAssertTrue(testCase, readCount == 21);
        for (int64_t i = 0; i < 21; i++) {
            BamChunkRead *read = stList_get(reads, i);
            switch (i) {
                case 0: // 410010     8M       ACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            (int64_t[]) {10, 11, 12, 13, 14, 15, 16, 17}, 8);
                    break;
                case 1: // 410010     2S8M2S   AAACGTACGTTT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            (int64_t[]) {10, 11, 12, 13, 14, 15, 16, 17}, 8);
                    break;
                case 2: // 410010     4S8M4S   AAAAACGTACGTTTTT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            (int64_t[]) {10, 11, 12, 13, 14, 15, 16, 17}, 8);
                    break;
                case 3: // 410012     8M       ACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            (int64_t[]) {12, 13, 14, 15, 16, 17, 18, 19}, 8);
                    break;
                case 4: // 410012     8M1I     ACGTACGTA
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            (int64_t[]) {12, 13, 14, 15, 16, 17, 18, 19}, 8);
                    break;
                case 5: // 410012     8M1D     ACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            (int64_t[]) {12, 13, 14, 15, 16, 17, 18, 19}, 8);
                    break;
                case 6: // 410012     7M2I     ACGTACGAA
                    CuAssertTrue(testCase, stString_eq("ACGTACGAA", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6}, (int64_t[]) {12, 13, 14, 15, 16, 17, 18},
                                            7);
                    break;
                case 7: // 410012     7M1D     ACGTACG
                    CuAssertTrue(testCase, stString_eq("ACGTACG", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6}, (int64_t[]) {12, 13, 14, 15, 16, 17, 18},
                                            7);
                    break;
                case 8: // 410012     2S8M2S   AAACGTACGTTT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            (int64_t[]) {12, 13, 14, 15, 16, 17, 18, 19}, 8);
                    break;
                case 9: // 410012     2S8M1I2S         AAACGTACGTATT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            (int64_t[]) {12, 13, 14, 15, 16, 17, 18, 19}, 8);
                    break;
                case 10: // 410012     2S8M1D2S         AAACGTACGTTT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            (int64_t[]) {12, 13, 14, 15, 16, 17, 18, 19}, 8);
                    break;
                case 11: // 410012    2S7M2I2S     AAACGTACGAATT
                    CuAssertTrue(testCase, stString_eq("ACGTACGAA", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6}, (int64_t[]) {12, 13, 14, 15, 16, 17, 18},
                                            7);
                    break;
                case 12:  // 410012    2S7M1D2S     AAACGTACGTT
                    CuAssertTrue(testCase, stString_eq("ACGTACG", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6}, (int64_t[]) {12, 13, 14, 15, 16, 17, 18},
                                            7);
                    break;
                case 13: // 410016     8M       ACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3}, (int64_t[]) {16, 17, 18, 19}, 4);
                    break;
                case 14: // 410016     3M1D4M   ACGACGT
                    CuAssertTrue(testCase, stString_eq("ACG", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2}, (int64_t[]) {16, 17, 18}, 3);
                    break;
                case 15: // 410016     3M2I4M   ACGCCTACG
                    CuAssertTrue(testCase, stString_eq("ACGCCT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 5}, (int64_t[]) {16, 17, 18, 19}, 4);
                    break;
                case 16: // 410016     2S8M2S   AAACGTACGTTT
                    CuAssertTrue(testCase, stString_eq("ACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3}, (int64_t[]) {16, 17, 18, 19}, 4);
                    break;
                case 17: // 410016     2S3M1D4M2S       AAACGACGTTT
                    CuAssertTrue(testCase, stString_eq("ACG", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2}, (int64_t[]) {16, 17, 18}, 3);
                    break;
                case 18: // 410016     2S3M2I4M2S       AAACGCCTACGTT
                    CuAssertTrue(testCase, stString_eq("ACGCCT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 5}, (int64_t[]) {16, 17, 18, 19}, 4);
                    break;
                case 19: // 410016     8M2S     ACGTACGTTT
                    CuAssertTrue(testCase, stString_eq("ACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3}, (int64_t[]) {16, 17, 18, 19}, 4);
                    break;
                case 20: // 410016     2S8M     AAACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3}, (int64_t[]) {16, 17, 18, 19}, 4);
                    break;
                case 21: // 410020     8M       ACGTACGT
                case 22: // 410020     2S8M2S   AAACGTACGTTT
                default:
                    CuAssertTrue(testCase, FALSE);
            }
        }
        stList_destruct(reads);
        stList_destruct(alignments);
    }
    CuAssertTrue(testCase, foundChunk);
    free(chunker->params);
    bamChunker_destruct(chunker);
}


static void test_readAlignmentsWithSoftclippingChunkEnd(CuTest *testCase) {
    BamChunker *chunker = bamChunker_construct(INPUT_BAM, getParameters(20, 0, TRUE));

    BamChunk *chunk = NULL;
    bool foundChunk = FALSE;
    for (int64_t chunkIdx = 0; chunkIdx < chunker->chunkCount; chunkIdx++) {
        chunk = bamChunker_getChunk(chunker, chunkIdx);
        if (!stString_eq(chunk->refSeqName, "contig_1") || !(chunk->chunkOverlapStart == 410000 &&
                                                             chunk->chunkOverlapEnd == 410020))
            continue;
        CuAssertTrue(testCase, !foundChunk);
        foundChunk = TRUE;

        // analyze reads and alignments
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *alignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
        uint32_t readCount = convertToReadsAndAlignments(chunk, NULL, reads, alignments, chunker->params);
        CuAssertTrue(testCase, readCount == 21);
        for (int64_t i = 0; i < 21; i++) {
            BamChunkRead *read = stList_get(reads, i);
            switch (i) {
                case 0: // 410010     8M       ACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            (int64_t[]) {10, 11, 12, 13, 14, 15, 16, 17}, 8);
                    break;
                case 1: // 410010     2S8M2S   AAACGTACGTTT
                    CuAssertTrue(testCase, stString_eq("AAACGTACGTTT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {2, 3, 4, 5, 6, 7, 8, 9},
                                            (int64_t[]) {10, 11, 12, 13, 14, 15, 16, 17}, 8);
                    break;
                case 2: // 410010     4S8M4S   AAAAACGTACGTTTTT
                    CuAssertTrue(testCase, stString_eq("AAAAACGTACGTTT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {4, 5, 6, 7, 8, 9, 10, 11},
                                            (int64_t[]) {10, 11, 12, 13, 14, 15, 16, 17}, 8);
                    break;
                case 3: // 410012     8M       ACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            (int64_t[]) {12, 13, 14, 15, 16, 17, 18, 19}, 8);
                    break;
                case 4: // 410012     8M1I     ACGTACGTA
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            (int64_t[]) {12, 13, 14, 15, 16, 17, 18, 19}, 8);
                    break;
                case 5: // 410012     8M1D     ACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6, 7},
                                            (int64_t[]) {12, 13, 14, 15, 16, 17, 18, 19}, 8);
                    break;
                case 6: // 410012     7M2I     ACGTACGAA
                    CuAssertTrue(testCase, stString_eq("ACGTACGAA", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6}, (int64_t[]) {12, 13, 14, 15, 16, 17, 18},
                                            7);
                    break;
                case 7: // 410012     7M1D     ACGTACG
                    CuAssertTrue(testCase, stString_eq("ACGTACG", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3, 4, 5, 6}, (int64_t[]) {12, 13, 14, 15, 16, 17, 18},
                                            7);
                    break;
                case 8: // 410012     2S8M2S   AAACGTACGTTT
                    CuAssertTrue(testCase, stString_eq("AAACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {2, 3, 4, 5, 6, 7, 8, 9},
                                            (int64_t[]) {12, 13, 14, 15, 16, 17, 18, 19}, 8);
                    break;
                case 9: // 410012     2S8M1I2S         AAACGTACGTATT
                    CuAssertTrue(testCase, stString_eq("AAACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {2, 3, 4, 5, 6, 7, 8, 9},
                                            (int64_t[]) {12, 13, 14, 15, 16, 17, 18, 19}, 8);
                    break;
                case 10: // 410012     2S8M1D2S         AAACGTACGTTT
                    CuAssertTrue(testCase, stString_eq("AAACGTACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {2, 3, 4, 5, 6, 7, 8, 9},
                                            (int64_t[]) {12, 13, 14, 15, 16, 17, 18, 19}, 8);
                    break;
                case 11: // 410012    2S7M2I2S     AAACGTACGAATT
                    CuAssertTrue(testCase, stString_eq("AAACGTACGAAT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {2, 3, 4, 5, 6, 7, 8, 9},
                                            (int64_t[]) {12, 13, 14, 15, 16, 17, 18}, 7);
                    break;
                case 12:  // 410012    2S7M1D2S     AAACGTACGTT
                    CuAssertTrue(testCase, stString_eq("AAACGTACG", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {2, 3, 4, 5, 6, 7, 8}, (int64_t[]) {12, 13, 14, 15, 16, 17, 18},
                                            7);
                    break;
                case 13: // 410016     8M       ACGTACGT
                    CuAssertTrue(testCase, stString_eq("ACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3}, (int64_t[]) {16, 17, 18, 19}, 4);
                    break;
                case 14: // 410016     3M1D4M   ACGACGT
                    CuAssertTrue(testCase, stString_eq("ACG", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2}, (int64_t[]) {16, 17, 18}, 3);
                    break;
                case 15: // 410016     3M2I4M   ACGCCTACG
                    CuAssertTrue(testCase, stString_eq("ACGCCT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 5}, (int64_t[]) {16, 17, 18, 19}, 4);
                    break;
                case 16: // 410016     2S8M2S   AAACGTACGTTT
                    CuAssertTrue(testCase, stString_eq("AAACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {2, 3, 4, 5}, (int64_t[]) {16, 17, 18, 19}, 4);
                    break;
                case 17: // 410016     2S3M1D4M2S       AAACGACGTTT
                    CuAssertTrue(testCase, stString_eq("AAACG", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {2, 3, 4}, (int64_t[]) {16, 17, 18}, 3);
                    break;
                case 18: // 410016     2S3M2I4M2S       AAACGCCTACGTT
                    CuAssertTrue(testCase, stString_eq("AAACGCCT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {2, 3, 4, 7}, (int64_t[]) {16, 17, 18, 19}, 4);
                    break;
                case 19: // 410016     8M2S     ACGTACGTTT
                    CuAssertTrue(testCase, stString_eq("ACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {0, 1, 2, 3}, (int64_t[]) {16, 17, 18, 19}, 4);
                    break;
                case 20: // 410016     2S8M     AAACGTACGT
                    CuAssertTrue(testCase, stString_eq("AAACGT", read->rleRead->rleString));
                    assertAlignmentMatching(testCase, stList_get(alignments, i),
                                            (int64_t[]) {2, 3, 4, 5}, (int64_t[]) {16, 17, 18, 19}, 4);
                    break;
                case 21: // 410020     8M       ACGTACGT
                case 22: // 410020     2S8M2S   AAACGTACGTTT
                default:
                    CuAssertTrue(testCase, FALSE);
            }
        }
        stList_destruct(alignments);
        stList_destruct(reads);
    }
    CuAssertTrue(testCase, foundChunk);
    free(chunker->params);
    bamChunker_destruct(chunker);
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
    Params *params = params_readParams(INPUT_PARAMS);
    char **chunks = st_calloc(4, sizeof(char *));
    chunks[0] = stString_copy("AAAAAAAACC");
    chunks[1] = stString_copy("AACCCCCCCCGG");
    chunks[2] = stString_copy("CCGGGGGGGGTT");
    chunks[3] = stString_copy("GGTTTTTTTT");
    ChunkToStitch **chunksToStitch = getChunksToStitchFromStrings(chunks, 4);
    ChunkToStitch *result = mergeContigChunkz(chunksToStitch, 0, 4, FALSE, params);
    CuAssertTrue(testCase, strcmp(result->seqHap1, "AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT") == 0);
}


void test_mergeContigChunksThreaded(CuTest *testCase) {
    Params *params = params_readParams(INPUT_PARAMS);
    params->polishParams->chunkBoundary = 2;
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


CuSuite *chunkingTestSuite(void) {
    CuSuite *suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_getRegionChunker);
    SUITE_ADD_TEST(suite, test_getChunksByChrom);
    SUITE_ADD_TEST(suite, test_getChunksBy100kb);
    SUITE_ADD_TEST(suite, test_getQualityScores);
    SUITE_ADD_TEST(suite, test_getChunksWithBoundary);
    SUITE_ADD_TEST(suite, test_getChunksWithoutBoundary);
    SUITE_ADD_TEST(suite, test_getReadsWithSoftClipping);
    SUITE_ADD_TEST(suite, test_getReadsWithoutSoftClipping);
    SUITE_ADD_TEST(suite, test_readAlignmentsWithoutSoftclippingChunkStart);
    SUITE_ADD_TEST(suite, test_readAlignmentsWithSoftclippingChunkStart);
    SUITE_ADD_TEST(suite, test_readAlignmentsWithoutSoftclippingChunkEnd);
    SUITE_ADD_TEST(suite, test_readAlignmentsWithSoftclippingChunkEnd);
    SUITE_ADD_TEST(suite, test_mergeContigChunks);
    SUITE_ADD_TEST(suite, test_mergeContigChunksThreaded);

    return suite;
}
