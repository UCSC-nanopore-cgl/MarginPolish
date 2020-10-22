/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "margin.h"


static char* PARAMS = "../params/ont/r9.4/allParams.np.human.r94-g344.json";
static char* VCF1 = "../tests/data/vcfTest/vcfTest1.vcf";
static char* VCF1_GZ = "../tests/data/vcfTest/vcfTest1.vcf.gz";
static char* VCF2 = "../tests/data/vcfTest/vcfTest2.vcf";
static char* VCF2_REF = "../tests/data/vcfTest/vcfTest2.ref.fa";

/*
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr20	1000	.	G	A	30	.	.	GT:AD	0/1:1,1
chr20	2000	.	T	CCC	30	.	.	AD:GT	1,1:0|1
chr20	3000	.	T	C,A	30	.	.	GT:AD	1|2:1,1
chr20	4000	.	T	C	30	.	.	GT:AD:SE	0|1:1,1:SomethingElse
chr20	5000	.	A	GATTACA	30	.	.	GT:AD	1|0:1,1
chr20	6000	.	TC	T	30	.	.	GT:AD	1/0:1,1
chr20	7000	.	G	.	30	.	.	GT:AD	0|0:1,1
chr20	8000	.	G	A	30	.	.	GT:AD	1|1:1,1
chr20	250000000	.	G	A	30	.	.	GT:AD	1|0:1,1
*/


void assertVcfEntryCorrect(CuTest *testCase, VcfEntry *entry, char *ref, int64_t pos, char *allele1, char *allele2, bool rle) {
    RleString *rleA1 = rle ? rleString_construct(allele1) : rleString_construct_no_rle(allele1);
    RleString *rleA2 = rle ? rleString_construct(allele2) : rleString_construct_no_rle(allele2);
    CuAssertTrue(testCase, stString_eq(ref, entry->refSeqName));
    CuAssertTrue(testCase, pos == entry->refPos);
    CuAssertTrue(testCase, rleString_eq(rleA1, getVcfEntryAlleleH1(entry)));
    CuAssertTrue(testCase, rleString_eq(rleA2, getVcfEntryAlleleH2(entry)));
    rleString_destruct(rleA1);
    rleString_destruct(rleA2);
}

void test_vcfParseRLE(CuTest *testCase) {
    bool RLE = TRUE;
    Params *params = params_readParams(PARAMS);
    params->polishParams->useRunLengthEncoding = RLE;
    params->phaseParams->includeHomozygousVCFEntries = FALSE;

    stHash *vcfEntryMap = parseVcf(VCF1, params);
    stList *vcfEntries = stHash_search(vcfEntryMap, "chr20");

    CuAssertTrue(testCase, stList_length(vcfEntries) == 7);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 0), "chr20", 1000, "G", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 1), "chr20", 2000, "T", "CCC", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 2), "chr20", 3000, "C", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 3), "chr20", 4000, "T", "C", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 4), "chr20", 5000, "GATTACA", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 5), "chr20", 6000, "T", "TC", RLE);
    //assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 6), "chr20", 7000, "G", "G", RLE); //homozygous
    //assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 7), "chr20", 7000, "A", "A", RLE); //homozygous
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 6), "chr20", 250000000, "A", "G", RLE);

    stHash_destruct(vcfEntryMap);
    params_destruct(params);
}

void test_vcfParseRLEGZ(CuTest *testCase) {
    bool RLE = TRUE;
    Params *params = params_readParams(PARAMS);
    params->polishParams->useRunLengthEncoding = RLE;
    params->phaseParams->includeHomozygousVCFEntries = FALSE;

    stHash *vcfEntryMap = parseVcf(VCF1_GZ, params);
    stList *vcfEntries = stHash_search(vcfEntryMap, "chr20");

    CuAssertTrue(testCase, stList_length(vcfEntries) == 7);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 0), "chr20", 1000, "G", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 1), "chr20", 2000, "T", "CCC", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 2), "chr20", 3000, "C", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 3), "chr20", 4000, "T", "C", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 4), "chr20", 5000, "GATTACA", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 5), "chr20", 6000, "T", "TC", RLE);
    //assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 6), "chr20", 7000, "G", "G", RLE); //homozygous
    //assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 7), "chr20", 7000, "A", "A", RLE); //homozygous
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 6), "chr20", 250000000, "A", "G", RLE);

    stHash_destruct(vcfEntryMap);
    params_destruct(params);
}

void test_vcfParseRLESNP(CuTest *testCase) {
    bool RLE = TRUE;
    Params *params = params_readParams(PARAMS);
    params->polishParams->useRunLengthEncoding = RLE;
    params->phaseParams->includeHomozygousVCFEntries = FALSE;
    params->phaseParams->onlyUseSNPVCFEntries = TRUE;

    stHash *vcfEntryMap = parseVcf(VCF1, params);
    stList *vcfEntries = stHash_search(vcfEntryMap, "chr20");

    CuAssertTrue(testCase, stList_length(vcfEntries) == 4);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 0), "chr20", 1000, "G", "A", RLE);
    //assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 1), "chr20", 2000, "T", "CCC", RLE); // indel
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 1), "chr20", 3000, "C", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 2), "chr20", 4000, "T", "C", RLE);
    //assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 4), "chr20", 5000, "GATTACA", "A", RLE); // indel
    //assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 5), "chr20", 6000, "T", "TC", RLE); // indel
    //assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 6), "chr20", 7000, "G", "G", RLE); //homozygous
    //assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 7), "chr20", 7000, "A", "A", RLE); //homozygous
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 3), "chr20", 250000000, "A", "G", RLE);

    stHash_destruct(vcfEntryMap);
    params_destruct(params);
}

void test_vcfParseRLEHOM(CuTest *testCase) {
    bool RLE = TRUE;
    Params *params = params_readParams(PARAMS);
    params->polishParams->useRunLengthEncoding = RLE;
    params->phaseParams->includeHomozygousVCFEntries = TRUE;

    stHash *vcfEntryMap = parseVcf(VCF1, params);
    stList *vcfEntries = stHash_search(vcfEntryMap, "chr20");

    CuAssertTrue(testCase, stList_length(vcfEntries) == 9);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 0), "chr20", 1000, "G", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 1), "chr20", 2000, "T", "CCC", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 2), "chr20", 3000, "C", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 3), "chr20", 4000, "T", "C", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 4), "chr20", 5000, "GATTACA", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 5), "chr20", 6000, "T", "TC", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 6), "chr20", 7000, "G", "G", RLE); //homozygous
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 7), "chr20", 8000, "A", "A", RLE); //homozygous
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 8), "chr20", 250000000, "A", "G", RLE);

    stHash_destruct(vcfEntryMap);
    params_destruct(params);
}

void test_vcfParseRAW(CuTest *testCase) {
    bool RLE = FALSE;
    Params *params = params_readParams(PARAMS);
    params->polishParams->useRunLengthEncoding = RLE;
    params->phaseParams->includeHomozygousVCFEntries = TRUE;

    stHash *vcfEntryMap = parseVcf(VCF1, params);
    stList *vcfEntries = stHash_search(vcfEntryMap, "chr20");

    CuAssertTrue(testCase, stList_length(vcfEntries) == 9);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 0), "chr20", 1000, "G", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 1), "chr20", 2000, "T", "CCC", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 2), "chr20", 3000, "C", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 3), "chr20", 4000, "T", "C", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 4), "chr20", 5000, "GATTACA", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 5), "chr20", 6000, "T", "TC", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 6), "chr20", 7000, "G", "G", RLE); //homozygous
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 7), "chr20", 8000, "A", "A", RLE); //homozygous
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 8), "chr20", 250000000, "A", "G", RLE);

    stHash_destruct(vcfEntryMap);
    params_destruct(params);
}

void assertAlleleSubstringsCorrect2(CuTest *testCase, stList *alleleList, char **expectedAlleles, int expectedAlleleCount,
        stIntTuple *refPositions, int64_t expectedRefPosStart, int64_t expectedRefPosEnd) {
    CuAssertTrue(testCase, stList_length(alleleList) == expectedAlleleCount);
    if (refPositions != NULL) {
        int64_t refPosStart = stIntTuple_get(refPositions, 0);
        int64_t refPosEnd = stIntTuple_get(refPositions, 1);
        // these +1 are because these positions are in POA (1-based) space, not ref space (0-based)
        CuAssertTrue(testCase, refPosStart == expectedRefPosStart+1);
        CuAssertTrue(testCase, refPosEnd == expectedRefPosEnd+1);
    }
    for (int i = 0; i < stList_length(alleleList); i++) {
        char *allele = rleString_expand(stList_get(alleleList, i));
        char *expectedAllele = expectedAlleles[i];
        CuAssertTrue(testCase, stString_eq(allele, expectedAllele));
        free(allele);
    }
}

void assertAlleleSubstringsCorrect(CuTest *testCase, stList *alleleList, char **expectedAlleles, int expectedAlleleCount) {
    assertAlleleSubstringsCorrect2(testCase, alleleList, expectedAlleles, expectedAlleleCount, NULL, 0, 0);
}


void test_vcfAlleleSubstrings(CuTest *testCase) {
    bool RLE = FALSE;
    Params *params = params_readParams(PARAMS);
    params->polishParams->useRunLengthEncoding = RLE;
    params->phaseParams->includeHomozygousVCFEntries = TRUE;
    params->phaseParams->onlyUseSNPVCFEntries = FALSE;
    params->polishParams->columnAnchorTrim = 2;

    stHash *vcfEntryMap = parseVcf(VCF2, params);
    stList *vcfEntries = stHash_search(vcfEntryMap, "vcfTest2");

    // just checking that we read the VCF correct
    CuAssertTrue(testCase, stList_length(vcfEntries) == 11);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 0), "vcfTest2", 0, "A", "G", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 1), "vcfTest2", 1, "A", "G", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 2), "vcfTest2", 32, "A", "G", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 3), "vcfTest2", 48, "A", "C", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 4), "vcfTest2", 64, "G", "GCC", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 5), "vcfTest2", 72, "A", "ACTG", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 6), "vcfTest2", 80, "GGG", "G", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 7), "vcfTest2", 88, "AGGG", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 8), "vcfTest2", 96, "CCC", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 9), "vcfTest2", 126, "A", "G", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 10), "vcfTest2", 127, "C", "A", RLE);

    // get substrings
    // this conversion is to put things in poa-space, which should be countered by allele substrings
    vcfEntries = getVcfEntriesForRegion(vcfEntryMap, NULL, "vcfTest2", 0, 128);
    stHash *referenceSequences = parseReferenceSequences(VCF2_REF);
    char *refSeq = stHash_search(referenceSequences, "vcfTest2");
    RleString *refRleString = RLE ? rleString_construct(refSeq) : rleString_construct_no_rle(refSeq);
    stList *allAlleleSubstrings = stList_construct3(0, (void (*)(void*)) stList_destruct);
    stList *allRefPositions = stList_construct3(0, (void (*)(void*)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(vcfEntries); i++) {
        int64_t start, end;
        stList_append(allAlleleSubstrings, getAlleleSubstrings(stList_get(vcfEntries, i), refRleString,
                                                               params, &start, &end));
        stList_append(allRefPositions, stIntTuple_construct2(start, end));
    }

    // assert correctness of strings
    char *alleles0[] = {"AAA", "GAA"};
    assertAlleleSubstringsCorrect2(testCase, stList_get(allAlleleSubstrings, 0), alleles0, 2,
            stList_get(allRefPositions, 0), 0, 3);
    char *alleles1[] = {"AAAA", "AGAA"};
    assertAlleleSubstringsCorrect2(testCase, stList_get(allAlleleSubstrings, 1), alleles1, 2,
                                  stList_get(allRefPositions, 1), 0, 4);
    char *alleles2[] = {"TTAGA", "TTGGA"};
    assertAlleleSubstringsCorrect2(testCase, stList_get(allAlleleSubstrings, 2), alleles2, 2,
                                   stList_get(allRefPositions, 2), 30, 35);
    char *alleles3[] = {"CGAAC", "CGCAC", "CGGAC", "CGTAC"};
    assertAlleleSubstringsCorrect2(testCase, stList_get(allAlleleSubstrings, 3), alleles3, 4,
                                   stList_get(allRefPositions, 3), 46, 51);
    char *alleles4[] = {"ATGAC", "ATGCCAC"};
    assertAlleleSubstringsCorrect2(testCase, stList_get(allAlleleSubstrings, 4), alleles4, 2,
                                   stList_get(allRefPositions, 4), 62, 67);
    char *alleles5[] = {"CCAGA", "CCACTGGA", "CCCCCGA"};
    assertAlleleSubstringsCorrect2(testCase, stList_get(allAlleleSubstrings, 5), alleles5, 3,
                                   stList_get(allRefPositions, 5), 70, 75);
    char *alleles6[] = {"ACGGGAG", "ACGAG"};
    assertAlleleSubstringsCorrect2(testCase, stList_get(allAlleleSubstrings, 6), alleles6, 2,
                                  stList_get(allRefPositions, 6), 78, 85);
    char *alleles7[] = {"CCAGGGGA", "CCAGA", "CCAGGA"};
    assertAlleleSubstringsCorrect2(testCase, stList_get(allAlleleSubstrings, 7), alleles7, 3,
                                   stList_get(allRefPositions, 7), 86, 94);
    char *alleles8[] = {"CACCCAA", "CAAAA", "CAGGAAA", "CACAGAGAGAAA"};
    assertAlleleSubstringsCorrect2(testCase, stList_get(allAlleleSubstrings, 8), alleles8, 4,
                                   stList_get(allRefPositions, 8), 94, 101);
    char *alleles9[] = {"ATAC", "ATGC"};
    assertAlleleSubstringsCorrect2(testCase, stList_get(allAlleleSubstrings, 9), alleles9, 2,
                                   stList_get(allRefPositions, 9), 124, 128);
    char *alleles10[] = {"TAC", "TAA"};
    assertAlleleSubstringsCorrect2(testCase, stList_get(allAlleleSubstrings, 10), alleles10, 2,
                                   stList_get(allRefPositions, 10), 125, 128);
    rleString_destruct(refRleString);

    // same for alleles in specific region
    stList *regionVcfEntries = getVcfEntriesForRegion(vcfEntryMap, NULL, "vcfTest2", 64, 128);
    char *regionRefSubstring = stString_getSubString(refSeq, 64, 64);
    refRleString = RLE ? rleString_construct(regionRefSubstring) : rleString_construct_no_rle(regionRefSubstring);
    stList *allRegionAlleleSubstrings = stList_construct3(0, (void (*)(void*)) stList_destruct);
    stList *allRegionRefPositions = stList_construct3(0, (void (*)(void*)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(regionVcfEntries); i++) {
        int64_t start, end;
        stList_append(allRegionAlleleSubstrings, getAlleleSubstrings(stList_get(regionVcfEntries, i),
                                                                     refRleString, params, &start, &end));
        stList_append(allRegionRefPositions, stIntTuple_construct2(start, end));
    }

    CuAssertTrue(testCase, stList_length(allRegionAlleleSubstrings) == 7);
    char *regionAlleles4[] = {"GAC", "GCCAC"};
    assertAlleleSubstringsCorrect2(testCase, stList_get(allRegionAlleleSubstrings, 0), regionAlleles4, 2,
                                   stList_get(allRegionRefPositions, 0), 0, 3);
    assertAlleleSubstringsCorrect2(testCase, stList_get(allRegionAlleleSubstrings, 1), alleles5, 3,
                                   stList_get(allRegionRefPositions, 1), 6, 11);
    assertAlleleSubstringsCorrect2(testCase, stList_get(allRegionAlleleSubstrings, 2), alleles6, 2,
                                   stList_get(allRegionRefPositions, 2), 14, 21);
    assertAlleleSubstringsCorrect2(testCase, stList_get(allRegionAlleleSubstrings, 3), alleles7, 3,
                                   stList_get(allRegionRefPositions, 3), 22, 30);
    assertAlleleSubstringsCorrect2(testCase, stList_get(allRegionAlleleSubstrings, 4), alleles8, 4,
                                   stList_get(allRegionRefPositions, 4), 30, 37);
    assertAlleleSubstringsCorrect2(testCase, stList_get(allRegionAlleleSubstrings, 5), alleles9, 2,
                                   stList_get(allRegionRefPositions, 5), 60, 64);
    assertAlleleSubstringsCorrect2(testCase, stList_get(allRegionAlleleSubstrings, 6), alleles10, 2,
                                   stList_get(allRegionRefPositions, 6), 61, 64);

    rleString_destruct(refRleString);
    stList_destruct(allRegionRefPositions);
    stList_destruct(allRegionAlleleSubstrings);
    stList_destruct(allRefPositions);
    stList_destruct(allAlleleSubstrings);
    stList_destruct(regionVcfEntries);
    stHash_destruct(vcfEntryMap);
    params_destruct(params);
}

void test_vcfBinarySearchForVcfEntryStartIdx(CuTest *testCase) {
    for (int64_t testIdx = 0; testIdx < 100; testIdx++) {
        stList *vcfEntries = stList_construct3(0, (void(*)(void*)) vcfEntry_destruct);
        int64_t refPos = 0;
        for (int64_t entryIdx = 0; entryIdx < st_randomInt64(32, 512); entryIdx++) {
            refPos += st_randomInt64(0, 16);
            stList_append(vcfEntries, vcfEntry_construct("", refPos, refPos, -1, NULL, -1, -1));
        }

        // ensure we test pos 0, posOutsideRange, and the rest
        int64_t desiredRefPos = testIdx == 0 ? 0 : testIdx == 1 ? refPos + 1 : st_randomInt64(0, refPos);

        int64_t correctIdx = -1;
        for (int64_t entryIdx = 0; entryIdx < stList_length(vcfEntries); entryIdx++) {
            VcfEntry *entry = stList_get(vcfEntries, entryIdx);
            if (entry->refPos >= desiredRefPos) {
                correctIdx = entryIdx;
                break;
            }
        }

        int64_t queryIdx = binarySearchVcfListForFirstIndexAfterRefPos(vcfEntries, desiredRefPos);

        if (queryIdx != correctIdx) {
            CuAssertTrue(testCase, queryIdx == correctIdx);
        }

        stList_destruct(vcfEntries);

    }
}

CuSuite *vcfTestSuite(void) {
    CuSuite *suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_vcfParseRLE);
    SUITE_ADD_TEST(suite, test_vcfParseRLEGZ);
    SUITE_ADD_TEST(suite, test_vcfParseRAW);
    SUITE_ADD_TEST(suite, test_vcfParseRLEHOM);
    SUITE_ADD_TEST(suite, test_vcfParseRLESNP);
    SUITE_ADD_TEST(suite, test_vcfAlleleSubstrings);
    SUITE_ADD_TEST(suite, test_vcfBinarySearchForVcfEntryStartIdx);

    return suite;
}
