/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "margin.h"


static char* PARAMS = "../params/ont/r9.4/allParams.np.human.r94-g344.json";
static char* VCF1 = "../tests/data/vcfTest/vcfTest1.vcf";

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
    CuAssertTrue(testCase, rleString_eq(rleA1, entry->allele1));
    CuAssertTrue(testCase, rleString_eq(rleA2, entry->allele2));
    rleString_destruct(rleA1);
    rleString_destruct(rleA2);
}

void test_vcfParseRLE(CuTest *testCase) {
    bool RLE = TRUE;
    Params *params = params_readParams(PARAMS);
    params->polishParams->useRunLengthEncoding = RLE;

    stList *vcfEntries = parseVcf(VCF1, params);

    CuAssertTrue(testCase, stList_length(vcfEntries) == 7);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 0), "chr20", 1000, "G", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 1), "chr20", 2000, "T", "CCC", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 2), "chr20", 3000, "C", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 3), "chr20", 4000, "T", "C", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 4), "chr20", 5000, "GATTACA", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 5), "chr20", 6000, "T", "TC", RLE);
    //assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 6), "chr20", 7000, "G", "A", TRUE); //homozygous
    //assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 7), "chr20", 7000, "G", "A", TRUE); //homozygous
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 6), "chr20", 250000000, "A", "G", RLE);

    stList_destruct(vcfEntries);
    params_destruct(params);
}

void test_vcfParseRLEHOM(CuTest *testCase) {
    bool RLE = TRUE;
    Params *params = params_readParams(PARAMS);
    params->polishParams->useRunLengthEncoding = RLE;
    params->phaseParams->includeHomozygousVCFEntries = TRUE;

    stList *vcfEntries = parseVcf(VCF1, params);

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

    stList_destruct(vcfEntries);
    params_destruct(params);
}

void test_vcfParseRAW(CuTest *testCase) {
    bool RLE = FALSE;
    Params *params = params_readParams(PARAMS);
    params->polishParams->useRunLengthEncoding = RLE;

    stList *vcfEntries = parseVcf(VCF1, params);

    CuAssertTrue(testCase, stList_length(vcfEntries) == 7);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 0), "chr20", 1000, "G", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 1), "chr20", 2000, "T", "CCC", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 2), "chr20", 3000, "C", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 3), "chr20", 4000, "T", "C", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 4), "chr20", 5000, "GATTACA", "A", RLE);
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 5), "chr20", 6000, "T", "TC", RLE);
    //assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 6), "chr20", 7000, "G", "A", TRUE); //homozygous
    //assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 7), "chr20", 7000, "G", "A", TRUE); //homozygous
    assertVcfEntryCorrect(testCase, stList_get(vcfEntries, 6), "chr20", 250000000, "A", "G", RLE);

    stList_destruct(vcfEntries);
    params_destruct(params);
}

CuSuite *vcfTestSuite(void) {
    CuSuite *suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_vcfParseRLE);
    SUITE_ADD_TEST(suite, test_vcfParseRAW);
    SUITE_ADD_TEST(suite, test_vcfParseRLEHOM);

    return suite;
}
