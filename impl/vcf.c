//
// Created by tpesout on 5/28/20.
//

#include "margin.h"

VcfEntry *vcfEntry_construct(char *refSeqName, int64_t refPos, int64_t rawRefPos, double phredQuality,
        RleString *allele1, RleString *allele2) {
    VcfEntry *vcfEntry = st_calloc(1, sizeof(VcfEntry));
    vcfEntry->refSeqName = stString_copy(refSeqName);
    vcfEntry->refPos = refPos;
    vcfEntry->rawRefPosInformativeOnly = rawRefPos;
    vcfEntry->phredQuality = phredQuality;
    vcfEntry->allele1 = allele1;
    vcfEntry->allele2 = allele2;
    return vcfEntry;
}

void vcfEntry_destruct(VcfEntry *vcfEntry) {
    if (vcfEntry->allele1 != NULL) rleString_destruct(vcfEntry->allele1);
    if (vcfEntry->allele2 != NULL) rleString_destruct(vcfEntry->allele2);
    free(vcfEntry->refSeqName);
    free(vcfEntry);
}

#define GQ_WEIGHT .8
#define  Q_WEIGHT .2
#define DEFAULT_MIN_VCF_QUAL -1

stList *parseVcf(char *vcfFile, Params *params) {
    return parseVcf2(vcfFile, NULL, params);
}
stList *parseVcf2(char *vcfFile, char *regionStr, Params *params) {
    stList *entries = stList_construct3(0, (void(*)(void*))vcfEntry_destruct);
    FILE *fp = safe_fopen(vcfFile, "r");
    if (fp == NULL) {
        st_errAbort("Could not open VCF %s\n", vcfFile);
    }
    char regionContig[128] = "";
    int regionStart = 0;
    int regionEnd = 0;
    if (regionStr != NULL) {
        int scanRet = sscanf(regionStr, "%[^:]:%d-%d", regionContig, &regionStart, &regionEnd);
        if (scanRet != 3 && scanRet != 1) {
            st_errAbort("Region in unexpected format (expected %%s:%%d-%%d or %%s)): %s", regionStr);
        } else if (regionStart < 0 || regionEnd < 0 || regionEnd < regionStart) {
            st_errAbort("Start and end locations in region must be positive, start must be less than end: %s", regionStr);
        }
        if (scanRet == 1) {
            regionStart = -1;
            regionEnd = -1;
        }
    }
    int64_t skippedForRegion = 0;
    int64_t skippedForNotPass = 0;
    int64_t skippedForHomozygous = 0;

    char *line = NULL;
    while ((line = stFile_getLineFromFile(fp)) != NULL) {
        if (strlen(line) == 0 || line[0] == '#') {
            free(line);
            continue;
        }
        stList *elements = stString_split(line);
        if (stList_length(elements) < 10) {
            st_errAbort("Malformed VCF line: %s\n", line);
        }

        // location data
        char *chrom = stList_get(elements, 0);
        int64_t pos = atol(stList_get(elements, 1));

        // quick fail
        if (regionStr != NULL && (!stString_eq(regionContig, chrom) || (regionStart >= 0 && !(regionStart <= pos && pos < regionEnd)))) {
            skippedForRegion++;
            free(line);
            stList_destruct(elements);
            continue;
        }
        if (params->phaseParams->onlyUsePassVCFEntries && !(
                stString_eqcase("PASS", stList_get(elements, 6)) || stString_eq(".", stList_get(elements, 6)))) {
            skippedForNotPass++;
            free(line);
            stList_destruct(elements);
            continue;
        }

        // get genotype idx
        stList *format = stString_splitByString(stList_get(elements, 8), ":");
        int64_t gtIdx = -1;
        int64_t gqIdx = -1;
        for (int64_t i = 0; i < stList_length(format); i++) {
            if (stString_eq("GT", stList_get(format, i))) {
                gtIdx = i;
            } else if (stString_eq("GQ", stList_get(format, i))) {
                gqIdx = i;
            }
        }
        if (gtIdx == -1) {
            stList_destruct(format);
            free(line);
            st_logInfo("  Missing GT record for format string '%s' in VCF line:\n\t\t%s\n", stList_get(elements, 8), line);
            continue;
        }

        // get genotype
        stList *sample = stString_splitByString(stList_get(elements, 9), ":");
        char *genotypeStr = stList_get(sample, gtIdx);
        if (strlen(genotypeStr) != 3) {
            stList_destruct(sample);
            stList_destruct(format);
            stList_destruct(elements);
            free(line);
            st_logInfo("  Unexpected genotype str '%s' in VCF line:\n\t\t%s\n", genotypeStr, line);
            continue;
        }

        // init gt to ref
        int64_t gt1 = 0;
        int64_t gt2 = 0;

        // fail
        if (genotypeStr[0] != '.') {
            gt1 = genotypeStr[0] - '0';
        }
        if (genotypeStr[2] != '.') {
            gt2 = genotypeStr[2] - '0';
        }
        assert(gt1 >= 0 && gt2 >=0);

        // get variant quality
        double variantQuality = -1;
        double genotypeQuality = -1;
        double quality = DEFAULT_MIN_VCF_QUAL;
        if (gqIdx != -1) {
            genotypeQuality = atof(stList_get(sample,gqIdx));
        }
        if (!stString_eq(".", stList_get(elements, 5))) {
            quality = atof(stList_get(elements, 5));
        }
        if (genotypeQuality > 0 && quality > 0) {
            variantQuality = toPhred( GQ_WEIGHT * fromPhred(genotypeQuality) + Q_WEIGHT * fromPhred(quality));
        } else if (genotypeQuality > 0) {
            variantQuality = genotypeQuality;
        } else if (quality > 0) {
            variantQuality = quality;
        }

        // get alleles
        stList *alleles = stList_construct(); //ref is freed in elements, alts are freed in altAlleles
        stList_append(alleles, stList_get(elements, 3)); //ref
        stList *altAlleles = stString_splitByString(stList_get(elements, 4), ",");
        stList_appendAll(alleles, altAlleles);
        assert(stList_length(alleles) >= (gt1 > gt2 ? gt1 : gt2));
        char *allele1 = stString_copy(stList_get(alleles, gt1));
        char *allele2 = stString_copy(stList_get(alleles, gt2));
        if (stString_eq(allele2, ".")) {
            free(allele2);
            allele2 = stString_copy(stList_get(alleles, 0)); //ref
        }

        // save it
        if (gt1 != gt2 || params->phaseParams->includeHomozygousVCFEntries) {
            RleString *rleAllele1 = params->polishParams->useRunLengthEncoding ? rleString_construct(allele1)
                    : rleString_construct_no_rle(allele1);
            RleString *rleAllele2 = params->polishParams->useRunLengthEncoding ? rleString_construct(allele2)
                    : rleString_construct_no_rle(allele2);
            VcfEntry *entry = vcfEntry_construct(chrom, pos, pos, variantQuality, rleAllele1, rleAllele2);
            stList_append(entries, entry);
        } else {
            skippedForHomozygous++;
        }

        // cleanup
        stList_destruct(alleles);
        stList_destruct(altAlleles);
        stList_destruct(format);
        stList_destruct(sample);
        stList_destruct(elements);
        free(line);
        free(allele1);
        free(allele2);

    }

    // cleanup
    fclose(fp);

    st_logCritical("> Parsed %"PRId64" %sVCF entries from %s; skipped %"PRId64" for region, %"PRId64" for not being "
                   "PASS, %"PRId64" for being homozygous\n",
            stList_length(entries), params->phaseParams->includeHomozygousVCFEntries ? " " : "HET ", vcfFile,
            skippedForRegion, skippedForNotPass, skippedForHomozygous);
    if (stList_length(entries) == 0) {
        st_errAbort("No valid VCF entries found!");
    }
    return entries;
}

stList *getVcfEntriesForRegion2(stList *vcfEntries, uint64_t *rleMap, char *refSeqName, int64_t startPos, int64_t endPos, double minQual) {
    stList *regionEntries = stList_construct3(0, (void(*)(void*))vcfEntry_destruct);
    int64_t qualityFilteredCount = 0;
    for (int64_t i = 0; i < stList_length(vcfEntries); i++) {
        VcfEntry *e = stList_get(vcfEntries, i);
        if (!stString_eq(refSeqName, e->refSeqName)) continue;
        if (startPos > e->refPos) continue;
        if (endPos <= e->refPos) continue;
        if (minQual > e->phredQuality) {
            qualityFilteredCount++;
            continue;
        }
        int64_t refPos = e->refPos - startPos;
        if (rleMap != NULL) {
            refPos = rleMap[refPos];
        }
        VcfEntry *copy = vcfEntry_construct(e->refSeqName, refPos, e->rawRefPosInformativeOnly, e->phredQuality,
                rleString_copy(e->allele1), rleString_copy(e->allele2));
        stList_append(regionEntries, copy);
    }
    if (minQual != DEFAULT_MIN_VCF_QUAL) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Filtered %"PRIu64" VCF records for quality below %.2f, with %"PRIu64" remaining.\n",
                   logIdentifier, minQual, stList_length(regionEntries));
        free(logIdentifier);
    }
    return regionEntries;
}


stList *getVcfEntriesForRegion(stList *vcfEntries, uint64_t *rleMap, char *refSeqName, int64_t startPos, int64_t endPos) {
    return getVcfEntriesForRegion2(vcfEntries, rleMap, refSeqName, startPos, endPos, DEFAULT_MIN_VCF_QUAL);
}