//
// Created by tpesout on 5/28/20.
//

#include "margin.h"
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>


VcfEntry *vcfEntry_construct(char *refSeqName, int64_t refPos, int64_t rawRefPos, double phredQuality,
        stList *alleles, int64_t gt1, int64_t gt2) {
    VcfEntry *vcfEntry = st_calloc(1, sizeof(VcfEntry));
    vcfEntry->refSeqName = stString_copy(refSeqName);
    vcfEntry->refPos = refPos;
    vcfEntry->rawRefPosInformativeOnly = rawRefPos;
    vcfEntry->phredQuality = phredQuality;
    vcfEntry->alleles = alleles == NULL ? stList_construct3(0, (void (*)(void*))rleString_destruct) : alleles;
    vcfEntry->gt1 = gt1;
    vcfEntry->gt2 = gt2;
    vcfEntry->alleleSubstrings = NULL;
    vcfEntry->refAlnStart = -1;
    vcfEntry->refAlnStopExcl = -1;
    return vcfEntry;
}

void vcfEntry_destruct(VcfEntry *vcfEntry) {
    stList_destruct(vcfEntry->alleles);
    if (vcfEntry->alleleSubstrings != NULL) stList_destruct(vcfEntry->alleleSubstrings);
    free(vcfEntry->refSeqName);
    free(vcfEntry);
}

RleString *getVcfEntryAlleleH1(VcfEntry *vcfEntry) {
    assert(vcfEntry->gt1 < stList_length(vcfEntry->alleles));
    return stList_get(vcfEntry->alleles, vcfEntry->gt1);
}
RleString *getVcfEntryAlleleH2(VcfEntry *vcfEntry) {
    assert(vcfEntry->gt2 < stList_length(vcfEntry->alleles));
    return stList_get(vcfEntry->alleles, vcfEntry->gt2);
}

#define GQ_WEIGHT .8
#define  Q_WEIGHT .2
#define DEFAULT_MIN_VCF_QUAL -1

stHash *parseVcf(char *vcfFile, Params *params) {
    return parseVcf2(vcfFile, NULL, params);
}


stHash *getContigToVariantMap(stList *vcfEntries) {
    stHash *map = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, (void(*)(void*))stList_destruct);
    for (int64_t i = 0; i < stList_length(vcfEntries); i++) {
        VcfEntry *vcfEntry = stList_get(vcfEntries, i);
        stList *contigList = stHash_search(map, vcfEntry->refSeqName);
        if (contigList == NULL) {
            contigList = stList_construct3(0, (void(*)(void*))vcfEntry_destruct);
            stHash_insert(map, stString_copy(vcfEntry->refSeqName), contigList);
        }
        stList_append(contigList, vcfEntry);
    }
    return map;
}

int vcfEntryCmpFn(const void *a, const void *b) {
    VcfEntry *A = (VcfEntry*) a;
    VcfEntry *B = (VcfEntry*) b;
    if (A->refPos == B->refPos) return 0;
    return A->refPos < B->refPos ? -1 : 1;
}

stHash *parseVcf2(char *vcfFile, char *regionStr, Params *params) {
    stHash *entries = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, (void(*)(void*))stList_destruct);

    //open vcf file
    htsFile *fp = hts_open(vcfFile,"rb");
    if (fp == NULL) {
        st_errAbort("Could not open VCF %s\n", vcfFile);
    }

    // region manage
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
    int64_t totalEntries = 0;
    int64_t skippedForRegion = 0;
    int64_t skippedForIndel = 0;
    int64_t skippedForNotPass = 0;
    int64_t skippedForHomozygous = 0;

    //read header
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    int nsmpl = bcf_hdr_nsamples(hdr);
    if (nsmpl > 1) {
        st_logCritical("> Got %d samples reading %s, will only take VCF records for the first\n", nsmpl, vcfFile);
    }


    bcf1_t *rec    = bcf_init();

    //save for each vcf record
    while ( bcf_read(fp, hdr, rec) >= 0 )
    {
        //unpack for read REF,ALT,INFO,etc
        bcf_unpack(rec, BCF_UN_ALL);
        totalEntries++;

        // location data
        char *chrom = bcf_hdr_id2name(hdr, rec->rid);
        int64_t pos = rec->pos;

        // quick fail
        if (regionStr != NULL && (!stString_eq(regionContig, chrom) || (regionStart >= 0 && !(regionStart <= pos && pos < regionEnd)))) {
            skippedForRegion++;
            continue;
        }
        if (params->phaseParams->onlyUsePassVCFEntries && !bcf_has_filter(hdr, rec, "PASS")) {
            skippedForNotPass++;
            continue;
        }
        if (params->phaseParams->onlyUseSNPVCFEntries && !bcf_is_snp(rec)) {
            skippedForIndel++;
            continue;
        }


        // genotype
        int32_t *gt_arr = NULL, ngt_arr = 0;
        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt<=0 || bcf_gt_is_missing(gt_arr[0]) || gt_arr[1] == bcf_int32_vector_end) {
            //todo
            return NULL; // GT not present
        }
        int gt1 = bcf_gt_allele(gt_arr[0]);
        int gt2 = bcf_gt_allele(gt_arr[1]);
        free(gt_arr);
        if (!params->phaseParams->includeHomozygousVCFEntries && gt1 == gt2) {
            skippedForHomozygous++;
            continue;
        }


        // todo variant quality
        double variantQuality = -1;
        double genotypeQuality = -1;
        double quality = DEFAULT_MIN_VCF_QUAL;
        if (genotypeQuality > 0 && quality > 0) {
            variantQuality = toPhred( GQ_WEIGHT * fromPhred(genotypeQuality) + Q_WEIGHT * fromPhred(quality));
        } else if (genotypeQuality > 0) {
            variantQuality = genotypeQuality;
        } else if (quality > 0) {
            variantQuality = quality;
        }

        // get alleles
        stList *alleles = stList_construct3(0, (void (*)(void*)) rleString_destruct);
        for (int i=0; i<rec->n_allele; ++i) {
            stList_append(alleles,
                    params->polishParams->useRunLengthEncoding ?
                    rleString_construct(rec->d.allele[i]) :
                    rleString_construct_no_rle(rec->d.allele[i]));
        }

        // save it
        VcfEntry *entry = vcfEntry_construct(chrom, pos, pos, variantQuality, alleles, gt1, gt2);
        stList *contigList = stHash_search(entries, entry->refSeqName);
        if (contigList == NULL) {
            contigList = stList_construct3(0, (void(*)(void*))vcfEntry_destruct);
            stHash_insert(entries, stString_copy(entry->refSeqName), contigList);
        }
        stList_append(contigList, entry);
    }

    // cleanup
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ( (ret=hts_close(fp)) ) {
        st_logCritical("> Failed to close VCF %s with code %d\n", vcfFile, ret);
    }

    // logging
    st_logCritical("> Parsed %"PRId64" %sVCF entries from %s; skipped %"PRId64" for region, %"PRId64" for not being "
                   "PASS, %"PRId64" for being homozygous, %"PRId64" for being INDEL\n",
                   totalEntries, params->phaseParams->includeHomozygousVCFEntries ? " " : "HET ", vcfFile,
                   skippedForRegion, skippedForNotPass, skippedForHomozygous, skippedForIndel);
    if (totalEntries == 0) {
        st_errAbort("No valid VCF entries found!");
    }

    // ensure sorted
    stHashIterator *itor = stHash_getIterator(entries);
    char *contigName = NULL;
    while ((contigName = stHash_getNext(itor)) != NULL) {
        stList *contigEntries = stHash_search(entries, contigName);
        stList_sort(contigEntries, vcfEntryCmpFn);
        assert(((VcfEntry*) stList_get(contigEntries, 0))->refPos <= ((VcfEntry*) stList_get(contigEntries,
                stList_length(contigEntries) - 1))->refPos);
    }
    stHash_destructIterator(itor);

    // finish
    return entries;
}

stList *copyListOfRleStrings(stList *toCopy) {
    stList *copy = stList_construct3(0, (void (*)(void*)) rleString_destruct);
    for (int i = 0; i < stList_length(toCopy); i++) {
        stList_append(copy, rleString_copy(stList_get(toCopy, i)));
    }
    return copy;
}

int64_t binarySearchVcfListForFirstIndexAfterRefPos2(stList *vcfEntries, int64_t desiredEntryPos, int64_t startPos, int64_t endPosIncl) {
    if (endPosIncl - startPos == 1) {
        return ((VcfEntry*)stList_get(vcfEntries, startPos))->refPos >= desiredEntryPos ? startPos : endPosIncl;
    }
    int64_t middlePos = startPos + (endPosIncl - startPos) / 2;
    VcfEntry *middleEntry = stList_get(vcfEntries, middlePos);
    if (middleEntry->refPos < desiredEntryPos) {
        return binarySearchVcfListForFirstIndexAfterRefPos2(vcfEntries, desiredEntryPos, middlePos, endPosIncl);
    } else /*if (middleEntry->refPos >= desiredEntryPos)*/ {
        return binarySearchVcfListForFirstIndexAfterRefPos2(vcfEntries, desiredEntryPos, startPos, middlePos);
    }

}
int64_t binarySearchVcfListForFirstIndexAfterRefPos(stList *vcfEntries, int64_t refPos) {
    if (stList_length(vcfEntries) == 0) return -1;
    if (((VcfEntry*)stList_get(vcfEntries, stList_length(vcfEntries) - 1))->refPos < refPos) return -1;
    if (((VcfEntry*)stList_get(vcfEntries, 0))->refPos > refPos) return 0;
    return binarySearchVcfListForFirstIndexAfterRefPos2(vcfEntries, refPos, 0, stList_length(vcfEntries) - 1);
}

stList *getVcfEntriesForRegion2(stHash *vcfEntryMap, uint64_t *rleMap, char *refSeqName, int64_t startPos, int64_t endPos, double minQual) {
    // get entries and sanity check
    stList *vcfEntries = stHash_search(vcfEntryMap, refSeqName);
    if (vcfEntries == NULL) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Reference Sequence %s not found in VCF entries\n", logIdentifier, refSeqName);
        free(logIdentifier);
        return NULL;
    }

    // the entries for this region we want
    stList *regionEntries = stList_construct3(0, (void(*)(void*))vcfEntry_destruct);

    // binary search through list for start pos
    int64_t startIdx = binarySearchVcfListForFirstIndexAfterRefPos(vcfEntries, startPos);
    if (startIdx == -1) {
        return regionEntries;
    }

    // get all entries from start until we're out
    int64_t qualityFilteredCount = 0;
    for (int64_t i = startIdx; i < stList_length(vcfEntries); i++) {
        VcfEntry *e = stList_get(vcfEntries, i);
        assert(stString_eq(refSeqName, e->refSeqName));
        assert(startPos <= e->refPos);
        if (endPos <= e->refPos) break;
        if (minQual > e->phredQuality) {
            qualityFilteredCount++;
            continue;
        }
        int64_t refPos = e->refPos - startPos + 1; // e->refPos is in 0-based space, poa is 1-based. convert here
        if (rleMap != NULL) {
            refPos = rleMap[refPos];
        }
        VcfEntry *copy = vcfEntry_construct(e->refSeqName, refPos, e->rawRefPosInformativeOnly, e->phredQuality,
                copyListOfRleStrings(e->alleles), e->gt1, e->gt2);
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


stList *getVcfEntriesForRegion(stHash *vcfEntryMap, uint64_t *rleMap, char *refSeqName, int64_t startPos, int64_t endPos) {
    return getVcfEntriesForRegion2(vcfEntryMap, rleMap, refSeqName, startPos, endPos, DEFAULT_MIN_VCF_QUAL);
}


stList *getAlleleSubstrings2(VcfEntry *entry, char *referenceSeq, int64_t refSeqLen, int64_t *refStartPos,
        int64_t *refEndPosExcl, int64_t expansion, bool useRunLengthEncoding) {
    stList *substrings = stList_construct3(0, (void (*)(void*)) rleString_destruct);
    assert(stList_length(entry->alleles) >= 1);

    // parameters for substringing
    int64_t pos = entry->refPos;
    pos--; // at this point, reference positions are in poa-space (1-based), need to convert to ref-space (0-based)

    //get ref info
    char *refAllele = rleString_expand(stList_get(entry->alleles, 0));
    int64_t refAlleleLen = strlen(refAllele);
    for (int64_t i = 0; i < refAlleleLen; i++) {
        assert(referenceSeq[pos + i] == refAllele[i]);
    }
    free(refAllele);

    // get prefix and suffix strings (from ref seq)
    int64_t pStart = pos - expansion;
    int64_t sStart = pos + refAlleleLen;
    *refStartPos = pStart < 0 ? 0 : pStart;
    *refEndPosExcl = sStart + expansion >= refSeqLen ? refSeqLen : sStart + expansion;
    char *prefix = stString_getSubString(referenceSeq, *refStartPos, pStart < 0 ? pos : expansion);
    char *suffix = stString_getSubString(referenceSeq, sStart, sStart + expansion >= refSeqLen ? refSeqLen - sStart : expansion);

    // get alleles
    for (int64_t i = 0; i < stList_length(entry->alleles); i++) {
        char *expandedAllele = rleString_expand(stList_get(entry->alleles, i));
        char *fullAlleleSubstring = stString_print("%s%s%s", prefix, expandedAllele, suffix);
        stList_append(substrings,
                useRunLengthEncoding ?
                rleString_construct(fullAlleleSubstring) :
                rleString_construct_no_rle(fullAlleleSubstring));
        free(expandedAllele);
        free(fullAlleleSubstring);
    }

    // cleanup
    free(prefix);
    free(suffix);

    // put refStartPos and endPos back in poa-space
    (*refStartPos)++;
    (*refEndPosExcl)++;
    return substrings;
}

stList *getAlleleSubstrings(VcfEntry *entry, RleString *referenceSeq, Params *params,
                            int64_t *refStartPos, int64_t *refEndPosExcl) {
    char *rawRefSeq = rleString_expand(referenceSeq);
    stList *alleleSubstrings = getAlleleSubstrings2(entry, rawRefSeq, referenceSeq->nonRleLength, refStartPos, refEndPosExcl,
                                                    params->polishParams->columnAnchorTrim,
                                                    params->polishParams->useRunLengthEncoding);
    free(rawRefSeq);
    return alleleSubstrings;

}

void updateVcfEntriesWithSubstringsAndPositions(stList *vcfEntries, char *referenceSeq, int64_t refSeqLen, Params *params) {
    for (int64_t i = 0; i < stList_length(vcfEntries); i++) {
        VcfEntry *vcfEntry = stList_get(vcfEntries, i);
        assert(vcfEntry->alleleSubstrings == NULL);
        vcfEntry->alleleSubstrings = getAlleleSubstrings2(vcfEntry, referenceSeq, refSeqLen, &vcfEntry->refAlnStart,
                &vcfEntry->refAlnStopExcl, params->polishParams->columnAnchorTrim, params->polishParams->useRunLengthEncoding);

    }
}
