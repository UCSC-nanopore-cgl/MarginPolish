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
    vcfEntry->quality = phredQuality;
    vcfEntry->alleles = alleles == NULL ? stList_construct3(0, (void(*)(void*))rleString_destruct) : alleles;
    vcfEntry->gt1 = gt1;
    vcfEntry->gt2 = gt2;
    vcfEntry->alleleSubstrings = NULL;
    vcfEntry->refAlnStart = -1;
    vcfEntry->refAlnStopIncl = -1;
    vcfEntry->alleleIdxToReads = stList_construct3(0, (void(*)(void*))stList_destruct);
    for (int i = 0; i < stList_length(vcfEntry->alleles); i++) {
        stList_append(vcfEntry->alleleIdxToReads, stList_construct());
    }
    vcfEntry->rootVcfEntry = NULL;
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

int vcfEntry_positionCmp(const void *a, const void *b) {
    VcfEntry *A = (VcfEntry*) a;
    VcfEntry *B = (VcfEntry*) b;
    if (A->refPos == B->refPos) return 0;
    return A->refPos < B->refPos ? -1 : 1;
}

int vcfEntry_qualityCmp(const void *a, const void *b) {
    VcfEntry *A = (VcfEntry*) a;
    VcfEntry *B = (VcfEntry*) b;
    if (A->quality == B->quality) return 0;
    return A->quality < B->quality ? -1 : 1; //ascending order, and we pop off end
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

        double quality = rec->qual;

        // get alleles
        stList *alleles = stList_construct3(0, (void (*)(void*)) rleString_destruct);
        for (int i=0; i<rec->n_allele; ++i) {
            stList_append(alleles,
                    params->polishParams->useRunLengthEncoding ?
                    rleString_construct(rec->d.allele[i]) :
                    rleString_construct_no_rle(rec->d.allele[i]));
        }

        // save it
        VcfEntry *entry = vcfEntry_construct(chrom, pos, pos, quality, alleles, gt1, gt2);
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
        stList_sort(contigEntries, vcfEntry_positionCmp);
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

stList *getVcfEntriesForRegion(stHash *vcfEntryMap, uint64_t *rleMap, char *refSeqName, int64_t startPos,
        int64_t endPos, Params *params) {
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
    stList *filteredEntries = stList_construct3(0, (void(*)(void*))vcfEntry_destruct);

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
        if (params->phaseParams->minVariantQuality > e->quality) {
            qualityFilteredCount++;
            continue;
        }
        int64_t refPos = e->refPos - startPos + 1; // e->refPos is in 0-based space, poa is 1-based. convert here
        if (rleMap != NULL) {
            refPos = rleMap[refPos];
        }

        // make variant
        VcfEntry *copy = vcfEntry_construct(e->refSeqName, refPos, e->rawRefPosInformativeOnly, e->quality,
                copyListOfRleStrings(e->alleles), e->gt1, e->gt2);
        copy->rootVcfEntry = e;

        // where to save it?
        if (params->phaseParams->useVariantSelectionAdaptiveSampling &&
            e->quality < params->phaseParams->variantSelectionAdaptiveSamplingPrimaryThreshold) {
            stList_append(filteredEntries, copy);
        } else {
            stList_append(regionEntries, copy);
        }
    }

    // do we need to keep some filtered variants?
    int64_t initiallyFilteredCount = stList_length(filteredEntries);
    int64_t currentVariantCount = stList_length(regionEntries);
    int64_t desiredVariantCount = (endPos - startPos) /
                                  params->phaseParams->variantSelectionAdaptiveSamplingDesiredBasepairsPerVariant;
    if (params->phaseParams->useVariantSelectionAdaptiveSampling && currentVariantCount < desiredVariantCount) {
        // shuffle so that ties are not broken by position
        stList_shuffle(filteredEntries);
        // sort by quality descending
        stList_sort(filteredEntries, vcfEntry_qualityCmp);

        // take from top until we have desired amount or are empty
        int64_t currentFilteredCount = initiallyFilteredCount;
        while (currentFilteredCount > 0 && currentVariantCount < desiredVariantCount) {
            VcfEntry *entry = stList_pop(filteredEntries);
            stList_append(regionEntries, entry);
            currentFilteredCount--;
            currentVariantCount++;
        }

        // resort entries
        stList_sort(regionEntries, vcfEntry_positionCmp);
    }

    char *logIdentifier = getLogIdentifier();
    st_logInfo(" %s Filtered %"PRIu64" VCF records for quality < %.2f, kept %"PRId64" variants with quality < %.2f, totalling %"PRIu64" (every %"PRId64"bp).\n",
               logIdentifier, qualityFilteredCount, params->phaseParams->minVariantQuality,
               initiallyFilteredCount - stList_length(filteredEntries),
               params->phaseParams->variantSelectionAdaptiveSamplingPrimaryThreshold,  stList_length(regionEntries),
               (int64_t) ((endPos - startPos) / stList_length(regionEntries)));

    stList_destruct(filteredEntries);
    free(logIdentifier);
    return regionEntries;
}


stList *getAlleleSubstrings2(VcfEntry *entry, char *referenceSeq, int64_t refSeqLen, int64_t *refStartPos,
        int64_t *refEndPosIncl, bool putRefPosInPOASpace, int64_t expansion, bool useRunLengthEncoding) {
    stList *substrings = stList_construct3(0, (void (*)(void*)) rleString_destruct);
    assert(stList_length(entry->alleles) >= 1);

    // parameters for substringing
    int64_t pos = entry->refPos;
    // at this point, reference positions are in poa-space (1-based), need to convert to ref-space (0-based)
    pos--;

    //get ref info
    char *refAllele = rleString_expand(stList_get(entry->alleles, 0));
    int64_t refAlleleLen = strlen(refAllele);
    for (int64_t i = 0; i < refAlleleLen; i++) {
        // if we have an insert at the end of a chunk
        if (pos+i >= refSeqLen) {
            assert(pos+1 == refSeqLen);
            refAlleleLen = i;
            break;
        }
        char refChar = referenceSeq[pos + i];
        char alleleChar = refAllele[i];
        assert(refChar == alleleChar);
    }
    free(refAllele);

    // get prefix and suffix strings (from ref seq)
    int64_t pStart = pos - expansion;
    int64_t sStart = pos + refAlleleLen;
    int64_t sLen = sStart + expansion >= refSeqLen ? refSeqLen - sStart : expansion;
    if (sStart >= refSeqLen) {
        assert(sStart == refSeqLen);
        sStart = refSeqLen - 1;
        sLen = 0;
    }
    *refStartPos = pStart < 0 ? 0 : pStart;
    *refEndPosIncl = sStart + expansion >= refSeqLen ? refSeqLen - 1 : sStart + expansion;

    char *prefix = stString_getSubString(referenceSeq, *refStartPos, pStart < 0 ? pos : expansion);
    char *suffix = stString_getSubString(referenceSeq, sStart, sLen);

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
    if (putRefPosInPOASpace) {
        (*refStartPos)++;
        (*refEndPosIncl)++;
    }
    return substrings;
}

stList *getAlleleSubstrings(VcfEntry *entry, RleString *referenceSeq, Params *params,
        int64_t *refStartPos, int64_t *refEndPosIncl, bool refPosInPOASpace) {
    char *rawRefSeq = rleString_expand(referenceSeq);
    stList *alleleSubstrings = getAlleleSubstrings2(entry, rawRefSeq, referenceSeq->nonRleLength, refStartPos, refEndPosIncl,
                                                    refPosInPOASpace, params->polishParams->columnAnchorTrim,
                                                    params->polishParams->useRunLengthEncoding);
    free(rawRefSeq);
    return alleleSubstrings;

}

void updateVcfEntriesWithSubstringsAndPositions(stList *vcfEntries, char *referenceSeq, int64_t refSeqLen,
                                                bool refPosInPOASpace, Params *params) {
    for (int64_t i = 0; i < stList_length(vcfEntries); i++) {
        VcfEntry *vcfEntry = stList_get(vcfEntries, i);
        assert(vcfEntry->alleleSubstrings == NULL);
        vcfEntry->alleleSubstrings = getAlleleSubstrings2(vcfEntry, referenceSeq, refSeqLen, &vcfEntry->refAlnStart,
                &vcfEntry->refAlnStopIncl, refPosInPOASpace, params->polishParams->columnAnchorTrim,
                params->polishParams->useRunLengthEncoding);

    }
}

//TODO
/*BubbleGraph *updateOriginalVcfEntriesWithBubbleData(BamChunk *bamChunk, stList *bamChunkReads, stHash *readIdToIdx,
        stGenomeFragment *gF, BubbleGraph *bg, stList *chunkVcfEntriesToBubbles, stSet *hap1Reads, stSet *hap2Reads,
        Params *params, char *logIdentifier) {

    // prep
    stHash *vcfEntryToReadSubstrings = buildVcfEntryToReadSubstringsMap(bamChunkReads, params);

    // loop over all primary bubbles
    for (uint64_t primaryBubbleIdx = 0; primaryBubbleIdx < gF->length; primaryBubbleIdx++) {

        // bubble and hap info
        Bubble *primaryBubble = &bg->bubbles[gF->refStart + primaryBubbleIdx];
        int64_t hap1AlleleNo = gF->haplotypeString1[primaryBubbleIdx];
        int64_t hap2AlleleNo = gF->haplotypeString2[primaryBubbleIdx];

        VcfEntry *chunkVcfEntry = stList_get(chunkVcfEntriesToBubbles, primaryBubbleIdx);
        VcfEntry *rootVcfEntry = chunkVcfEntry->rootVcfEntry;
        assert(rootVcfEntry != NULL);

        char *hap1 = rleString_expand(primaryBubble->alleles[hap1AlleleNo]);
        char *hap2 = rleString_expand(primaryBubble->alleles[hap2AlleleNo]);

        // anchor positions
        uint64_t refStart = primaryBubble->refStart;

        // make allele list from primary haplotype alleles
        stList *alleles = stList_construct3(0, free);
        stList_append(alleles, rleString_expand(hap1));
        stList_append(alleles, rleString_expand(hap2));

        // get read substrings
        stList *readSubstrings = stHash_search(vcfEntryToReadSubstrings, vcfEntry);

        // nothing to phase with
        if(stList_length(readSubstrings) == 0) {
            stList_destruct(readSubstrings);
            continue;
        }

        Bubble *b = st_malloc(sizeof(Bubble)); // Make a bubble and add to list of bubbles
        b->variantPositionOffsets = NULL;

        // Set the coordinates
        b->refStart = (uint64_t) refStart;

        // The reference allele
        b->refAllele = params->polishParams->useRunLengthEncoding ?
                       rleString_construct(stList_get(vcfEntry->alleles, 0)) :
                       rleString_construct_no_rle(stList_get(vcfEntry->alleles, 0));

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
            b->alleles[j] = params->polishParams->useRunLengthEncoding ?
                            rleString_construct(stList_get(alleles, j)) : rleString_construct_no_rle(stList_get(alleles, j));
        }

        // Get allele supports
        b->alleleReadSupports = st_calloc(b->readNo * b->alleleNo, sizeof(float));

        stList *anchorPairs = stList_construct(); // Currently empty

        SymbolString alleleSymbolStrings[b->alleleNo];
        for (int64_t j = 0; j < b->alleleNo; j++) {
            alleleSymbolStrings[j] = rleString_constructSymbolString(b->alleles[j], 0,
                                                                     b->alleles[j]->length,
                                                                     params->polishParams->alphabet,
                                                                     params->polishParams->useRepeatCountsInAlignment,
                                                                     maximumRepeatLengthExcl);
        }

        // get alignment likelihoods
        stHash *cachedScores = stHash_construct3(rleString_stringKey, rleString_expandedStringEqualKey,
                                                 (void (*)(void *)) rleString_destruct, free);
        for (int64_t k = 0; k < b->readNo; k++) {
            RleString *readSubstring = bamChunkReadSubstring_getRleString(b->reads[k]);
            SymbolString rS = rleString_constructSymbolString(readSubstring, 0, readSubstring->length,
                                                              params->polishParams->alphabet,
                                                              params->polishParams->useRepeatCountsInAlignment,
                                                              maximumRepeatLengthExcl);
            StateMachine *sM = b->reads[k]->read->forwardStrand
                               ? params->polishParams->stateMachineForForwardStrandRead
                               : params->polishParams->stateMachineForReverseStrandRead;

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
                            (float) computeForwardProbability(alleleSymbolStrings[j], rS, anchorPairs,
                                                              params->polishParams->p, sM, 0, 0);
                }
            }

            symbolString_destruct(rS);
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

        }

        // cleanup
        stHash_destruct(cachedScores);
        for (int64_t j = 0; j < b->alleleNo; j++) {
            symbolString_destruct(alleleSymbolStrings[j]);
        }
        stList_destruct(anchorPairs);
        stList_destruct(alleles);
        bubble_destruct(*b);
        free(b);
    }

    // get scores and save to appropriate sets
    double totalNoScoreVariantsSpanned = 0.0;
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
                totalNoScoreVariantsSpanned += stList_length(bcr->bamChunkReadVcfEntrySubstrings->vcfEntries);
                noScoreCount++;
            }
            unclassifiedCount++;
        }
    }

    // loggit
    int64_t length = stList_length(bamChunkReads);
    st_logInfo(" %s Of %"PRId64" filtered reads: %"PRId64" (%.2f) were hap1, %"PRId64" (%.2f) were hap2, %"PRId64" (%.2f) were unclassified with %"PRId64" (%.2f) having no score (avg spanned variants %.2f).\n",
               logIdentifier, length, hap1Count, 1.0*hap1Count/length, hap2Count, 1.0*hap2Count/length,
               unclassifiedCount, 1.0*unclassifiedCount/length, noScoreCount,
               1.0*noScoreCount/(unclassifiedCount == 0 ? 1 : unclassifiedCount),
               totalNoScoreVariantsSpanned / (noScoreCount == 0 ? 1 : noScoreCount));


    // other cleanup
    stHash_destruct(totalReadScore_hap1);
    stHash_destruct(totalReadScore_hap2);
    stHash_destruct(vcfEntryToReadSubstrings);
}*/
