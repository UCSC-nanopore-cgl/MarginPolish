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
    vcfEntry->alleleIdxToReads = stList_construct3(0, (void(*)(void*))stSet_destruct);
    for (int i = 0; i < stList_length(vcfEntry->alleles); i++) {
        stList_append(vcfEntry->alleleIdxToReads, stSet_construct());
    }
    vcfEntry->rootVcfEntry = NULL;
    vcfEntry->genotypeProb = -1.0;
    vcfEntry->haplotype1Prob = -1.0;
    vcfEntry->haplotype2Prob = -1.0;
    return vcfEntry;
}

void vcfEntry_destruct(VcfEntry *vcfEntry) {
    stList_destruct(vcfEntry->alleles);
    if (vcfEntry->alleleSubstrings != NULL) stList_destruct(vcfEntry->alleleSubstrings);
    if (vcfEntry->alleleIdxToReads != NULL) stList_destruct(vcfEntry->alleleIdxToReads);
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
        int gt1 = -1;
        int gt2 = -1;
        int32_t *gt_arr = NULL, ngt_arr = 0;
        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt>0 && !bcf_gt_is_missing(gt_arr[0])  && gt_arr[1] != bcf_int32_vector_end) {
            gt1 = bcf_gt_allele(gt_arr[0]);
            gt2 = bcf_gt_allele(gt_arr[1]);
        }
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

int64_t binarySearchVcfListForFirstIndexAtOrAfterRefPos2(stList *vcfEntries, int64_t desiredEntryPos, int64_t startPos,
                                                         int64_t endPosIncl) {
    if (endPosIncl - startPos == 1) {
        return ((VcfEntry*)stList_get(vcfEntries, startPos))->refPos >= desiredEntryPos ? startPos : endPosIncl;
    }
    int64_t middlePos = startPos + (endPosIncl - startPos) / 2;
    VcfEntry *middleEntry = stList_get(vcfEntries, middlePos);
    if (middleEntry->refPos < desiredEntryPos) {
        return binarySearchVcfListForFirstIndexAtOrAfterRefPos2(vcfEntries, desiredEntryPos, middlePos, endPosIncl);
    } else /*if (middleEntry->refPos >= desiredEntryPos)*/ {
        return binarySearchVcfListForFirstIndexAtOrAfterRefPos2(vcfEntries, desiredEntryPos, startPos, middlePos);
    }

}
int64_t binarySearchVcfListForFirstIndexAtOrAfterRefPos(stList *vcfEntries, int64_t refPos) {
    if (stList_length(vcfEntries) == 0) return -1;
    if (((VcfEntry*)stList_get(vcfEntries, stList_length(vcfEntries) - 1))->refPos < refPos) return -1;
    if (((VcfEntry*)stList_get(vcfEntries, 0))->refPos >= refPos) return 0;
    return binarySearchVcfListForFirstIndexAtOrAfterRefPos2(vcfEntries, refPos, 0, stList_length(vcfEntries) - 1);
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
    int64_t startIdx = binarySearchVcfListForFirstIndexAtOrAfterRefPos(vcfEntries, startPos);
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

    int64_t totalKeptEntries = stList_length(regionEntries);
    char *logIdentifier = getLogIdentifier();
    st_logInfo(" %s Filtered %"PRIu64" VCF records for quality < %.2f, kept %"PRId64" variants with quality < %.2f, totalling %"PRIu64" (every %"PRId64"bp).\n",
               logIdentifier, qualityFilteredCount, params->phaseParams->minVariantQuality,
               initiallyFilteredCount - stList_length(filteredEntries),
               params->phaseParams->variantSelectionAdaptiveSamplingPrimaryThreshold,  stList_length(regionEntries),
               (int64_t) ((endPos - startPos) / (totalKeptEntries == 0 ? -1 : totalKeptEntries)));

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


stHash *buildVcfEntryToBamChunkReadMap(stList *bamChunkReads) {
    // get hash of allele to list of substrings
    stHash *vcfEntriesToReads = stHash_construct2(NULL, (void(*)(void*))stList_destruct);
    for (int64_t i = 0; i < stList_length(bamChunkReads); i++) {
        BamChunkRead *bcr = stList_get(bamChunkReads, i);
        BamChunkReadVcfEntrySubstrings *bcrves = bcr->bamChunkReadVcfEntrySubstrings;
        for (uint64_t j = 0; j < stList_length(bcrves->vcfEntries); j++) {
            VcfEntry *vcfEntry = stList_get(bcrves->vcfEntries, j);

            // save
            stList *bcrList = stHash_search(vcfEntriesToReads, vcfEntry);
            if (bcrList == NULL) {
                bcrList = stList_construct();
                stHash_insert(vcfEntriesToReads, vcfEntry, bcrList);
            }
            stList_append(bcrList, bcr);
        }
    }

    return vcfEntriesToReads;
}


void updateOriginalVcfEntriesWithBubbleData(BamChunk *bamChunk, stList *bamChunkReads, stHash *readIdToIdx,
        stGenomeFragment *gF, BubbleGraph *bg, stList *chunkVcfEntriesToBubbles, stSet *hap1Reads, stSet *hap2Reads,
        char *logIdentifier) {

    // prep
    stHash *vcfEntryToReads = buildVcfEntryToBamChunkReadMap(bamChunkReads);

    // loop over all primary bubbles in the actual chunk boundaries
    for (uint64_t primaryBubbleIdx = 0; primaryBubbleIdx < gF->length; primaryBubbleIdx++) {
        // bubble and hap info
        Bubble *primaryBubble = &bg->bubbles[gF->refStart + primaryBubbleIdx];
        int64_t hap1AlleleNo = gF->haplotypeString1[primaryBubbleIdx];
        int64_t hap2AlleleNo = gF->haplotypeString2[primaryBubbleIdx];
        // TODO these probs are in log space: convert to [0-1]
        float genotypeProb = gF->genotypeProbs[primaryBubbleIdx];
        float haplotype1Prob = gF->haplotypeProbs1[primaryBubbleIdx];
        float haplotype2Prob = gF->haplotypeProbs2[primaryBubbleIdx];

        // vcf data
        VcfEntry *chunkVcfEntry = stList_get(chunkVcfEntriesToBubbles, primaryBubbleIdx);
        VcfEntry *rootVcfEntry = chunkVcfEntry->rootVcfEntry;

        // sanity checks
        assert(rootVcfEntry != NULL);
        assert(stList_length(chunkVcfEntry->alleles) == primaryBubble->alleleNo);

        // not in chunk
        if (rootVcfEntry->refPos < bamChunk->chunkStart || rootVcfEntry->refPos >= bamChunk->chunkEnd) {
            continue;
        }

        // get read substrings
        stList *entryBCRs = stHash_search(vcfEntryToReads, chunkVcfEntry);

        // nothing to phase with, make no updates
        if (stList_length(entryBCRs) == 0) {
            rootVcfEntry->gt1 = -1;
            rootVcfEntry->gt2 = -1;
            rootVcfEntry->genotypeProb = 0;
            rootVcfEntry->haplotype1Prob = 0;
            rootVcfEntry->haplotype2Prob = 0;
            continue;
        }

        // update genotypes and probs
        rootVcfEntry->gt1 = hap1AlleleNo;
        rootVcfEntry->gt2 = hap2AlleleNo;
        rootVcfEntry->genotypeProb = fromLog(genotypeProb);
        rootVcfEntry->haplotype1Prob = fromLog(haplotype1Prob);
        rootVcfEntry->haplotype2Prob = fromLog(haplotype2Prob);

        // update vcf read indices
        int64_t unMatchedReads = 0;
        stSet *hap1RootVcfEntryReadIndices = stList_get(rootVcfEntry->alleleIdxToReads, hap1AlleleNo);
        stSet *hap2RootVcfEntryReadIndices = stList_get(rootVcfEntry->alleleIdxToReads, hap2AlleleNo);
        for (int64_t i = 0; i < stList_length(entryBCRs); i++) {
            BamChunkRead *bcr = stList_get(entryBCRs, i);
            int64_t readIdx = (int64_t) stHash_search(readIdToIdx, bcr->readName);
            assert(readIdx != 0);
            if (stSet_search(hap1Reads, bcr) != NULL) {
                stSet_insert(hap1RootVcfEntryReadIndices, (void*) readIdx);
            } else if (stSet_search(hap2Reads, bcr) != NULL) {
                stSet_insert(hap2RootVcfEntryReadIndices, (void*) readIdx);
            } else {
                unMatchedReads++;
            }
        }
        if (unMatchedReads == stList_length(entryBCRs)) {
            st_logInfo(" %s No reads (out of %"PRId64") were aligned to VCF entry at pos %s:%"PRId64"\n",
                    logIdentifier, unMatchedReads, rootVcfEntry->refSeqName, rootVcfEntry->rawRefPosInformativeOnly);
        }

    }

    stHash_destruct(vcfEntryToReads);
}


void updateHaplotypeSwitchingInVcfEntries(BamChunker *chunker, bool *chunkWasSwitched, stHash *vcfEntryMap) {
    // trakcing
    int64_t totalSwitchedVcfEntries = 0;
    int64_t totalVcfEntries = 0;

    // for iteration
    char *currContig = NULL;
    stList *currVcfEntries = NULL;
    int64_t currVcfEntryIdx = 0;
    for (int64_t i = 0; i < chunker->chunkCount; i++) {
        BamChunk *chunk = stList_get(chunker->chunks, i);
        if (currContig == NULL || !stString_eq(currContig, chunk->refSeqName)) {
            currContig = chunk->refSeqName;
            currVcfEntries = stHash_search(vcfEntryMap, currContig);
            if (currVcfEntries == NULL) {
                // no entries in vcf entry map, we can continue
                continue;
            }
            currContig = chunk->refSeqName;
            currVcfEntryIdx = binarySearchVcfListForFirstIndexAtOrAfterRefPos(currVcfEntries, chunk->chunkStart);
            if (currVcfEntryIdx < 0) {
                // no entries for this contig
                currContig = NULL;
                continue;
            }
        }

        VcfEntry *vcfEntry;
        while (currVcfEntryIdx < stList_length(currVcfEntries) &&
                (vcfEntry = stList_get(currVcfEntries, currVcfEntryIdx))->refPos < chunk->chunkEnd) {
            if (vcfEntry->refPos < chunk->chunkStart) {
                // should not happen
                st_logInfo("  While switching haplotypes on VCF entries got entry starting before chunk (%s %"PRId64" < %"PRId64")\n",
                           vcfEntry->refSeqName, vcfEntry->refPos, chunk->chunkStart);
                currVcfEntryIdx++;
                continue;
            }

            // update
            if (chunkWasSwitched[i]) {
                int64_t tmpi = vcfEntry->gt1;
                vcfEntry->gt1 = vcfEntry->gt2;
                vcfEntry->gt2 = tmpi;
                float tmpf = vcfEntry->haplotype1Prob;
                vcfEntry->haplotype1Prob = vcfEntry->haplotype2Prob;
                vcfEntry->haplotype2Prob = tmpf;

                totalSwitchedVcfEntries++;
            }
            totalVcfEntries++;
            currVcfEntryIdx++;
        }
    }
    st_logInfo("  Switched %"PRId64"/%"PRId64" (%.2f) VCF entry haplotypes after stitching\n",
            totalSwitchedVcfEntries, totalVcfEntries, 1.0 * totalSwitchedVcfEntries / totalVcfEntries);
}


void writeUnphasedVariant(bcf_hdr_t *hdr, bcf1_t *rec, int32_t origGt1, int32_t origGt2) {
    int32_t *tmpia = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int));
    tmpia[0] = bcf_gt_unphased(origGt1);
    tmpia[1] = bcf_gt_unphased(origGt2);
    bcf_update_genotypes(hdr, rec, tmpia, 2);
    free(tmpia);
}

void recordPhaseSet(int32_t phaseSet, VcfEntry *prevHetVcfEntry, stList *phaseSetLengths, char *reasonStr, FILE *phaseSetBedOut) {
    if (phaseSet != -1 && prevHetVcfEntry != NULL) {
        int64_t psLength = prevHetVcfEntry->refPos - phaseSet;
        stList_append(phaseSetLengths, (void*) psLength);
        if (phaseSetBedOut != NULL) {
            fprintf(phaseSetBedOut, "%s\t%"PRId32"\t%"PRId64"\t%s\n", prevHetVcfEntry->refSeqName, phaseSet,
                    prevHetVcfEntry->refPos, reasonStr);
        }
    }
}

int cmpint64(const void *I, const void *J) {
    int64_t i = (int64_t) I;
    int64_t j = (int64_t) J;
    return i > j ? 1 : i < j ? -1 : 0;
}

void writePhasedVcf(char *inputVcfFile, char *regionStr, char *outputVcfFile, char *phaseSetBedFile,
        stHash *vcfEntryMap, Params *params) {

    //open files
    htsFile *fpIn = hts_open(inputVcfFile,"rb");
    if (fpIn == NULL) {
        st_logCritical("Could not open input VCF for reading %s\n", inputVcfFile);
        return;
    }
    htsFile *fpOut = hts_open(outputVcfFile, "w");
    if (fpOut == NULL) {
        st_logCritical("Could not open output VCF for writing %s\n", outputVcfFile);
        return;
    }
    FILE *phaseSetBedOut = phaseSetBedFile == NULL ? NULL : fopen(phaseSetBedFile, "w");
    if (phaseSetBedOut == NULL && phaseSetBedFile != NULL) {
        st_logCritical("Could not open phase set BED file for writing %s\n", phaseSetBedFile);
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

    //read header
    bcf_hdr_t *hdr = bcf_hdr_read(fpIn);
    int nsmpl = bcf_hdr_nsamples(hdr);
    if (nsmpl > 1) {
        st_logCritical("> Got %d samples reading %s, will only take VCF records for the first\n", nsmpl, inputVcfFile);
    }

    // ensure these are present
    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set Identifier\">");
    if (params->phaseParams->updateAllOutputVCFFormatFields) {
        bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
        bcf_hdr_append(hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
        bcf_hdr_append(hdr, "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">");
        bcf_hdr_append(hdr, "##FORMAT=<ID=HD,Number=2,Type=Integer,Description=\"Haplotype Depth\">");
        bcf_hdr_append(hdr, "##FORMAT=<ID=HCPV,Number=2,Type=Integer,Description=\"Haplotype Concordance with Previous Variant\">");
        bcf_hdr_append(hdr, "##FORMAT=<ID=HDPV,Number=2,Type=Integer,Description=\"Haplotype Discordance with Previous Variant\">");
    }


    // write header
    bcf_hdr_write(fpOut, hdr);
    bcf1_t *rec    = bcf_init();

    // tracking total entries
    int64_t totalEntries = 0;
    int64_t totalPhasedWritten = 0;
    int64_t skippedBecasueNotConsidered = 0;
    int64_t notPhasedBecauseMarginCalledHomozygous = 0;
    int64_t notPhasedBecauseMarginCalledHetDifferentFromInputVCF = 0;

    // tracking vcf entries
    VcfEntry *prevHetVcfEntry = NULL;
    VcfEntry *currVcfEntry = NULL;
    int32_t phaseSet = -1;
    int64_t nextVcfEntryIdx = 0;
    char *currChrom = NULL;
    stList *currChromVcfEntries = NULL;

    // phase set data
    stList *phaseSetLengths = stList_construct();

    // iterate
    while ( bcf_read(fpIn, hdr, rec) >= 0 )
    {
        //unpack for read REF,ALT,INFO,etc
        bcf_unpack(rec, BCF_UN_ALL);
        totalEntries++;

        // location data
        const char *chrom = bcf_hdr_id2name(hdr, rec->rid);
        int64_t pos = rec->pos;

        // skipped cases
        bool skipVariantAnalysis = FALSE;
        if (regionStr != NULL && (!stString_eq(regionContig, chrom) || (regionStart >= 0 && !(regionStart <= pos && pos < regionEnd)))) {
            skipVariantAnalysis = TRUE;
        }
        if (params->phaseParams->onlyUsePassVCFEntries && !bcf_has_filter(hdr, rec, "PASS")) {
            skipVariantAnalysis = TRUE;
        }
        if (params->phaseParams->onlyUseSNPVCFEntries && !bcf_is_snp(rec)) {
            skipVariantAnalysis = TRUE;
        }

        // genotype
        int32_t origGt1 = -1;
        int32_t origGt2 = -1;
        int32_t *gt_arr = NULL, ngt_arr = 0;
        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt>0 && !bcf_gt_is_missing(gt_arr[0])  && gt_arr[1] != bcf_int32_vector_end) {
            origGt1 = bcf_gt_allele(gt_arr[0]);
            origGt2 = bcf_gt_allele(gt_arr[1]);
        }
        free(gt_arr);
        if (!params->phaseParams->includeHomozygousVCFEntries && origGt1 == origGt2) {
            skipVariantAnalysis = TRUE;
        }

        // all skipped variants are written this way
        if (skipVariantAnalysis) {
            skippedBecasueNotConsidered++;
            writeUnphasedVariant(hdr, rec, origGt1, origGt2);
            continue;
        }

        // setup our vcf entries for new chrom
        if (currChrom == NULL || !stString_eq(currChrom, chrom)) {
            // handle phase set
            recordPhaseSet(phaseSet, prevHetVcfEntry, phaseSetLengths, "ContigEnd\t", phaseSetBedOut);
            // free old value
            if (currChrom != NULL) free(currChrom);
            // init chrom and entries
            currChrom = stString_copy(chrom);
            currChromVcfEntries = stHash_search(vcfEntryMap, currChrom);
            assert(currChromVcfEntries != NULL);
            // prep
            prevHetVcfEntry = NULL;
            currVcfEntry = NULL;
            nextVcfEntryIdx = 0;
            phaseSet = -1;
        }

        // find new current vcf entry
        VcfEntry *nextVcfEntry = NULL;
        int64_t skippedVcfEntries = 0;
        while (nextVcfEntryIdx < stList_length(currChromVcfEntries)) {
            nextVcfEntry = stList_get(currChromVcfEntries, nextVcfEntryIdx);
            if (nextVcfEntry->refPos == pos) {
                // found it
                nextVcfEntryIdx++;
                break;
            } else if (nextVcfEntry->refPos > pos) {
                //we have missed our variant, should not happen
                nextVcfEntry = NULL;
                break;
            } else if (nextVcfEntry->refPos < pos) {
                // we aren't at our variant yet, this should not happen
                skippedVcfEntries++;
            } else {
                assert(FALSE);
            }
            nextVcfEntryIdx++;
        }
        // this is only error case where something unexpected has happened
        if (skippedVcfEntries > 0) {
            st_logCritical("  Skipped %"PRId64" considered VCF entries searching %ssuccessfully for entry at %s:%"PRId64"\n",
                    skippedVcfEntries, nextVcfEntry == NULL ? "un":"", chrom, pos);
        }

        // handle case where we did not find this variant (bad)
        if (nextVcfEntry == NULL) {
            //loggit
            if (nextVcfEntryIdx < stList_length(currChromVcfEntries)) {
                nextVcfEntry = stList_get(currChromVcfEntries, nextVcfEntryIdx);
            }
            st_logCritical("  When writing VCF entries at %s:%"PRId64", did not find existing entry (prev %s:%"PRId64", next %s:%"PRId64")\n",
                       chrom, pos, prevHetVcfEntry != NULL ? prevHetVcfEntry->refSeqName : "NULL",
                       prevHetVcfEntry != NULL ? prevHetVcfEntry->refPos : -1, nextVcfEntry->refSeqName,
                       nextVcfEntry->refPos);
            // write variant
            writeUnphasedVariant(hdr, rec, origGt1, origGt2);
            continue;
        }

        // handle case where we found this variant, but it was for some reason filtered out (ok)
        if (nextVcfEntry->genotypeProb == -1.0) {
            skippedBecasueNotConsidered++;
            writeUnphasedVariant(hdr, rec, origGt1, origGt2);
            continue;
        }

        // iterate
        if (currVcfEntry != NULL && currVcfEntry->gt1 != currVcfEntry->gt2) {
            // prev must be het
            prevHetVcfEntry = currVcfEntry;
        }
        currVcfEntry = nextVcfEntry;

        // get variant data from margin analysis
        // genotype
        int gt1 = (int) currVcfEntry->gt1;
        int gt2 = (int) currVcfEntry->gt2;
        // probs
        int32_t gtProb = (int32_t) toPhred(currVcfEntry->genotypeProb);
        int32_t hp1Prob = (int32_t) toPhred(currVcfEntry->haplotype1Prob);
        int32_t hp2Prob = (int32_t) toPhred(currVcfEntry->haplotype2Prob);
        // depths
        int32_t depth = 0;
        int32_t hap1Depth = -1;
        int32_t hap2Depth = -1;
        for (int i = 0; i < stList_length(currVcfEntry->alleleIdxToReads); i++) {
            int32_t hpDepth = (int32_t) stSet_size(stList_get(currVcfEntry->alleleIdxToReads, i));
            depth += hpDepth;
            if (i == gt1) hap1Depth = hpDepth;
            if (i == gt2) hap2Depth = hpDepth;
        }
        // read concordancy with with previous het variant
        int32_t hcpv1 = -1;
        int32_t hcpv2 = -1;
        int32_t hdpv1 = -1;
        int32_t hdpv2 = -1;
        bool determinedHetConcordancy = FALSE;
        if (prevHetVcfEntry != NULL && gt1 != gt2 && prevHetVcfEntry->gt1 >= 0 && currVcfEntry->gt1 >= 0) {
            stSet *prevH1 = stList_get(prevHetVcfEntry->alleleIdxToReads, prevHetVcfEntry->gt1);
            stSet *prevH2 = stList_get(prevHetVcfEntry->alleleIdxToReads, prevHetVcfEntry->gt2);
            stSet *currH1 = stList_get(currVcfEntry->alleleIdxToReads, currVcfEntry->gt1);
            stSet *currH2 = stList_get(currVcfEntry->alleleIdxToReads, currVcfEntry->gt2);
            hcpv1 = (int32_t) stSet_sizeOfIntersection(prevH1, currH1);
            hcpv2 = (int32_t) stSet_sizeOfIntersection(prevH2, currH2);
            hdpv1 = (int32_t) stSet_sizeOfIntersection(prevH2, currH1);
            hdpv2 = (int32_t) stSet_sizeOfIntersection(prevH1, currH2);
            determinedHetConcordancy = TRUE;
        }

        // new phase set consideration
        bool newPhaseSet = FALSE;
        char *newPhaseSetReason = NULL;
        if (prevHetVcfEntry == NULL) {
            newPhaseSet = TRUE;
            st_logInfo("  Calling new phase set at %s:%"PRId64" because no previous HET\n", chrom, pos);
            newPhaseSetReason = stString_print("NoHet\t");
        } else if (determinedHetConcordancy) {
            //TODO switched to AND not OR for primary phasing detection
            if (hcpv1 == 0 && hcpv2 == 0) {
                newPhaseSet = TRUE;
                st_logInfo("  Calling new phase set at %s:%"PRId64" because missing concordancy (H1:%"PRId32", H2:%"PRId32")\n",
                        chrom, pos, hcpv1, hcpv2);
                newPhaseSetReason = stString_print("MissingConcordancy\tH1-%"PRId32"_H2-%"PRId32, hcpv1, hcpv2);
            } else if (binomialPValue(hcpv1 + hcpv2, hcpv1) < params->phaseParams->phasesetMinBinomialReadSplitLikelihood) {
                newPhaseSet = TRUE;
                st_logInfo("  Calling new phase set at %s:%"PRId64" because unlikely concordancy (H1:%"PRId32", H2:%"PRId32", prob:%.8f)\n",
                        chrom, pos, hcpv1, hcpv2, binomialPValue(hcpv1 + hcpv2, hcpv1));
                newPhaseSetReason = stString_print("UnlikelyConcordancy\tH1-%"PRId32"_H2-%"PRId32"_Prob-%.8f", hcpv1, hcpv2, binomialPValue(hcpv1 + hcpv2, hcpv1));
            } else if (1.0 * (hdpv1 + hdpv2) / (hcpv1 + hcpv2 + hdpv1 + hdpv2) > params->phaseParams->phasesetMaxDiscordantRatio) {
                newPhaseSet = TRUE;
                st_logInfo("  Calling new phase set at %s:%"PRId64" because of discordancy (H1D:%"PRId32"+H2D:%"PRId32" / H1C:%"PRId32"+H2C:%"PRId32"+H1D:%"PRId32"+H2D:%"PRId32" = %.4f)\n",
                        chrom, pos, hdpv1, hdpv2, hcpv1, hcpv2, hdpv1, hdpv2, 1.0 * (hdpv1 + hdpv2) / (hcpv1 + hcpv2 + hdpv1 + hdpv2));
                newPhaseSetReason = stString_print("Discordancy\tH1D-%"PRId32"_H2D-%"PRId32"_H1C-%"PRId32"_H2C-%"PRId32"_ratio-%.4f",
                        hcpv1, hcpv2, hdpv1, hdpv2, 1.0 * (hdpv1 + hdpv2) / (hcpv1 + hcpv2 + hdpv1 + hdpv2));
            }
        }

        if (newPhaseSet) {
            recordPhaseSet(phaseSet, prevHetVcfEntry, phaseSetLengths, newPhaseSetReason, phaseSetBedOut);
            free(newPhaseSetReason);
            phaseSet = (int32_t) pos;
        }
        bool writePhaseSet;
        if (gt1 != gt2) {
            writePhaseSet = TRUE;
        } else {
            writePhaseSet = FALSE;
            notPhasedBecauseMarginCalledHomozygous++;
        }
        bool isPhased = !newPhaseSet && writePhaseSet;

        // write values
        int32_t *tmpia = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int));
        if (params->phaseParams->updateAllOutputVCFFormatFields) {
            // write everything, it is ok to clobber existing data
            // write genotype
            if (isPhased) {
                tmpia[0] = gt1 < 0 ? bcf_gt_missing : bcf_gt_phased(gt1);
                tmpia[1] = gt1 < 0 ? bcf_gt_missing : bcf_gt_phased(gt2);
            } else {
                tmpia[0] = gt1 < 0 ? bcf_gt_missing : bcf_gt_unphased(gt1);
                tmpia[1] = gt1 < 0 ? bcf_gt_missing : bcf_gt_unphased(gt2);
            }
            bcf_update_genotypes(hdr, rec, tmpia, 2);
            // write quality info
            tmpia[0] = gtProb;
            bcf_update_format_int32(hdr, rec, "GQ", tmpia, 1);
            tmpia[0] = hp1Prob;
            tmpia[1] = hp2Prob;
            bcf_update_format_int32(hdr, rec, "HQ", tmpia, 2);
            // write depth info
            tmpia[0] = depth;
            bcf_update_format_int32(hdr, rec, "DP", tmpia, 1);
            tmpia[0] = hap1Depth;
            tmpia[1] = hap2Depth;
            bcf_update_format_int32(hdr, rec, "HD", tmpia, 2);
            // write read concordancy (only makes sense with het variants)
            if (gt1 != gt2) {
                tmpia[0] = hcpv1;
                tmpia[1] = hcpv2;
                bcf_update_format_int32(hdr, rec, "HCPV", tmpia, 2);
                tmpia[0] = hdpv1;
                tmpia[1] = hdpv2;
                bcf_update_format_int32(hdr, rec, "HDPV", tmpia, 2);
            }
        } else {
            // only write gt and phase set, not ok to clobber existing data
            // only update genotype (and phase set) if we match the called genotype
            if ( !( (gt1 == origGt1 && gt2 == origGt2) || (gt1 == origGt2 && gt2 == origGt1) ) ) {
                // we have not found the same genotypes as we originally got, phasing cannot be trusted
                isPhased = FALSE;
                writePhaseSet = FALSE;
                if (gt1 != gt2) {
                    notPhasedBecauseMarginCalledHetDifferentFromInputVCF++;
                }
            }

            // write GT, either phased with MP or unphased
            if (isPhased) {
                tmpia[0] = bcf_gt_phased(gt1);
                tmpia[1] = bcf_gt_phased(gt2);
            } else {
                tmpia[0] = bcf_gt_unphased(origGt1);
                tmpia[1] = bcf_gt_unphased(origGt1);
                assert(writePhaseSet == FALSE || newPhaseSet);
            }
            bcf_update_genotypes(hdr, rec, tmpia, 2);
        }

        // only update phase set on hets called by margin
        if (writePhaseSet) {
            tmpia[0] = phaseSet;
            bcf_update_format_int32(hdr, rec, "PS", tmpia, 1);
            totalPhasedWritten++;
        }

        // save it
        bcf_write(fpOut, hdr, rec);
        free(tmpia);
    }

    // loggit
    if (params->phaseParams->updateAllOutputVCFFormatFields) {
        st_logCritical("  Wrote %"PRId64" variants: %"PRId64" were phased; skipped %"PRId64" for not being analyzed, %"PRId64" for being homozygous\n",
                totalEntries, totalPhasedWritten, skippedBecasueNotConsidered, notPhasedBecauseMarginCalledHomozygous);
    } else {
        st_logCritical("  Wrote %"PRId64" variants: %"PRId64" were phased; skipped %"PRId64" for not being analyzed, %"PRId64" for being homozygous, %"PRId64" for disagreement with margin\n",
                       totalEntries, totalPhasedWritten, skippedBecasueNotConsidered,
                       notPhasedBecauseMarginCalledHomozygous, notPhasedBecauseMarginCalledHetDifferentFromInputVCF);
    }

    // finish phase sets
    recordPhaseSet(phaseSet, prevHetVcfEntry, phaseSetLengths, "ContigEnd\t", phaseSetBedOut);
    stList_sort(phaseSetLengths, cmpint64);
    int64_t minPhaseSetLen = INT64_MAX;
    int64_t maxPhaseSetLen = 0;
    int64_t totalPhaseSetLength = 0;
    for (int64_t i = 0; i < stList_length(phaseSetLengths); i++) {
        int64_t psLen = (int64_t) stList_get(phaseSetLengths, i);
        totalPhaseSetLength += psLen;
        minPhaseSetLen = psLen < minPhaseSetLen ? psLen : minPhaseSetLen;
        maxPhaseSetLen = maxPhaseSetLen < psLen ? psLen : maxPhaseSetLen;
    }
    int64_t avgPhaseSetLen = totalPhaseSetLength / stList_length(phaseSetLengths);
    int64_t n50PhaseSetLen = 0;
    int64_t n50tmp = 0;
    for (int64_t i = 0; i < stList_length(phaseSetLengths); i++) {
        int64_t psLen = (int64_t) stList_get(phaseSetLengths, i);
        n50tmp += psLen;
        if (n50tmp > totalPhaseSetLength / 2) {
            n50PhaseSetLen = psLen;
            break;
        }
    }
    st_logCritical("  Identified %"PRId64" phase sets with lengths avg:%"PRId64", min:%"PRId64", max:%"PRId64", N50:%"PRId64"\n",
            stList_length(phaseSetLengths), avgPhaseSetLen, minPhaseSetLen, maxPhaseSetLen, n50PhaseSetLen);


    // cleanup
    if (currChrom != NULL) free(currChrom);
    stList_destruct(phaseSetLengths);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ( (ret=hts_close(fpIn)) ) {
        st_logCritical("  Failed to close input VCF %s with code %d!\n", inputVcfFile, ret);
    }
    if ( (ret=hts_close(fpOut)) ) {
        st_logCritical("  Failed to close output VCF %s with code %d!\n", outputVcfFile, ret);
    }
    if (phaseSetBedOut != NULL) {
        fclose(phaseSetBedOut);
    }
}
