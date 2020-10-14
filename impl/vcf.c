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
    return vcfEntry;
}

void vcfEntry_destruct(VcfEntry *vcfEntry) {
    stList_destruct(vcfEntry->alleles);
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

stList *parseVcf(char *vcfFile, Params *params) {
    return parseVcf2(vcfFile, NULL, params);
}

stList *parseVcf2(char *vcfFile, char *regionStr, Params *params) {
    stList *entries = stList_construct3(0, (void(*)(void*))vcfEntry_destruct);

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
        stList_append(entries, entry);
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
                   stList_length(entries), params->phaseParams->includeHomozygousVCFEntries ? " " : "HET ", vcfFile,
                   skippedForRegion, skippedForNotPass, skippedForHomozygous, skippedForIndel);
    if (stList_length(entries) == 0) {
        st_errAbort("No valid VCF entries found!");
    }
    return entries;
}

stList *copyListOfRleStrings(stList *toCopy) {
    stList *copy = stList_construct3(0, (void (*)(void*)) rleString_destruct);
    for (int i = 0; i < stList_length(toCopy); i++) {
        stList_append(copy, rleString_copy(stList_get(toCopy, i)));
    }
    return copy;
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


stList *getVcfEntriesForRegion(stList *vcfEntries, uint64_t *rleMap, char *refSeqName, int64_t startPos, int64_t endPos) {
    return getVcfEntriesForRegion2(vcfEntries, rleMap, refSeqName, startPos, endPos, DEFAULT_MIN_VCF_QUAL);
}



void read_vcf(char *fname)
{
	//open vcf file
    htsFile *fp    = hts_open(fname,"rb");

    //read header
    bcf_hdr_t *hdr = bcf_hdr_read(fp);

    bcf1_t *rec    = bcf_init();

    //save for each vcf record
    while ( bcf_read(fp, hdr, rec)>=0 )
    {
    	//unpack for read REF,ALT,INFO,etc
        bcf_unpack(rec, BCF_UN_STR);
        bcf_unpack(rec, BCF_UN_INFO);

        //read CHROM
        printf("%s\n", bcf_hdr_id2name(hdr, rec->rid));

        //read POS
        printf("%lu\n", (unsigned long)rec->pos);

        //read ID(rsid)
        printf("%s\n", rec->d.id);

        //rec->d.allele[0] is REF other is ALT
        printf("%lu\n", (unsigned long)rec->n_allele);
        for (int i=0; i<rec->n_allele; ++i)
            printf("%s\n", rec->d.allele[i]);

        //check if is snp
        printf("ttt %d\n", bcf_is_snp(rec));
    }



    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ( (ret=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",fname,ret);
        exit(ret);
    }
}