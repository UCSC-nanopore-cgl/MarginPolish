/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"
#include "bedidx.h"

#define DEFAULT_ALIGNMENT_SCORE 10

BamChunkRead *bamChunkRead_construct() {
    return bamChunkRead_construct2(NULL, NULL, NULL, TRUE, NULL);
}

BamChunkRead *bamChunkRead_construct2(char *readName, char *nucleotides,
                                      uint8_t *qualities, bool forwardStrand, bool useRunLengthEncoding) {
    return bamChunkRead_construct3(readName, nucleotides, qualities, forwardStrand, 0, useRunLengthEncoding);
}

BamChunkRead *bamChunkRead_construct3(char *readName, char *nucleotides, uint8_t *qualities, bool forwardStrand,
                                      int64_t fullReadLength, bool useRunLengthEncoding) {
    BamChunkRead *r = calloc(1, sizeof(BamChunkRead));
    r->readName = stString_copy(readName);
    r->forwardStrand = forwardStrand;
    r->fullReadLength = fullReadLength;
    assert(nucleotides != NULL);
    r->rleRead = useRunLengthEncoding ? rleString_construct(nucleotides) : rleString_construct_no_rle(nucleotides);
    if (qualities != NULL) {
        r->qualities = rleString_rleQualities(r->rleRead, qualities);
    }

    return r;
}

BamChunkRead *bamChunkRead_constructCopy(BamChunkRead *copy) {
    BamChunkRead *r = calloc(1, sizeof(BamChunkRead));
    r->readName = stString_copy(copy->readName);
    r->forwardStrand = copy->forwardStrand;
    r->fullReadLength = copy->fullReadLength;
    r->rleRead = rleString_copy(copy->rleRead);
    if (copy->qualities != NULL) {
        r->qualities = st_calloc(r->rleRead->length, sizeof(uint8_t));
        for (int64_t i = 0; i < r->rleRead->length; i++) {
            r->qualities[i] = copy->qualities[i];
        }
    }

    return r;
}

void bamChunkRead_destruct(BamChunkRead *r) {
    if (r->readName != NULL) free(r->readName);
    if (r->rleRead != NULL) rleString_destruct(r->rleRead);
    if (r->qualities != NULL) free(r->qualities);
    free(r);
}

BamChunkReadVcfEntrySubstrings *bamChunkReadVcfEntrySubstrings_construct(char *readName, bool forwardStrand,
                                                                         int64_t fullReadLength) {
    return bamChunkReadVcfEntrySubstrings_construct2(readName, forwardStrand, fullReadLength,
                                                     stList_construct3(0, free),
                                                     stList_construct3(0, free), stList_construct());
}

BamChunkReadVcfEntrySubstrings *bamChunkReadVcfEntrySubstrings_construct2(char *readName, bool forwardStrand,
                                                                          int64_t fullReadLength,
                                                                          stList *readSubstrings,
                                                                          stList *readSubstringQualities,
                                                                          stList *vcfEntries) {
    BamChunkReadVcfEntrySubstrings *bcrs = calloc(1, sizeof(BamChunkReadVcfEntrySubstrings));
    bcrs->readName = stString_copy(readName);
    bcrs->forwardStrand = forwardStrand;
    bcrs->fullReadLength = fullReadLength;
    bcrs->readSubstrings = readSubstrings;
    bcrs->readSubstringQualities = readSubstringQualities;
    bcrs->vcfEntries = vcfEntries;
    return bcrs;
}

void bamChunkReadVcfEntrySubstrings_destruct(BamChunkReadVcfEntrySubstrings *bcrs) {
    stList_destruct(bcrs->readSubstrings);
    stList_destruct(bcrs->vcfEntries);
    free(bcrs);
}

void bamChunkReadVcfEntrySubstrings_saveSubstring(BamChunkReadVcfEntrySubstrings *bcrs, char *read, uint8_t *qualities,
                                                  VcfEntry *vcfEntry) {
    stList_append(bcrs->readSubstrings, read);
    stList_append(bcrs->vcfEntries, vcfEntry);
    // ensure we can free, even if these are null
    stList_append(bcrs->readSubstringQualities, qualities == NULL ? st_calloc(1, sizeof(uint8_t)) : qualities);
}
