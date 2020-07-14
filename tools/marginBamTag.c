/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <memory.h>
#include <hashTableC.h>
#include <unistd.h>
#include <time.h>
#include "marginVersion.h"

#include "margin.h"
#include "htsIntegration.h"
#include "helenFeatures.h"


/*
 * Main functions
 */

void usage() {
    fprintf(stderr, "usage: marginTagBam <IN_BAM_FILE> <TAG_INFO_FILE> <OUT_BAM_FILE> \n");
    fprintf(stderr, "Version: %s \n\n", MARGIN_POLISH_VERSION_H);
    fprintf(stderr, "Tag IN_BAM_FILE with haplotype info in TAG_INFO_FILE, writing out to OUT_BAM_FILE.\n");

    fprintf(stderr, "\nRequired arguments:\n");
    fprintf(stderr, "    IN_BAM_FILE input file in BAM format\n");
    fprintf(stderr, "    TAG_INFO_FILE file in format \"read_id\\t[none|H1|H2]\\t...\n");
    fprintf(stderr, "    OUT_BAM_FILE output file in BAM format\n");
    fprintf(stderr, "\n");
}


int main(int argc, char *argv[]) {

    // Parameters / arguments
    char *bamInFile = NULL;
    char *readInfoFile = NULL;
    char *bamOutFile = NULL;

    if (argc < 4) {
        usage();
        return 0;
    }

    bamInFile = stString_copy(argv[1]);
    readInfoFile = stString_copy(argv[2]);
    bamOutFile = stString_copy(argv[3]);

    // for logging
    st_setLogLevel(info);

    // sanity check (verify files are accessible)
    if (access(bamInFile, R_OK) != 0) {
        st_errAbort("Could not read from input bam file: %s\n", bamInFile);
    }
    if (access(readInfoFile, R_OK) != 0) {
        st_errAbort("Could not read from read info file: %s\n", readInfoFile);
    }

    // get read info
    stHash *readInfo = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
    FILE *fp = fopen(readInfoFile, "r");
    if (fp == NULL) {
        st_errAbort("Could not open read info file %s\n", readInfoFile);
    }
    char *line = NULL;
    int linenr = -1;
    int hap1Count = 0;
    int hap2Count = 0;
    while ((line = stFile_getLineFromFile(fp)) != NULL) {
        linenr++;
        // header
        if (strlen(line) == 0 || line[0] == '#') {
            free(line);
            continue;
        }
        // get line parts
        stList *elements = stString_split(line);
        assert(stList_length(elements) >= 2);
        // read name
        char *readName = stString_copy(stList_get(elements, 0));
        // haplotag
        int64_t ht = -1;
        char *htInfo = stList_get(elements, 1);
        if (stString_eq(htInfo, "H1")) {
            ht = 1;
            hap1Count++;
        } else if (stString_eq(htInfo, "H2")) {
            ht = 2;
            hap2Count++;
        } else if (!stString_eq(htInfo, "none")) {
            st_errAbort("Unexpected haplotag descriptor, expected \"H1\", \"H2\", or \"none\": %s\n\tline %d: \"%s\"",
                    htInfo, linenr, line);
        }
        // save
        st_logInfo("  Saving %s with tag value %"PRId64"\n", readName, ht);
        stHash_insert(readInfo, readName, (void*)ht);
        // cleanup
        stList_destruct(elements);
        free(line);
    }
    fclose(fp);
    st_logCritical("Read %"PRId64" read haplotags, with %d H1 and %d H2\n", stHash_size(readInfo), hap1Count, hap2Count);


    // counting
    int64_t h1Count = 0;
    int64_t h2Count = 0;
    int64_t h0Count = 0;
    int64_t hUnCount = 0;

    // input file management
    samFile *in = hts_open(bamInFile, "r");
    st_logCritical("Reading from BAM: %s \n", bamInFile);
    hts_itr_multi_t *iter = NULL;
    if (in == NULL) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamInFile);
    }
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();
    int r;

    // open bam to write to
    st_logCritical("Writing to BAM: %s \n", bamOutFile);
    samFile *out = hts_open(bamOutFile, "wb");
    r = sam_hdr_write(out, bamHdr);

    // fetch alignments
    while (sam_read1(in,bamHdr,aln) >= 0) {
        // get haplotag value
        char *readName = bam_get_qname(aln);
        int32_t hp = 0;
        if (stHash_search(readInfo, readName) == NULL) {
            st_logInfo("Did not find %s in read info\n", readName);
            hUnCount++;
        } else {
            int64_t foundTag = (int64_t) stHash_search(readInfo, readName);
            if (foundTag == -1) {
                h0Count++;
            } else if (foundTag == 1) {
                hp = 1;
                h1Count++;
            } else if (foundTag == 2) {
                hp = 2;
                h2Count++;
            } else {
                st_errAbort("Got haplotag value %"PRId64" for read %s!", foundTag, readName);
            }
        }

        st_logInfo("Writing read %s with tag %"PRId32"\n", readName, hp);
        if (bam_aux_get(aln, "HP") != NULL) {
            bam_aux_update_int(aln, "HP", hp);
        } else {
            bam_aux_append(aln, "HP", 'i', sizeof(hp), (uint8_t*) &hp);
        }
        r = sam_write1(out, bamHdr, aln);
    }
    st_logCritical("Wrote reads with divisions: H1 %"PRId64", H2 %"PRId64", and H0 %"PRId64"\n",
               h1Count, h2Count, h0Count);
    st_logCritical("Found %"PRId64" reads which were not annotated in info file (tagged as H0, but not counted above).\n", hUnCount);

    // Cleanup
    bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
    sam_close(in);
    sam_close(out);
    stHash_destruct(readInfo);

    st_logInfo("Fin.");

    return 0;
}

