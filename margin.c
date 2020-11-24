/*
 * Copyright (C) 2020 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <stdio.h>
#include <string.h>
#include "marginVersion.h"

int polish_main(int argc, char *argv[]);
int phase_main(int argc, char *argv[]);

void usage() {
    fprintf(stderr, "Program: margin (tools for analysis of long read data)\n");
    fprintf(stderr, "Version: %s \n\n", MARGIN_POLISH_VERSION_H);

    fprintf(stderr, "usage: margin <command> [options]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Commands:\n");
    fprintf(stderr, "    polish             Polishes a reference sequence using read data\n");
    fprintf(stderr, "    phase              Haplotags reads and phases variants using read data and VCF\n");
    fprintf(stderr, "\n");
}


int main(int argc, char *argv[]) {

    if (argc < 2 || strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
        usage();
        return argc < 2 ? 1 : 0;
    }

    if (strcmp(argv[1], "polish") == 0) {
        return polish_main(argc-1, argv+1);
    } else if (strcmp(argv[1], "phase") == 0) {
        return phase_main(argc-1, argv+1);
    } else if (strcmp(argv[1], "version") == 0) {
        fprintf(stderr, "%s\n", MARGIN_POLISH_VERSION_H);
        return 0;
    } else {
        usage();
        fprintf(stderr, "\nunrecognized command '%s'\n", argv[1]);
        return 1;
    }
}
