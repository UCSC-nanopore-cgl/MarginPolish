#!/usr/bin/env python3
import argparse
import sys
import json
import glob
import matplotlib.pyplot as plt
import pysam
import collections
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.stats import pearsonr


def parse_args(args = None):
    parser = argparse.ArgumentParser("Compares phasing for reads haplotyped by margin")
    parser.add_argument('--input', '-i', dest='input', default=None, required=True, type=str,
                        help='JSON file describing bubbles and phasing ')
    parser.add_argument('--truth_hap1', '-1', dest='truth_hap1', default=None, required=True, type=str,
                        help='Truth Hap1 readset')
    parser.add_argument('--truth_hap2', '-2', dest='truth_hap2', default=None, required=True, type=str,
                        help='Truth Hap2 readset')

    parser.add_argument('--het_vcf', '-v', dest='het_vcf', default=None, required=False, type=str,
                        help='VCF containing only heterozygous sites')
    parser.add_argument('--result_vcf', '-r', dest='result_vcf', default=None, required=False, type=str,
                        help='VCF containing TP/FP/FN classified sites')
    parser.add_argument('--chunks', '-c', dest='chunks', default=None, required=False, type=str,
                        help='File describing chunk positions')
    parser.add_argument('--figure_name', '-f', dest='figure_name', default=None, required=False, type=str,
                        help='Figure name (will save if set)')
    parser.add_argument('--untagged_only', '-u', dest='untagged_only', default=False, required=False, action='store_true',
                        help='Only plot untagged reads')

    return parser.parse_args() if args is None else parser.parse_args(args)

def log(msg):
    print(msg, file=sys.stderr)


PRIMARY="primary"
FILTERED="filtered"
READS="reads"
NAME="name"
STRAND="strand"
QUAL="qual"
HAP="hap"
HAP_SUPPORT="hapSupport"
REF_POS="refPos"
STRAND_SKEW="strandSkew"

def main():
    args = parse_args()

    # get json file
    with open(args.input, 'r') as istream:
        json_doc = json.load(istream)

    # get truth reads
    truth_h1 = set()
    truth_h2 = set()
    with open(args.truth_hap1, 'r') as fin:
        for line in fin:
            truth_h1.add(line.split(',')[0].strip())
    with open(args.truth_hap2, 'r') as fin:
        for line in fin:
            truth_h2.add(line.split(',')[0].strip())
    log("Found {} truth H1 reads and {} truth H2 reads".format(len(truth_h1), len(truth_h2)))

    # parse json data
    primary_reads = set()
    filtered_reads = set()
    read_data = {}
    for bubble in json_doc[PRIMARY]:
        for read in bubble[READS]:
            primary_reads.add(read[NAME])
    for bubble in json_doc[FILTERED]:
        for read in bubble[READS]:
            filtered_reads.add(read[NAME])
    for read in json_doc[READS]:
        read[PRIMARY] = read[NAME] in primary_reads
        read[FILTERED] = read[NAME] in filtered_reads
        # save read
        read_data[read[NAME]] = read

    print(read_data)





if __name__ == "__main__":
    main()

