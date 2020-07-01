#!/usr/bin/env python3
import argparse
import sys
import json
import glob
import matplotlib.pyplot as plt
import math
from collections import defaultdict
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
    parser.add_argument('--spacer_size', '-z', dest='spacer_size', default=500, required=False, type=int,
                        help='Size in BP to plot spacers between bubble sites')
    parser.add_argument('--verbose', '-v', dest='verbose', default=False, required=False, action='store_true',
                        help='Print extra information on plot')
    parser.add_argument('--figure_name', '-f', dest='figure_name', default=None, required=False, type=str,
                        help='Figure name (will save if set)')

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
HAP_SUPPORT_H1="hapSupportH1"
HAP_SUPPORT_H2="hapSupportH2"
STRAND_SKEW="strandSkew"
REF_POS="refPos"
START_POS="startPos"
END_POS="endPos"
TYPE="type"


def log_filtered_read_types(read_data, xpos_list, primary_reads, filtered_reads, filtered_bubbles_merged):

    xpos_type_hap_values = defaultdict(lambda: defaultdict(lambda: [[],[],[]]))

    # write points to plot
    for read in read_data.values():
        name = read[NAME]

        # then iterate over sites (xpos)
        for x,p in enumerate(xpos_list):
            # only plot points where reads overlap
            if read[START_POS] <= p < read[END_POS]:

                # primary reads
                if name in primary_reads:
                    continue

                # filtered reads
                elif name in filtered_reads:
                    if p in filtered_bubbles_merged and name in filtered_bubbles_merged[p][READS]:
                        bub_read_summary = filtered_bubbles_merged[p][READS][name]
                        for bub_read in bub_read_summary[READS]:
                            xpos_type_hap_values[p][bub_read[TYPE]][bub_read[HAP]].append(bub_read[HAP_SUPPORT])

    for pos in sorted(list(xpos_type_hap_values.keys())):
        type_hap_values = xpos_type_hap_values[pos]
        pos_total = sum(map(lambda x: sum(x[1]) + sum(x[2]), type_hap_values.values()))
        hap1_total = sum(map(lambda x: sum(x[1]), type_hap_values.values()))
        hap2_total = sum(map(lambda x: sum(x[2]), type_hap_values.values()))
        log("{}".format(pos))
        log("  H1: {:11d} ({:.3f})".format(int(hap1_total), hap1_total / pos_total))
        log("  H2: {:11d} ({:.3f})".format(int(hap2_total), hap2_total / pos_total))
        for type in sorted(list(xpos_type_hap_values[pos].keys())):
            hap1_values = type_hap_values[type][1]
            hap2_values = type_hap_values[type][2]
            hap1_total = sum(hap1_values)
            hap2_total = sum(hap2_values)
            type_total = max(1, hap1_total + hap2_total)
            log("    {}: {:9d} ({:.3f})".format(type.upper(), int(type_total), type_total / pos_total))
            log("        H1:  total {:9d} ({:.3f}), avg {:9d}".format(int(hap1_total), hap1_total/type_total,
                                                                    0 if len(hap1_values) == 0 else int(np.mean(hap1_values))))
            log("        H2:  total {:9d} ({:.3f}), avg {:9d}".format(int(hap2_total), hap2_total/type_total,
                                                                    0 if len(hap2_values) == 0 else int(np.mean(hap2_values))))



def main():
    args = parse_args()
    ss = args.spacer_size

    # get json file
    log("Parsing {} for read and phasing info".format(args.input))
    with open(args.input, 'r') as istream:
        json_doc = json.load(istream)
    if FILTERED not in json_doc:
        json_doc[FILTERED] = []

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
    primary_bubbles = dict()
    filtered_bubbles = dict()
    primary_reads = set()
    filtered_reads = set()
    hap0_reads = set()
    hap1_reads = set()
    hap2_reads = set()
    read_data = {}
    for bubble in json_doc[PRIMARY]:
        primary_bubbles[bubble[REF_POS]] = bubble
        bubble_reads = dict()
        for read in bubble[READS]:
            primary_reads.add(read[NAME])
            assert(read[NAME] not in bubble_reads)
            bubble_reads[read[NAME]] = read
        bubble[READS] = bubble_reads
    for bubble in json_doc[FILTERED]:
        pos=bubble[REF_POS]
        assert(pos not in filtered_bubbles)
        filtered_bubbles[pos] = bubble
        bubble_reads = dict()
        for read in bubble[READS]:
            assert(read[NAME] not in bubble_reads)
            filtered_reads.add(read[NAME])
            bubble_reads[read[NAME]] = read
        bubble[READS] = bubble_reads
    for read in json_doc[READS]:
        read[PRIMARY] = read[NAME] in primary_reads
        read[FILTERED] = read[NAME] in filtered_reads
        if read[HAP] == 1: hap1_reads.add(read[NAME])
        if read[HAP] == 2: hap2_reads.add(read[NAME])
        else: hap0_reads.add(read[NAME])
        # save read
        read_data[read[NAME]] = read
    log("Parsed {} primary bubbles, {} filtered bubbles, and {} reads".format(len(json_doc[PRIMARY]), len(json_doc[FILTERED]),
                                                                      len(read_data)))
    prev_primary_read_count = len(primary_reads)
    primary_reads = primary_reads.difference(filtered_reads)
    log("Removed {} filtered reads from primary set, currently {} primary and {} filtered.".format(
        prev_primary_read_count - len(primary_reads), len(primary_reads), len(filtered_reads)))

    # figure out "correctness", potentially switch haplotypes
    h1_th1 = len(hap1_reads.intersection(truth_h1))
    h1_th2 = len(hap1_reads.intersection(truth_h2))
    h2_th1 = len(hap2_reads.intersection(truth_h1))
    h2_th2 = len(hap2_reads.intersection(truth_h2))
    right = h1_th1 + h2_th2
    rong = h1_th2 + h2_th1
    total = right + rong
    log("Got {} ({:.3f}) right and {} ({:.3f}) wrong haplotyped reads, {}.".format(right, right/total, rong, rong/total,
        "maintaining orientation" if right > rong else "switching orientation"))
    if right < rong:
        tmp = truth_h1
        truth_h1 = truth_h2
        truth_h2 = tmp

    # get ordering of x positions based on primary bubbles and the spacing size
    primary_bubble_positions = list(primary_bubbles.keys())
    primary_bubble_positions.sort()
    xpos_list = []
    xpos_ticks = []
    xpos_labels = []
    idx = 0
    prev_p = None
    for p in primary_bubble_positions:
        if prev_p is None:
            xpos_list.append(p-ss)
            idx += 1
        else:
            prev_dist = p - prev_p
            prev_dist_increment = int(float(prev_dist) / (math.ceil(float(prev_dist) / ss)) + 1)
            prev_p += prev_dist_increment
            while prev_p < p:
                xpos_list.append(prev_p)
                idx += 1
                prev_p += prev_dist_increment
        xpos_ticks.append(idx)
        xpos_list.append(p)
        xpos_labels.append(p)
        idx += 1
        prev_p = p
    xpos_list.append(prev_p + ss)

    # how to find the closest based-on-primary-bubbles positions
    def closest_position(pos):
        diff = 10000000
        best_cpos = None
        for cpos in xpos_list:
            curr_diff = abs(pos-cpos)
            if curr_diff < diff:
                diff = curr_diff
                best_cpos = cpos
            else:
                break
        return best_cpos

    # get ordering of primary and filtered reads
    primary_reads_l = list(primary_reads)
    primary_reads_l.sort(key=lambda x: read_data[x][START_POS]*300000000+read_data[x][END_POS])
    primary_reads_h1 = list(filter(lambda x: x in truth_h1, primary_reads_l))
    primary_reads_h2 = list(filter(lambda x: x in truth_h2, primary_reads_l))
    filtered_reads_l = list(filtered_reads)
    filtered_reads_l.sort(key=lambda x: read_data[x][START_POS]*300000000+read_data[x][END_POS])
    filtered_reads_h1 = list(filter(lambda x: x in truth_h1 and x not in primary_reads, filtered_reads_l))
    filtered_reads_h2 = list(filter(lambda x: x in truth_h2 and x not in primary_reads, filtered_reads_l))

    # helper for getting non-overlapping read lists
    def generate_non_overlapping_read_list(input_read_names):
        output_read_name_lists = list()
        while len(input_read_names) > 0:
            read_name_list = list()
            output_read_name_lists.append(read_name_list)
            read_name = input_read_names.pop(0)
            read_name_list.append(read_name)
            while True:
                curr_read_name = read_name_list[-1]
                read_end_pos = read_data[curr_read_name][END_POS]
                closest_read_end_pos = closest_position(read_end_pos)
                xpos_idx = xpos_list.index(closest_read_end_pos)
                if xpos_idx + 2 >= len(xpos_list):
                    break
                next_read_start_pos = xpos_list[xpos_idx+2]
                next_read_name = None
                for irn in input_read_names:
                    if read_data[irn][START_POS] >= next_read_start_pos:
                        next_read_name = irn
                        break
                if next_read_name is not None:
                    input_read_names.remove(next_read_name)
                    read_name_list.append(next_read_name)
                else:
                    break
        return output_read_name_lists

    # loggit
    # log_filtered_read_types(read_data, xpos_list, primary_reads, filtered_reads, filtered_bubbles_merged)

    # helper to generate ypos_map
    ypos_map = {}
    def update_ypos_map(read_list_list, ypos):
        for read_list in read_list_list:
            for read in read_list:
                ypos_map[read] = -1 * ypos
            ypos += 1
        return ypos

    # save read positioning
    ypos = update_ypos_map(generate_non_overlapping_read_list(primary_reads_h1), 0)
    ypos_map["spacer_p"] = -1 * ypos
    ypos = update_ypos_map(generate_non_overlapping_read_list(primary_reads_h2), ypos + 1)
    ypos_map["spacer_pf"] = -1 * ypos
    ypos = update_ypos_map(generate_non_overlapping_read_list(filtered_reads_h1), ypos + 1)
    ypos_map["spacer_f"] = -1 * ypos
    update_ypos_map(generate_non_overlapping_read_list(filtered_reads_h2), ypos + 1)

    # sizing config
    SPACER_POINT_SIZE=256
    SPACER_POINT_ALPHA=.5

    # write points to plot
    # plt.figure(num=None, figsize=(.15 * len(set(ypos_map.values())), .1 * len(xpos_list)), dpi=80, facecolor='w', edgecolor='k')
    plt.figure(num=None, figsize=(18,12))
    for read in read_data.values():
        # first iterate over reads (ypos)
        name = read[NAME]
        if name not in ypos_map:
            continue
        y = ypos_map[name]

        # then iterate over sites (xpos)
        firstReadPos = True
        for x,p in enumerate(xpos_list):
            # only plot points where reads overlap
            if read[START_POS] <= p < read[END_POS]:
                # verbose means we log read names
                if firstReadPos and args.verbose:
                    plt.annotate(name if len(name) < 8 else name[:8], (x-1, y-.5), fontfamily='monospace', fontsize=6)
                    firstReadPos = False

                # plot primary reads
                if name in primary_reads:
                    if p in primary_bubbles:
                        color="grey"
                        marker="_"
                        alpha=1.0
                        if name in primary_bubbles[p][READS]:
                            bub_read = primary_bubbles[p][READS][name]
                            supportH1 = bub_read[HAP_SUPPORT_H1]
                            supportH2 = bub_read[HAP_SUPPORT_H2]
                            if supportH1 > supportH2:
                                color = "blue"
                                marker = "o"
                                alpha = -1 / supportH1
                            elif supportH2 > supportH1:
                                color = "red"
                                marker = "o"
                                alpha = -1 / supportH2
                        plt.scatter(x,y,marker=marker,color=color,alpha=alpha,edgecolors="black")
                    else:
                        color = "darkblue" if read[HAP] == 1 else ("darkred" if read[HAP] == 2 else "dimgrey")
                        plt.scatter(x,y, marker="4" if read[STRAND] == "+" else "3", s=8, color=color, alpha=SPACER_POINT_ALPHA)

                # plot filtered reads
                elif name in filtered_reads:
                    if p in filtered_bubbles:
                    # if p in filtered_bubbles_merged:
                        color="grey"
                        marker="_"
                        alpha=1.0
                        if name in filtered_bubbles[p][READS]:
                            bub_read = filtered_bubbles[p][READS][name]
                            supportH1 = bub_read[HAP_SUPPORT_H1]
                            supportH2 = bub_read[HAP_SUPPORT_H2]
                            if supportH1 > supportH2:
                                color = "blue"
                                marker = "o"
                                alpha = -1 / supportH1
                            elif supportH2 > supportH1:
                                color = "red"
                                marker = "o"
                                alpha = -1 / supportH2
                        plt.scatter(x,y,marker=marker,color=color,alpha=alpha,edgecolors="black")
                    else:
                        color = "darkblue" if read[HAP] == 1 else ("darkred" if read[HAP] == 2 else "dimgrey")
                        plt.scatter(x,y, marker="4" if read[STRAND] == "+" else "3", s=8, color=color, alpha=SPACER_POINT_ALPHA)

    # separate primary, filtered, and haplotypes
    plt.hlines(y=ypos_map['spacer_p'], color="black", linestyle=':', xmin=-1, xmax=len(xpos_list), alpha=.5)
    plt.hlines(y=ypos_map['spacer_f'], color="black", linestyle=':', xmin=-1, xmax=len(xpos_list), alpha=.5)
    plt.hlines(y=ypos_map['spacer_pf'], color="black", linestyle='--', xmin=-1, xmax=len(xpos_list), alpha=.5)

    # describe het sites
    plt.xticks(xpos_ticks, xpos_labels, rotation=45, fontsize=10)
    plt.yticks([],[])
    plt.tight_layout()
    if args.figure_name is not None:
        plt.savefig(args.figure_name, dpi=360)
    plt.show()
    plt.close()






if __name__ == "__main__":
    main()

