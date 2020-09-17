#!/usr/bin/env python3
import argparse
import sys
import json
import glob
import matplotlib.pyplot as plt
from matplotlib import cm
import math
from collections import defaultdict
import pysam
import collections
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.stats import pearsonr


PRIMARY="primary"
REF_POS="refPos"

# HET_DIST_RANGES = [0, 1000, 2000, 4000, 8000, 12000, 16000, 20000, 26000, 32000]
HET_DIST_RANGES = [2000, 4000, 8000, 12000, 16000, 24000, 32000]
# HET_DIST_RANGES = [2000, 4000, 8000, 16000, 24000]


class Chunk:
    def __init__(self, idx, contig, overlap_start, overlap_end, chunk_start, chunk_end):
        self.idx = idx
        self.contig = contig
        self.overlap_start = overlap_start
        self.overlap_end = overlap_end
        self.chunk_start = chunk_start
        self.chunk_end = chunk_end
        self.het_sites = list()
        self.boundary_sites = list()
        self.phasing_file = None

    def add_het_site(self, site_pos):
        if self.chunk_start <= site_pos < self.chunk_end:
            self.het_sites.append(site_pos)
            return True
        elif self.overlap_start <= site_pos <= self.overlap_end:
            self.boundary_sites.append(site_pos)
            return False
        else:
            log("Site {} in chunk {} outside of boundaries {}, {}".format(
                site_pos, self.idx, self.overlap_start, self.overlap_end))
            return None


class ChunkData:
    def __init__(self):
        self.chunk_by_idx = dict()
        self.contig_chunk_start = dict()
        self.contig_chunk_end = dict()
        self.contigs = set()

    def save_chunk(self, idx, contig, overlap_start, overlap_end, chunk_start, chunk_end):
        chunk = Chunk(idx, contig, overlap_start, overlap_end, chunk_start,chunk_end)
        self.chunk_by_idx[idx] = chunk
        if contig in self.contigs:
            if chunk_start < self.contig_chunk_start[contig].chunk_start:
                self.contig_chunk_start[contig] = chunk
            if chunk_end > self.contig_chunk_end[contig].chunk_end:
                self.contig_chunk_end[contig] = chunk
        else:
            self.contigs.add(contig)
            self.contig_chunk_start[contig] = chunk
            self.contig_chunk_end[contig] = chunk


class RegionOfInterest:
    def __init__(self, contig, start, end, bucket_size):
        self.contig = contig
        self.start = start
        self.end = end
        self.chunks = set()
        self.hets_by_distance = dict()

        self.bucket_size = bucket_size
        self.bucket_count = int(math.ceil((end - start) / bucket_size))
        for dist in HET_DIST_RANGES:
            self.hets_by_distance[dist] = [0 for _ in range(self.bucket_count)]

    def desc(self):
        return "{}-{}-{}".format(self.contig, self.start, self.end)

    def bucket_idx(self, pos):
        bucket = (pos - self.start) // self.bucket_size
        if bucket < 0 or bucket >= self.bucket_count:
            return None
        return bucket

    def handle_het_site(self, pos):
        distance_from_pos = 0
        while distance_from_pos < max(HET_DIST_RANGES):
            bucket = self.bucket_idx(pos + distance_from_pos)
            if bucket is not None:
                for dist in HET_DIST_RANGES:
                    if distance_from_pos > dist: continue
                    self.hets_by_distance[dist][bucket] += 1
            bucket = self.bucket_idx(pos - distance_from_pos)
            if distance_from_pos != 0 and bucket is not None:
                for dist in HET_DIST_RANGES:
                    if distance_from_pos > dist: continue
                    self.hets_by_distance[dist][bucket] += 1
            distance_from_pos += self.bucket_size





def parse_args(args = None):
    parser = argparse.ArgumentParser("Finds homozygous regions where phasing is uncertain")
    parser.add_argument('--phasing_file_glob', '-i', dest='phasing_file_glob', default=None, required=True, type=str,
                        help='Glob matching describing bubbles and phasing JSON files')
    parser.add_argument('--chunking_file', '-c', dest='chunking_file', default=None, required=True, type=str,
                        help='File describing chunk boundaries')
    parser.add_argument('--region', '-r', dest='region', default=None, required=False, type=str,
                        help='Comma-separated regions of interest, in "chr" or "chr:start-end" format')
    parser.add_argument('--bucket_size', '-s', dest='bucket_size', default=256, required=False, type=int,
                        help='Size in BP for heterozyous buckets')
    parser.add_argument('--plot', '-p', dest='plot', default=False, required=False, action='store_true',
                        help='Produce a plot of the data')
    parser.add_argument('--figure_name', '-f', dest='figure_name', default=None, required=False, type=str,
                        help='Figure name (will save only if set), requires --plot flag')

    return parser.parse_args() if args is None else parser.parse_args(args)


def log(msg):
    print(msg, file=sys.stderr)


def has_overlap(s1, e1, s2, e2):
    return s2 <= s1 < e2 or s1 <= s2 < e1


def get_chunk_idx_from_filename(filename):
    vf_parts = filename.split(".")
    assert(len(vf_parts) > 5)
    chunk_part = vf_parts[-4]
    assert(chunk_part[0] == "C")
    return int(chunk_part[1:])


def get_chunk_info(chunk_file):
    log("Reading chunk info from {}".format(chunk_file))
    chunk_data = ChunkData()
    with open(chunk_file) as cin:
        cidx = 0
        for line in cin:
            parts = line.split(",")
            assert(len(parts) == 5)
            chunk_data.save_chunk(cidx, parts[0], int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]))
            cidx += 1
    return chunk_data


def save_basic_variant_info(phasing_file, chunk_data :ChunkData):
    # get chunk data
    cidx = get_chunk_idx_from_filename(phasing_file)
    if cidx not in chunk_data.chunk_by_idx:
        log("\nChunk idx {} for file {} not in chunk data!  Existing keys: \n\t{}".format(
            cidx, phasing_file, chunk_data.chunk_by_idx.keys()))
        sys.exit(1)
    chunk = chunk_data.chunk_by_idx[cidx]
    chunk.phasing_file = phasing_file

    # get json file
    log("  Parsing read and phasing info for chunk {} from {}".format(cidx, phasing_file))
    with open(phasing_file, 'r') as istream:
        json_doc = json.load(istream)

    chunk_het_count = 0
    boundary_het_count = 0
    for bubble in json_doc[PRIMARY]:
        pos = int(bubble[REF_POS])
        pos_in_chunk = chunk.add_het_site(pos)
        if pos_in_chunk:
            chunk_het_count += 1
        else:
            boundary_het_count += 1

    log("    Got {} HETs in chunk and {} out".format(chunk_het_count, boundary_het_count))


def handle_region_of_interest(roi :RegionOfInterest, chunk_data :ChunkData, args):
    # handle each het pos
    for chunk_idx in roi.chunks:
        for site in chunk_data.chunk_by_idx[chunk_idx].het_sites:
            roi.handle_het_site(site)

    # plottit
    plt.figure(num=None, figsize=(18,4))
    colormap = cm.get_cmap('viridis', max(HET_DIST_RANGES))
    x = [(roi.start + b*roi.bucket_size) // 1000 for b in range(roi.bucket_count)]

    # linear
    for dist in reversed(HET_DIST_RANGES):
        plt.plot(x, roi.hets_by_distance[dist], color=colormap(dist))
    plt.title(roi.desc())
    if args.figure_name is not None:
        plt.savefig("{}.{}.linear.png".format(args.figure_name, roi.desc()))
    plt.show()
    plt.close()


def main():
    args = parse_args()

    # get chunks
    chunk_data = get_chunk_info(args.chunking_file)

    # get region
    regions_of_interest = defaultdict(lambda: list())
    chunks_of_interest = set()
    if args.region is not None:
        regions = args.region.split(",")
        for region in regions:
            region_parts = region.split(":")
            ctg = region_parts[0]
            start = 0 if len(region_parts) == 1 else int(region_parts[1].split("-")[0])
            end = sys.maxsize if len(region_parts) == 1 else int(region_parts[1].split("-")[-1])
            regions_of_interest[ctg].append(RegionOfInterest(ctg, start, end, args.bucket_size))
        for chunk in chunk_data.chunk_by_idx.values():
            if chunk.contig not in regions_of_interest: continue
            for roi in regions_of_interest[chunk.contig]:
                if has_overlap(chunk.chunk_start, chunk.chunk_end, roi.start, roi.end):
                    chunks_of_interest.add(chunk.idx)
                    roi.chunks.add(chunk.idx)
    else:
        for contig in chunk_data.contigs:
            regions_of_interest[contig].append(RegionOfInterest(
                contig, chunk_data.contig_chunk_start[contig], chunk_data.contig_chunk_end[contig], args.bucket_size))
        chunks_of_interest.intersection_update(chunk_data.chunk_by_idx.values())


    # get variant info
    phasing_files = glob.glob(args.phasing_file_glob)
    log("Reading {} phasing files matching {}".format(len(phasing_files), args.phasing_file_glob))
    for phasing_file in phasing_files:
        chunk_idx = get_chunk_idx_from_filename(phasing_file)
        if chunk_idx not in chunks_of_interest:
            continue
        save_basic_variant_info(phasing_file, chunk_data)

    # look at contigs we care about
    all_regions_of_interest = list()
    for contig_regions in regions_of_interest.values():
        for roi in contig_regions:
            all_regions_of_interest.append(roi)
    log("Examining {} regions of interest".format(len(regions_of_interest)))
    for roi in all_regions_of_interest:
        log("Examining region {}".format(roi.desc()))
        handle_region_of_interest(roi, chunk_data, args)







if __name__ == "__main__":
    main()

