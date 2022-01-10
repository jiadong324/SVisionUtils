#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/1/10
'''

import pysam
import os
import pandas as pd
from intervaltree import IntervalTree

from Configs import *


def run_filter_region_graph(svision_vcf, svision_exact_graph, ref_fasta, exclude_file, exclude_graphid, min_sr, max_sv_size, outdir):

    filtered_prefix = '.'.join(os.path.basename(svision_vcf).split('.')[0:-1])

    exclude_dict = parse_exclude_regions(exclude_file)
    ref_file = pysam.FastaFile(ref_fasta)
    incomplete_graphs, complete_graphs = get_incomplete_graphID(svision_exact_graph)

    csv_by_id = {}
    all_sv_num = 0
    filtered_vcf_writer = open(outdir + f'/{filtered_prefix}.filtered.vcf', 'w')

    with open(svision_vcf, 'r') as f:
        for line in f:
            if "#" in line:
                continue
            entries = line.strip().split("\t")
            chrom = entries[0]
            id = entries[2]
            if "_" in id:
                id = id.split("_")[0]

            if chrom in AUTOSOMES:
                start = int(entries[1])

                info_tokens = entries[7].split(";")
                info_dict = {}

                for token in info_tokens:
                    info_dict[token.split("=")[0]] = token.split("=")[1]

                this_chrom_exclude_tree = exclude_dict[chrom]
                if this_chrom_exclude_tree.overlaps(start, int(info_dict['END'])):
                    continue

                if contains_gaps(chrom, start, int(info_dict['END']), ref_file):
                    continue

                if int(info_dict["SVLEN"]) < max_sv_size and int(info_dict['SUPPORT']) >= min_sr:
                    all_sv_num += 1
                    if info_dict['GraphID'] != "-1":
                        ## Modified v1.3.5
                        if info_dict['GraphID'] in exclude_graphid:
                            continue

                        if info_dict['GraphID'] not in incomplete_graphs and entries[6] != 'Uncovered':
                            if id in csv_by_id:
                                csv_by_id[id].append((chrom, start, info_dict['END'], entries[2], info_dict['GraphID'], int(info_dict['SUPPORT'])))
                            else:
                                csv_by_id[id] = [(chrom, start, info_dict['END'], entries[2], info_dict['GraphID'], int(info_dict['SUPPORT']))]

                            filtered_vcf_writer.write(line)
                    else:
                        filtered_vcf_writer.write(line)

    high_conf_csv_list = list()
    for id, csvs in csv_by_id.items():

        sorted_csvs = sorted(csvs, key=lambda x: x[-1], reverse=True)
        selected_csv = sorted_csvs[0]
        high_conf_csv_list.append((selected_csv[0], int(selected_csv[1]), int(selected_csv[2]), selected_csv[3], selected_csv[4]))

    df_high_conf_csvs = pd.DataFrame(high_conf_csv_list, columns=['chrom', 'start', 'end', 'id', 'graphid'])
    sorter_index = dict(zip(AUTOSOMES, range(len(AUTOSOMES))))
    df_high_conf_csvs['chrom_rank'] = df_high_conf_csvs['chrom'].map(sorter_index)

    # Write high confident csvs
    df_high_conf_csvs.sort_values(by=['chrom_rank', 'start'], inplace=True)
    df_high_conf_csvs.drop('chrom_rank', 1, inplace=True)
    df_high_conf_csvs.to_csv(outdir + f"/{filtered_prefix}.Raw-CSVs.tsv", sep="\t", header=True, index=False)

    print('SV after filtering {0}, containing {1} CSVs'.format(all_sv_num, len(high_conf_csv_list)))


def parse_exclude_regions(exclude):
    exclude_dict = {}
    with open(exclude, 'r') as f:
        for line in f:
            entries = line.strip().split("\t")
            start = int(entries[1])
            end = int(entries[2])
            if entries[0] not in exclude_dict:
                exclude_dict[entries[0]] = IntervalTree()
                exclude_dict[entries[0]][start:end] = (start, end)
            else:
                exclude_dict[entries[0]][start:end] = (start, end)
    return exclude_dict

def contains_gaps(chrom, start, end, ref):
    seq = ref.fetch(chrom, start - 2000, end + 2000)
    in_gap = False
    for base in seq:
        if base == 'N':
            in_gap = True
            break
    return in_gap

def get_incomplete_graphID(graph):
    incomplete_graphs = []
    complete_graphs = []
    with open(graph, 'r') as f:
        for line in f:
            if ">" in line:
                entries = line.strip().split("\t")
                id = entries[0].split("=")[1]
                path = entries[4].split("=")[1]

                if path[-3] != "S":
                    incomplete_graphs.append(id)
                else:
                    complete_graphs.append(id)

    return incomplete_graphs, complete_graphs

def overlap(a,b,c,d):
    r = 0 if a==c and b==d else min(b,d)-max(a,c)
    if r>=0: return r