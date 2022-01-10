#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/1/10
'''

import os
import pandas as pd


from Configs import *

def annot_csvs_repeats(high_conf_csvs, outdir):

    prefix = os.path.basename(high_conf_csvs).split('.')[0]
    rmsk_overlap_file = os.path.join(outdir, prefix + '.rmsk.overlap.tsv')
    rmsk_overlap = '{0} intersect -wa -wb -a {1} -b {2} > {3}'.format(BEDTOOLS, high_conf_csvs, RMSK, rmsk_overlap_file)
    os.system(rmsk_overlap)

    trf_overlap_file = os.path.join(outdir, prefix + '.trf.overlap.tsv')
    trf_overlap = '{0} intersect -wa -wb -a {1} -b {2} > {3}'.format(BEDTOOLS, high_conf_csvs, TRF, trf_overlap_file)
    os.system(trf_overlap)


def assign_csv_repeats(workdir, csvs):

    prefix = os.path.basename(csvs).split('.')[0]
    rmsk_overlaps = os.path.join(workdir, prefix + '.rmsk.overlap.tsv')
    trf_overlaps = os.path.join(workdir, prefix + '.trf.overlap.tsv')
    annot_out = os.path.join(workdir, prefix + '.repannot.tsv')

    df_rmsk = pd.read_csv(rmsk_overlaps, sep="\t", names=['chrom', 'start', 'end', 'id', 'graphid', 'rpchrom', 'rpstart', 'rpend', 'rpsubtype', 'rptype'])

    annot_by_csv = {}

    for idx, row in df_rmsk.iterrows():
        if row['rptype'] == 'Simple_repeat':
            continue
        csv_id = "{0}-{1}-{2}-{3}".format(row['chrom'], row['start'], row['end'], row['id'])
        overlap_size = overlap(int(row['start']), int(row['end']), int(row['rpstart']), int(row['rpend']))
        if csv_id in annot_by_csv:
            annot_by_csv[csv_id].append((row['rpsubtype'], row['rptype'], overlap_size))
        else:
            annot_by_csv[csv_id] = [(row['rpsubtype'], row['rptype'], overlap_size)]


    df_trf = pd.read_csv(trf_overlaps, sep="\t", names=['chrom', 'start', 'end', 'id', 'graphid', 'rpchrom', 'rpstart', 'rpend', 'motif'])

    for idx, row in df_trf.iterrows():
        csv_id = "{0}-{1}-{2}-{3}".format(row['chrom'], row['start'], row['end'], row['id'])
        overlap_size = overlap(int(row['start']), int(row['end']), int(row['rpstart']), int(row['rpend']))
        subtype = "STR"
        if len(row['motif']) >= 7:
            subtype = 'VNTR'

        if csv_id not in annot_by_csv:
            annot_by_csv[csv_id] = [(subtype, subtype, overlap_size)]
        else:
            annot_by_csv[csv_id].append((subtype, subtype, overlap_size))

    df_csvs = pd.read_csv(csvs, sep="\t", names=['chrom', 'start', 'end', 'id', 'graphid'], header=None, skiprows=1)

    csv_annots = list()
    count_by_rptype = {}
    for idx, row in df_csvs.iterrows():
        csv_id = "{0}-{1}-{2}-{3}".format(row['chrom'], row['start'], row['end'], row['id'])
        rptype = 'None'
        rep_overlaps = 0
        if csv_id in annot_by_csv:

            annots = annot_by_csv[csv_id]
            if len(annots) == 1:
                rptype = annots[0][1]
                rep_overlaps = annots[0][2]
            else:
                sorted_annots_by_size = sorted(annots, key=lambda x:x[2], reverse=True)
                rptype = sorted_annots_by_size[0][1]
                rep_overlaps = sorted_annots_by_size[0][2]

        msk_pcrt = 100 * rep_overlaps / (int(row['end']) - int(row['start']))

        csv_annots.append((row['chrom'], row['start'], row['end'], row['id'], row['graphid'], round(msk_pcrt, 2), rptype))

        if rptype in count_by_rptype:
            count_by_rptype[rptype] += 1
        else:
            count_by_rptype[rptype] = 1


    df_csv_annots = pd.DataFrame(csv_annots, columns=['chrom', 'start', 'end', 'id', 'graphid', 'pcrt', 'rptype'])

    sorter_index = dict(zip(AUTOSOMES, range(len(AUTOSOMES))))
    df_csv_annots['chrom_rank'] = df_csv_annots['chrom'].map(sorter_index)
    df_csv_annots.drop('chrom_rank', 1, inplace=True)
    df_csv_annots.to_csv(annot_out, index=False, header=False, sep="\t")


def overlap(a,b,c,d):
    r = 0 if a==c and b==d else min(b,d)-max(a,c)
    if r>=0: return r


def filter_csvs_by_trs(rep_annot, output):

    output_writer = open(output, 'w')
    simple_repeats = ['VNTR', 'STR']

    tr_calls = 0
    df_csv_reps = pd.read_csv(rep_annot, sep='\t', names=['chrom', 'start', 'end', 'id', 'graphid', 'pcrt', 'reptype'])

    for idx, row in df_csv_reps.iterrows():
        if row['reptype'] in simple_repeats:
            tr_calls += 1
            continue

        output_str = '\t'.join([str(val) for val in row.tolist()])
        print(output_str, file=output_writer)

    output_writer.close()
    print('CSVs outside TRs: ', len(df_csv_reps) - tr_calls)


def run_filter_rep(workdir, svision_vcf):

    filtered_prefix = '.'.join(os.path.basename(svision_vcf).split('.')[0:-1])

    raw_csv = f'{workdir}/{filtered_prefix}.Raw-CSVs.tsv'
    ## Step 1 repeat annotation
    annot_csvs_repeats(raw_csv, workdir)

    ## Step 2 assign repeat element of each CSV
    assign_csv_repeats(workdir, raw_csv)

    ## Step 3 filter CSV by VNTR/STR
    prefix = os.path.basename(raw_csv).split('.')[0]
    filter_csvs_by_trs(f'{workdir}/{prefix}.repannot.tsv', f'{workdir}/{prefix}.HQ-CSVs.tsv')