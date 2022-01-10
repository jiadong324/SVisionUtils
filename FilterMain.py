#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/1/10
'''

from FilterRG import run_filter_region_graph
from FilterRep import run_filter_rep


def main():

    ## Filtering for the results in the paper
    svision_vcf = '/Users/jiadonglin/SVision/HG00733/svision/v136/HG00733_docker_test.svision.s5.graph.vcf'

    svision_exact_graph = '/Users/jiadonglin/SVision/HG00733/svision/v136/HG00733_docker_test.graph_exactly_match.txt'
    ref_fasta = "/Users/jiadonglin/Data/genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    exclude_file = "/Users/jiadonglin/Data/genome/grch38.exclude_regions_cen.bed"
    outdir = '/Users/jiadonglin/SVision/HG00733/svision/v136'

    exclude_graphid = ['0', '3']

    # svision_vcf = options.vcf
    # svision_exact_graph = options.graph
    # ref_fasta = options.ref
    # exclude_file = options.region
    # outdir = options.outdir
    # exclude_graphid = options.id.split(',')
    # min_sr = options.min_sr
    # max_sv_size = options.max_sv_size

    print('===== Step 1 filter CSVs by region and graph structures =====')
    run_filter_region_graph(svision_vcf, svision_exact_graph, ref_fasta, exclude_file, exclude_graphid, 5, 100000, outdir)

    print('===== Step2 filter CSVs by simple repeats =====')
    run_filter_rep(outdir, svision_vcf)


if __name__ == '__main__':
    main()