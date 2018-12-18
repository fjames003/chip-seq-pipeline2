#!/usr/bin/env python

# Spike in Calibration wrapper
# Author: Frankie James (fjames003)

import sys
import os
import re
import argparse
import multiprocessing
from encode_common_genomic import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC STAR aligner.',
                                        description='')
    parser.add_argument('genome_bam', type=str,
                        help='Path to bam file aligned to primary genome')
    parser.add_argument('spike_in_bam', type=str,
                        help='Path to bam file aligned to spike in genome')
    parser.add_argument('chrom_sizes', type=str,
                        help='Path to file containing chromosome sizes')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--out-dir', default='', type=str,
                            help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    # check if fastqs have correct dimension
    if not os.path.exists(args.genome_bam):
        raise argparse.ArgumentTypeError('Need a bam file aligned to primary genome.')
    if not os.path.exists(args.spike_in_bam):
        raise argparse.ArgumentTypeError('Need a bam file aligned to spike in genome.')
    if not os.path.exists(args.chrom_sizes):
        raise argparse.ArgumentTypeError('Need to provide file with chrom sizes.')

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def generate_beds(genome_bam, spike_in_bam, out_dir):
    g_basename = "Aligned_to_primary_genome"
    g_prefix = os.path.join(out_dir, g_basename)

    s_basename = "Aligned_to_spikeIn_genome"
    s_prefix = os.path.join(out_dir, s_basename)

    genome_bed = "{}.bed".format(g_prefix)
    spike_in_bed = "{}.bed".format(s_prefix)

    cmd =  "bedtools bamtobed -i {0} | "
    cmd += "awk -v OFS=\"\\t\" "
    cmd += "'{{len = $3 - $2; print $1, $2, $3, len }}' > {1}"
    cmd = cmd.format(genome_bam, genome_bed)
    run_shell_cmd(cmd)

    cmd1 = "bedtools bamtobed -i {0} > {1}"
    cmd1 = cmd1.format(spike_in_bam, spike_in_bed)
    run_shell_cmd(cmd1)

    return genome_bed, spike_in_bed

def run_calibration(genome_bed, spike_bed, chrom_sizes, out_dir):
    basename = os.path.basename(strip_ext_bed(genome_bed))
    prefix = os.path.join(out_dir, basename)

    bed_graph = "{}.bedgraph"

    cmd = "bash spike_in_calibration.sh {0} {1} 1000 bg {2} 1 1000"
    cmd = cmd.format(genome_bed, spike_bed, chrom_sizes)

    run_shell_cmd(cmd)

    return bed_graph

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    genome_bed, spike_in_bed = generate_beds(args.genome_bam, args.spike_in_bam, args.out_dir)

    bed_graph = run_calibration(genome_bed, spike_in_bed, args.chrom_sizes, args.out_dir)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()
