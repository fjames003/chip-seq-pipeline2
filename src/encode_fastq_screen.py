#!/usr/bin/env python

# ENCODE DCC fastq_screen wrapper
# Author: Frankie James (fjames003)

import sys
import os
import re
import argparse
import multiprocessing
from encode_common_genomic import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC Fastq Screen.',
                                        description='')
    parser.add_argument('fastqs', nargs='+', type=str,
                        help='List of FASTQs (R1 and R2). \
                            FASTQs must be compressed with gzip (with .gz).')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end FASTQs.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--conf', type=str,
                        help='Configuration file for Fastq Screen.')
    parser.add_argument('--out-dir', default='', type=str,
                            help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    # check if fastqs have correct dimension
    if args.paired_end and len(args.fastqs)!=2:
        raise argparse.ArgumentTypeError('Need 2 fastqs for paired end.')
    if not args.paired_end and len(args.fastqs)!=1:
        raise argparse.ArgumentTypeError('Need 1 fastq for single end.')

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def screen_fastqs(fastqs, nth, conf, out_dir):
    basename1 = os.path.basename(strip_ext_fastq(fastqs[0]))
    prefix1 = os.path.join(out_dir, basename1)

    cmd =  'fastq_screen --threads {} '
    cmd += '--conf {} '
    cmd += '--aligner bowtie2 '
    cmd += '--outdir {} '
    cmd += '{}'
    cmd = cmd.format(nth, conf, " ".join(fastqs))

    run_shell_cmd(cmd)

    return ""

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # declare temp arrays
    temp_files = [] # files to deleted later at the end

    # STAR
    log.info('Running Fastq Screen...')
    temp = screen_fastqs(args.fastqs, args.nth, args.conf, args.out_dir)


    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()
