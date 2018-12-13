#!/usr/bin/env python

# ENCODE DCC Cutadapt wrapper
# Author: Frankie James (fjames003)

import sys
import os
import re
import argparse
import multiprocessing
from encode_common_genomic import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC Cutadapt adapter trimmer.',
                                        description='')
    parser.add_argument('fastqs', nargs='+', type=str,
                        help='List of FASTQs (R1 and R2). \
                            FASTQs must be compressed with gzip (with .gz).')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end FASTQs.')
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
    if args.paired_end and len(args.fastqs)!=2:
        raise argparse.ArgumentTypeError('Need 2 fastqs for paired end.')
    if not args.paired_end and len(args.fastqs)!=1:
        raise argparse.ArgumentTypeError('Need 1 fastq for single end.')

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def cutadapt_se(fastq, nth, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq))
    prefix = os.path.join(out_dir, basename)
    fastq_out = "{}.trimmed.fastq.gz".format(prefix)

    cmd =  'cutadapt -j {} '
    cmd += '--length 36 --minimum-length 35 '
    cmd += '--trim-n -e 0.1 -q 30 '
    cmd += '-a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT '
    cmd += '-o {} '
    cmd += '<(zcat {})'
    cmd = cmd.format(nth, fastq_out, fastq1)

    run_shell_cmd(cmd)

    return fastq_out

def cutadapt_pe(fastq1, fastq2, nth, out_dir):
    basename1 = os.path.basename(strip_ext_fastq(fastq1))
    prefix1 = os.path.join(out_dir, basename1)
    fastq_out1 = "{}.trimmed.fastq.gz".format(prefix1)

    basename2 = os.path.basename(strip_ext_fastq(fastq2))
    prefix2 = os.path.join(out_dir, basename2)
    fastq_out2 = "{}.trimmed.fastq.gz".format(prefix2)

    cmd =  'cutadapt -j {} '
    cmd += '--length 36 --minimum-length 35 '
    cmd += '--trim-n -e 0.1 -q 30 '
    cmd += '-a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT '
    cmd += '-o {} -p {} '
    cmd += '<(zcat {}) '
    cmd += '<(zcat {})'
    cmd = cmd.format(nth, fastq_out1, fastq_out2, fastq1, fastq2)

    run_shell_cmd(cmd)
    
    return fastq_out1, fastq_out2

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # declare temp arrays
    temp_files = [] # files to deleted later at the end

    # STAR
    log.info('Running Cutadapt...')
    if args.paired_end:
        fastq_out1, fastq_out2 = cutadapt_pe(args.fastqs[0], args.fastqs[1], args.nth, args.out_dir)
    else:
        fastq_out1 = cutadapt_se(args.fastqs[0], args.nth, args.out_dir)


    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()
