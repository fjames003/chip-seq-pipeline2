#!/usr/bin/env python

# ENCODE DCC STAR wrapper
# Author: Jin Lee (leepc12@gmail.com), Daniel Kim

import sys
import os
import re
import argparse
import multiprocessing
from encode_common_genomic import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC STAR aligner.',
                                        description='')
    parser.add_argument('STAR_index_prefix_or_tar', type=str,
                        help='Path for prefix (or a tarball .tar) \
                            for reference STAR index. \
                            Prefix must be like [PREFIX].sa. \
                            Tar ball must be packed without compression \
                            and directory by using command line \
                            "tar cvf [TAR] [TAR_PREFIX].*".')
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

# WDL glob() globs in an alphabetical order
# so R1 and R2 can be switched, which results in an
# unexpected behavior of a workflow.
# so we already prepended merge_fastqs_'end'_ (R1 or R2)
# to the basename of original filename in 'trim_adapter' task.
# now it's time to strip it.
def strip_merge_fastqs_prefix(fastq):
    return re.sub(r'^merge\_fastqs\_R\d\_','',str(fastq))

def make_read_length_file(fastq, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq))
    prefix = os.path.join(out_dir,
        strip_merge_fastqs_prefix(basename))
    txt = '{}.read_length.txt'.format(prefix)
    read_length = get_read_length(fastq)
    with open(txt,'w') as fp:
        fp.write(str(read_length))
    return txt


def STAR_se(fastq, nth, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq))
    prefix = os.path.join(out_dir,
        strip_merge_fastqs_prefix(basename))
    bam = '{}Aligned.sortedByCoord.out.bam'.format(prefix)

    cmd = 'STAR '
    cmd += '--runThreadN {0} '
    cmd += '--genomeDir /home/segil_lab/genomes/gencode_mm10_genomeDir '
    cmd += '--readFilesCommand zcat '
    cmd += '--readFilesIn {1} '
    cmd += '--outFileNamePrefix {2} '
    cmd += '--outSAMtype BAM SortedByCoordinate '
    cmd += '--alignEndsType EndToEnd '
    cmd += '--alignIntronMax 1 '
    cmd += '--outFilterMultimapNmax 1 '
    cmd += '--outSAMmultNmax 1 '
    cmd += '--outMultimapperOrder Random '
    cmd += '--outFilterMismatchNmax 5 '
    cmd += '--outFilterMismatchNoverLmax 0.3 '
    cmd += '--outFilterMatchNminOverLread 0.4 '
    cmd += '--outFilterScoreMinOverLread 0 '
    cmd += '--limitBAMsortRAM 26843545600 '
    cmd += '--genomeLoad LoadAndKeep '
    cmd += '--outBAMsortingThreadN {0} '
    cmd += '--alignMatesGapMax 0 '
    cmd += '--outSAMunmapped Within '
    cmd += '--seedSplitMin 10'
    cmd = cmd.format(nth, fastq, prefix)
    run_shell_cmd(cmd)

    return bam

def STAR_pe(fastq1, fastq2, nth, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq1))
    prefix = os.path.join(out_dir,
        strip_merge_fastqs_prefix(basename))
    # sam = '{}.sam'.format(prefix)
    # badcigar = '{}.badReads'.format(prefix)
    bam = '{}Aligned.sortedByCoord.out.bam'.format(prefix)

    cmd = 'STAR '
    cmd += '--runThreadN {} '
    cmd += '--genomeDir /home/segil_lab/genomes/gencode_mm10_genomeDir '
    cmd += '--readFilesCommand zcat '
    cmd += '--readFilesIn {} {} '
    cmd += '--outFileNamePrefix {} '
    cmd += '--outSAMtype BAM SortedByCoordinate '
    cmd += '--alignEndsType EndToEnd '
    cmd += '--alignIntronMax 1 '
    cmd += '--outFilterMultimapNmax 1 '
    cmd += '--outSAMmultNmax 1 '
    cmd += '--outMultimapperOrder Random '
    cmd += '--outFilterMismatchNmax 5 '
    cmd += '--outFilterMismatchNoverLmax 0.3 '
    cmd += '--outFilterMatchNminOverLread 0.4 '
    cmd += '--outFilterScoreMinOverLread 0 '
    cmd += '--limitBAMsortRAM 26843545600 '
    cmd += '--genomeLoad LoadAndKeep '
    cmd += '--outBAMsortingThreadN 10 '
    cmd += '--alignMatesGapMax 0 '
    cmd += '--outSAMunmapped Within '
    cmd += '--seedSplitMin 10'
    cmd = cmd.format(nth, fastq1, fastq2, prefix)
    # # multiprocessing for STAR_aln
    # pool = multiprocessing.Pool(2)
    # ret_val1 = pool.apply_async(
    #     STAR_aln, (fastq1, ref_index_prefix, max(1,nth/2), out_dir))
    # ret_val2 = pool.apply_async(
    #     STAR_aln, (fastq2, ref_index_prefix, max(1,nth/2), out_dir))
    # sai1 = ret_val1.get(BIG_INT)
    # sai2 = ret_val2.get(BIG_INT)
    # pool.close()
    # pool.join()

    # cmd = 'STAR sampe {} {} {} {} {} | gzip -nc > {}'.format(
    #     ref_index_prefix,
    #     sai1,
    #     sai2,
    #     fastq1,
    #     fastq2,
    #     sam)
    # run_shell_cmd(cmd)
    #
    # cmd2 = 'zcat -f {} | awk \'BEGIN {{FS="\\t" ; OFS="\\t"}} ! /^@/ && $6!="*" '
    # cmd2 += '{{ cigar=$6; gsub("[0-9]+D","",cigar); n = split(cigar,vals,"[A-Z]"); s = 0; '
    # cmd2 += 'for (i=1;i<=n;i++) s=s+vals[i]; seqlen=length($10); '
    # cmd2 += 'if (s!=seqlen) print $1"\\t"; }}\' | '
    # cmd2 += 'sort | uniq > {}'
    # cmd2 = cmd2.format(
    #     sam,
    #     badcigar)
    # run_shell_cmd(cmd2)
    #
    # # Remove bad CIGAR read pairs
    # if get_num_lines(badcigar)>0:
    #     cmd3 = 'zcat -f {} | grep -v -F -f {} | '
    #     cmd3 += 'samtools view -Su - | samtools sort - {}'
    #     cmd3 = cmd3.format(
    #         sam,
    #         badcigar,
    #         prefix)
    # else:
    #     cmd3 = 'samtools view -Su {} | samtools sort - {}'
    #     cmd3 = cmd3.format(
    #         sam,
    #         prefix)
    # run_shell_cmd(cmd3)
    #
    # rm_f([sai1, sai2, sam])
    return bam

def chk_STAR_index(prefix):
    index_sa = '{}.sa'.format(prefix)
    if not os.path.exists(index_sa):
        raise Exception("STAR index does not exists. "+
            "Prefix = {}".format(prefix))

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # declare temp arrays
    temp_files = [] # files to deleted later at the end

    # generate read length file
    log.info('Generating read length file...')
    R1_read_length_file = make_read_length_file(
                            args.fastqs[0], args.out_dir)
    if args.paired_end:
        R2_read_length_file = make_read_length_file(
                            args.fastqs[1], args.out_dir)

    # STAR
    log.info('Running STAR...')
    if args.paired_end:
        bam = STAR_pe(args.fastqs[0], args.fastqs[1], args.nth, args.out_dir)
    else:
        bam = STAR_se(args.fastqs[0], args.nth, args.out_dir)

    # initialize multithreading
    log.info('Initializing multi-threading...')
    num_process = min(2,args.nth)
    log.info('Number of threads={}.'.format(num_process))
    pool = multiprocessing.Pool(num_process)

    log.info('Running samtools index...')
    ret_val1 = pool.apply_async(
        samtools_index, (bam, args.out_dir))

    log.info('Running samtools flagstat...')
    ret_val2 = pool.apply_async(
        samtools_flagstat, (bam, args.out_dir))

    bai = ret_val1.get(BIG_INT)
    flagstat_qc = ret_val2.get(BIG_INT)

    log.info('Closing multi-threading...')
    pool.close()
    pool.join()

    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()
