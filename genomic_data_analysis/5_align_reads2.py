import argparse
import os
import pandas as pd
import re
import subprocess

"""
Sonal Singhal
created on 22 July 2020
Written assuming:
    * GATK 3.6
"""

def get_args():
    parser = argparse.ArgumentParser(
        description="Align reads to sample, step 2. "
                    " Uses GATK 4",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

    # lineage
    parser.add_argument(
        '--sample',
        type=str,
        default=None,
        help='sample for which to run script.'
        )

    # file
    parser.add_argument(
        '--file',
        type=str,
        default=None,
        help='File with sample info.'
        )
                
    # basedir
    parser.add_argument(
        '--dir',
        type=str,
        default=None,
        help="Full path to base dir with reads & assemblies "
             "everything else."
        )

    # GATK
    parser.add_argument(
        '--gatk',
        type=str,
        default=None,
        help='GATK executable, full path.'
        )
    
    # memory
    parser.add_argument(
        '--mem',
        type=int,
        default=1,
        help='Memory available, as an int, in terms of Gb.'
       )

    # qual
    parser.add_argument(
        '--qual',
        type=int,
        default=20,
        help='Minimum quality to retain variant for '
            'creating validated call set.'
        )

    # depth
    parser.add_argument(
        '--dp',
        type=int,
        default=10,
        help='Minimum depth required per individual to retain '
             'variant for creating validated call set.'
        )

    # samtools
    parser.add_argument(
        '--samtools',
        type=str,
            default = None,
        help='path to samtools'
        )
        
    
    # outdir
    parser.add_argument(
        '--outdir',
        type=str,
        default=None,
        help='Output directory for alignments, only needed '
             'if not running in context of pipeline.'
        )
                
    # bamfiles
    parser.add_argument(
        '--bamfile',
        type=str,
        default=None,
        help="Full path to BAM file"
        )

    # PRG
    parser.add_argument(
        '--prg',
        type=str,
        default=None,
        help="Full path to pseudoref genome if "
             "you aren't running in context of pipeline."
        )

    return parser.parse_args()


def get_files(args):
    # gets the bam files
    if args.bamfile:
        bam = args.bamfile
    else:
        bam = os.path.join(args.dir, 'alignments', 
            '%s.dup.rg.mateFixed.sorted.bam' % args.sample)
    subprocess.call("%s index %s" % (args.samtools, bam), shell = True)

    # gets the genome
    if args.prg:
        genome = args.prg
    else:
        genome = os.path.join(args.dir, 'ref_bias', '%s.fasta' % args.sample)

    # gets the outdir
    if args.outdir:
        outdir = args.outdir
    else:
        outdir = os.path.join(args.dir, 'alignments')
    
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    return bam, genome, outdir


def get_qual(args, bam, genome, dir):
    raw_vcf = os.path.join(dir, '%s.raw.vcf' % args.sample)
    filt_vcf = os.path.join(dir, '%s.filt.vcf' % args.sample)

    # makes the raw VCFs, only outputting variant SNPs
    subprocess.call("%s --java-options \"-Xmx%sg\" HaplotypeCaller  "
        "-R %s "
        "-I %s "
        "-O %s "
        "--output-mode EMIT_VARIANTS_ONLY" % 
        (args.gatk, args.mem, genome, bam, raw_vcf), shell = True)

    f = open(raw_vcf, 'r')
    o = open(filt_vcf, 'w')

    for l in f:
        if re.match('^#', l):
            o.write(l)
        else:
            d = re.split('\t', l.rstrip())
            
            # check if indel
            alleles = [d[3]] + re.split(',', d[4])
            indel = False
            for a in alleles:
                if len(a) > 1:
                    indel = True
            
            # depth cover
            n_inds = len(d[9:])
            dp = args.dp * n_inds
            snp_depth = int(re.search('DP=(\d+)', d[7]).group(1))

            if not indel and float(d[5]) >= args.qual and snp_depth >= dp:
                o.write(l)

    f.close()
    o.close()

    os.remove(raw_vcf)
    os.remove(raw_vcf + '.idx')

    subprocess.call("%s IndexFeatureFile --input %s" % (args.gatk, filt_vcf), shell = True)

    return filt_vcf


def recalibrate(args, bam, genome, vcf, dir):

    stem = bam.replace('.bam', '')
    out = stem + '.recal.bam'
    recal = '%s.recal.table' % stem
    
    # generate recal table
    subprocess.call("%s --java-options \"-Xmx%sg\" BaseRecalibrator  "
        "-R %s "
        "-I %s "
        "-O %s "
        "--known-sites %s" % 
        (args.gatk, args.mem, genome, bam, recal, vcf), shell = True)

    # print the recal reads
    subprocess.call("%s --java-options \"-Xmx%sg\" ApplyBQSR  "
        "-R %s "
        "-I %s "
        "--bqsr-recal-file %s "
        "-O %s" %
        (args.gatk, args.mem, genome, bam, recal, out), shell = True)

    os.remove(recal)
    os.remove(bam)
    os.remove(bam + '.bai')
    
    os.remove(vcf)
    os.remove(vcf + '.idx')


def main():
    args = get_args()
    files, genome, outdir = get_files(args) 
    vcf = get_qual(args, files, genome, outdir)
    recalibrate(args, files, genome, vcf, outdir)


if __name__ == "__main__":
    main()
