import argparse
import os
import gzip
import pandas as pd
import re
import subprocess

"""
Sonal Singhal
created on 22 June 2016
Written assuming:
    * GATK 3.6
"""

def get_args():
    parser = argparse.ArgumentParser(
        description="Call SNPs. Assumes bcftools",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # sample
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
        help="Full path to base dir with reads & assemblies & "
             "everything else."
        )
        
    # bcftools
    parser.add_argument(
        '--bcftools',
        type=str,
        default=None,
        help='bcftools executable, full path.'
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
                     'creating final call set.'
        )

    # depth
    parser.add_argument(
        '--dp',
        type=int,
        default=10,
        help='Minimum depth to retain variant for '
                     'creating final call set.'
        )
        
    # outdir
    parser.add_argument(
        '--outdir',
        type=str,
        default=None,
        help='Output directory for variants, '
             ' only needed if running out of pipeline.'
        )

    # bamfiles
    parser.add_argument(
        '--bamfile',
        type=str,
        default=None,
        help="Full path to BAM file, "
             "if running not in context of pipeline. "
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

    # get the bam files
    if args.bamfile:
        file = args.bamfile
    else:
        d = pd.read_csv(args.file)
        sample = args.sample
        file = os.path.join(args.dir, 'alignments', 
                        '%s.dup.rg.mateFixed.sorted.recal.bam' % sample)

    # get the prg associated with sample
    if args.prg:
        genome = args.prg
    else:
        genome = os.path.join(args.dir, 'ref_bias', '%s.fasta' % args.sample)

    # find the outdir
    if args.outdir:
        outdir = args.outdir
    else:
        outdir = os.path.join(args.dir, 'variants')

    # make the outdir if it doesn't exist
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    
    return file, genome, outdir


def get_vcf(args, bam, seq, dir):
    vcf_out1 = os.path.join(dir, '%s.raw.bcf' % args.sample)
    vcf_out2 = os.path.join(dir, '%s.raw.vcf.gz' % args.sample)

    subprocess.call("%s mpileup -ABI -a 'DP' -O b -f %s -o %s %s" % 
            (args.bcftools, seq, vcf_out1, bam), shell=True)
    subprocess.call("%s call -mO z -o %s %s" % 
            (args.bcftools, vcf_out2, vcf_out1), shell=True)
    
    return vcf_out2


def depth_filter(args, infile, dir):
    out = os.path.join(dir, '%s.qual_filtered%s.cov_filtered%s.vcf' % (args.sample, args.qual, args.dp))

    f = gzip.open(infile, 'r')
    o = open(out, 'w')

    for l in f:
        l = l.decode('utf-8')
        if re.search('^#', l):
            o.write(l)
        else:
            d = re.split('\t', l.rstrip())
            # only retain HQ sites

            qual = float(d[5])
            if qual >= args.qual:
                # the depth tag moves around
                # so find out where it is
                tags = re.split(':', d[8])
                depth = tags.index('DP')
                genos = d[9:]
                miss = 0
                for ix, gen in enumerate(genos):
                    info = re.split(':', gen)
                    # some sites will be missing already
                    if info[0] == './.':
                        miss += 1
                    # if too low, set to missing
                    elif int(info[depth]) < args.dp:
                        d[ix + 9] = re.sub('^\S/\S', './.', d[ix + 9])
                        miss += 1
                # only retain the site if someone was genotyped at it
                if miss < len(genos):
                    o.write('\t'.join(d) + '\n')
    f.close()
    o.close()

    # os.remove(infile)
    # os.remove(infile + '.idx')
    
    # gzip the file
    subprocess.call("gzip %s" % (out), shell=True)

def main():
    args = get_args()
    files, seq, dir = get_files(args)
    # get vcf and filter for qual
    out = get_vcf(args, files, seq, dir)
    # filter vcf for depth
    depth_filter(args, out, dir)


if __name__ == '__main__':
    main()
