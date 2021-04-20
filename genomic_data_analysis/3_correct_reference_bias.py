import argparse
import os
import pandas as pd
import re
import glob
import subprocess
import gzip

"""
Sonal Singhal
created on 22 March 2017
Written assuming:
    * samtools 1.3.1
    * picard 2.4.1
    * bwa
"""

def get_args():
    parser = argparse.ArgumentParser(
        description="Reference correct genome. Written assuming "
                    " samtools 1.3.1, picard 2.4.1, and"
                    " bwa",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

    # sample
    parser.add_argument(
        '--sample',
        type=str,
        default=None,
        help='sample for which to run script.'
        )

    # aligner
    parser.add_argument(
        '--bwa',
        type=str,
        default=None,
        help='bwa executable, full path.'
        )

    # samtools
    parser.add_argument(
        '--samtools',
        type=str,
        default=None,
        help='samtools executable, full path.'
        )

    # bcftools
    parser.add_argument(
        '--bcftools',
        type=str,
        default=None,
        help='bcftools executable, full path.'
        )

    # picard
    parser.add_argument(
        '--picard',
        type=str,
        default=None,
        help='picard executable, full path.'
        )
    
    # CPUs
    parser.add_argument(
        '--CPU',
        type=int,
        default=1,
        help='# of CPUs to use in alignment.'
        )

    # memory
    parser.add_argument(
        '--mem',
        type=int,
        default=1,
        help='Memory available, as an int, in terms of Gb.'
        )
               
    # outdir
    parser.add_argument(
        '--outdir',
        type=str,
        default=None,
        help='Output directory for alignments.'
        )

    # readdir
    parser.add_argument(
        '--readdir',
        type=str,
        default=None,
        help="Full path to files with reads."
        )

    # PRG
    parser.add_argument(
        '--prg',
        type=str,
        default=None,
        help="Full path to pseudoref genome."
        )

    return parser.parse_args()


def get_info(args):
    # get the genome
    genome = args.prg

    outdir = args.outdir
    
    # makes the outdir
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    
    # get the reads
    r1 = os.path.join(args.readdir, '%s_R1.final.fq.gz' % args.sample)
    r2 = os.path.join(args.readdir, '%s_R2.final.fq.gz' % args.sample)
    un = os.path.join(args.readdir, '%s_unpaired.final.fq.gz' % args.sample)
    reads = [r1, r2, un]

    return reads, genome, outdir


def prepare_seq(args, genome):
    # does all the prep necessary for the PRG
    if not os.path.isfile(genome + '.fai'):
        subprocess.call("%s faidx %s" % (args.samtools, genome), shell=True)
    out = re.sub('\.fasta', '.dict', genome)
    if not os.path.isfile(out):
        subprocess.call("java -jar %s CreateSequenceDictionary R=%s O=%s" %
                        (args.picard, genome, out), shell=True)
    subprocess.call("bwa index %s" % genome, shell = True)


def align_seq(args, sample, r, seq, ix, dir):
    
    out1a = os.path.join(dir, '%s_1.sam' % sample)
    out1as = os.path.join(dir, '%s_1.bam' % sample)
    out1b = os.path.join(dir, '%s_2.sam' % sample)
    out2a = os.path.join(dir, '%s.mateFixed.bam' % sample)
    out2b = os.path.join(dir, '%s.bam' % sample)
    out3a = os.path.join(dir, '%s_1.mateFixed.sorted.bam' %  sample)
    out3b = os.path.join(dir, '%s_2.mateFixed.sorted.bam' % sample)
    out4a = os.path.join(dir, '%s_1.rg.mateFixed.sorted.bam' % sample)
    out4b = os.path.join(dir, '%s_2.rg.mateFixed.sorted.bam' % sample)
    out5a = os.path.join(dir, '%s_1.dup.rg.mateFixed.sorted.bam' % sample)
    out5b = os.path.join(dir, '%s_2.dup.rg.mateFixed.sorted.bam' % sample)
    out6 = os.path.join(dir, '%s.dup.rg.mateFixed.sorted.%s.bam' % (sample, ix))
    m_1 = os.path.join(dir, '%s_1.intervals' % sample)
    m_2 = os.path.join(dir, '%s_2.intervals' % sample)
    
    # need a tmpdir for when sorting BAM files
    tmpdir = os.path.join(dir, args.sample)
    if not os.path.isdir(tmpdir):
        os.mkdir(tmpdir)

    # align
    subprocess.call("%s mem -t %s %s %s %s > %s" % 
        (args.bwa, args.CPU, seq, r[0], r[1], out1a), shell=True)
    subprocess.call("%s mem -t %s %s %s > %s" % 
        (args.bwa, args.CPU, seq, r[2], out1b), shell=True)
    # fixmate
    # note that had used samtools, appears to not work properly
    subprocess.call("%s view -@ %s -b %s > %s" % 
        (args.samtools, args.CPU, out1a, out1as), shell=True)
    subprocess.call("java -jar %s FixMateInformation I=%s O=%s" % 
        (args.picard, out1as, out2a), shell=True)
    subprocess.call("%s view -@ %s -b %s > %s" % 
        (args.samtools, args.CPU,  out1b, out2b), shell=True)
    # sorted
    subprocess.call("%s sort -@ %s -O bam -o %s -T %s %s" % 
        (args.samtools, args.CPU, out3a, tmpdir, out2a), shell=True)
    subprocess.call("%s sort -@ %s -O bam -o %s -T %s %s" % 
        (args.samtools, args.CPU, out3b, tmpdir, out1b), shell=True)
    # readgroup
    subprocess.call("java -jar %s AddOrReplaceReadGroups INPUT=%s "
                        "OUTPUT=%s RGLB=%s RGPL=Illumina RGPU=%s RGSM=%s" % 
                          (args.picard, out3a, out4a, sample,
                           sample, sample), shell=True)
    subprocess.call("java -jar %s AddOrReplaceReadGroups INPUT=%s "
                        "OUTPUT=%s RGLB=%s RGPL=Illumina RGPU=%s RGSM=%s" % 
                          (args.picard, out3b, out4b, sample,
                           sample, sample), shell=True)
    # mark read duplicates
    subprocess.call("java -jar %s MarkDuplicates I=%s O=%s ASO=coordinate METRICS_FILE=%s" % 
                       (args.picard, out4a, out5a, m_1), shell=True)
    subprocess.call("java -jar %s MarkDuplicates I=%s O=%s ASO=coordinate METRICS_FILE=%s" % 
                       (args.picard, out4b, out5b, m_2), shell=True)
    subprocess.call("%s merge -@ %s %s %s %s" %
                        (args.samtools, args.CPU, out6, out5a, out5b), shell=True)
    
    # remove the files
    [os.remove(x) for x in [out1a, out1as, out1b, out2a, out2b, out3a, 
                                out3b, out4a, out4b, out5a, out5b, m_1, m_2]]
    
    # remove the dir
    os.rmdir(tmpdir)

    return out6


def get_seq(seqfile):
    id = ''
    seq = {}
        
    f = open(seqfile, 'r')
    for l in f:
        if re.search('>', l):
            id = re.search('>(\S+)', l).group(1)
            seq[id] = ''
        else:
            seq[id] += l.upper().rstrip()
    f.close()

    for id, s in seq.items():
        seq[id] = list(s)

    return seq


def call_snps(args, bam, genome, ix, outdir):
    
    out1 = os.path.join(outdir, '%s.%s.bcf' % (args.sample, ix))
    out2 = os.path.join(outdir, '%s.%s.vcf.gz' % (args.sample, ix))

    if not os.path.isfile(out2):
        subprocess.call("%s mpileup -ABI -f %s -o %s %s" % 
            (args.bcftools, genome, out1, bam), shell=True)
        # output variant sites only
        subprocess.call("%s call -mvO z -o %s %s" % 
            (args.bcftools, out2, out1), shell=True)

    # get in seq
    seq = get_seq(genome)

    # mutate seq
    i = gzip.open(out2, 'r')
    for l in i:
        l = l.decode('utf-8')
        if not re.search('#', l):
            d = re.split('\t', l.rstrip())
            # check biallelic
            if d[3] in ['A', 'T', 'C', 'G'] and d[4] in ['A', 'T', 'C', 'G']:
                # get af
                genos = d[9:]
                genos = [re.search('(\S/\S)', x).group(1) for x in genos]
                genos = [re.split('/', x) for x in genos]
                genos = [x for ind in genos for x in ind]
                genos = [x for x in genos if x != '.']

                if len(genos) > 0:
                    af = genos.count('1') / float(len(genos))
                    # mutate seq
                    if af >= 0.5:
                        pos = int(d[1]) - 1
                        seq[d[0]][pos] = d[4]
    i.close()

    out = os.path.join(outdir, '%s.%s.fasta' % (args.sample, ix))
    o = open(out, 'w')
    for id, s in seq.items():
        o.write('>%s\n%s\n' % (id, ''.join(s)))
    o.close()

    os.remove(out1)
    os.remove(out2)
    
    return out


def main():
    # get arguments
    args = get_args()
    reads, genome, outdir = get_info(args)

    ind = args.sample
    # round 1
    # prep sequence
    prepare_seq(args, genome)
    # do the alignments
    bamout = align_seq(args, ind, reads, genome, 1, outdir)
    genome1 = call_snps(args, bamout, genome, 1, outdir)

    # round2
    # prep sequence
    prepare_seq(args, genome1)
    # do the alignments
    bamout = align_seq(args, ind, reads, genome1, 2, outdir)
    genome2 = call_snps(args, bamout, genome1, 2, outdir)

    # round2
    # prep sequence
    prepare_seq(args, genome2)
    # do the alignments
    bamout = align_seq(args, ind, reads, genome2, 3, outdir)
    genome3 = call_snps(args, bamout, genome2, 3, outdir)

    # round4
    # prep sequence
    prepare_seq(args, genome3)
    # do the alignments
    bamout = align_seq(args, ind, reads, genome3, 4, outdir)
    genome4 = call_snps(args, bamout, genome3, 4, outdir)

    # rename final
    final = re.sub('\.4\.', '.', genome4)
    subprocess.call("cp %s %s" % (genome4, final), shell=True)

    # probably need to do cleanup
    for end in ['dict', 'amb', 'ann', 'bwt', 'fai', 'pac', 'sa']:
        ixfiles = glob.glob('%s/%s*%s' % (outdir, ind, end))
        [os.remove(f) for f in ixfiles]
    bams = glob.glob('%s/%s*bam' % (outdir, ind))
    [os.remove(f) for f in bams]
    [os.remove(f) for f in [genome1, genome2, genome3, genome4]]

if __name__ == "__main__":
    main()
