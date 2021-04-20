import re
import gzip
import glob
import pandas as pd
import os

basedir = '/home/babs/Desktop/birds/'
# get sample to dir mapping
mdir = os.path.join(basedir, 'metadata')
vcf1 = os.path.join(mdir, 'unique_vcf_files_nodup.csv')
v = pd.read_csv(vcf1)
vcf = {}
for ix, row in v.iterrows():
    vcf[row['sample_name']] = re.sub('[^/]+$', '', row['file_path'])
outdir = os.path.join(basedir, 'PRG')
    
# mixed mapping files
mmdir = os.path.join(mdir, 'mappings')
maps = glob.glob(mmdir + '/*csv')
for mapf in maps:
    sample = re.sub('^.*\/', '', mapf)
    sample = re.sub('_mixed_vcf.*', '', sample)
    if sample not in vcf:
        vcf[sample] = mapf
    else:
        print('%s is both unique & mixed?!' % sample)

# get base directory for each individual
basedir2 = os.path.join(basedir, 'diploid-processing')

def get_seq(seqfile):
    seq = {}
    id = ''
    f = gzip.open(seqfile, 'r')
    for l in f:
        l = l.decode('utf-8') 
        if re.search('>', l):
            id = re.search('>(\S+)\|', l).group(1)
            seq[id] = ''
        else:
            seq[id] += l.rstrip()
    f.close()
    return seq

                    
for sp in vcf:
    if not re.search('csv', vcf[sp]):
        id = re.sub('^[^/]+', '', vcf[sp])
        id = re.sub('\/', '', id)
        seqf = os.path.join(basedir2, vcf[sp], '%s_contigs.fasta.gz' % id)
        if os.path.isfile(seqf):
            indseq = get_seq(seqf)
        else:
            indseq = {}
    else:
        # otherwise ...
        d = pd.read_csv(vcf[sp])
        # what are my unique outputs
        uniqvcf = d.file_path.unique().tolist()
        uniqseq = [re.sub('_SNPs_phased.vcf', '_contigs.fasta.gz', x) for x in uniqvcf] 

        d['seq_file'] = [re.sub('_SNPs_phased.vcf', '_contigs.fasta.gz', x) for x in d.file_path]
        mid = dict(zip(d.locus_name, d.seq_file))

        indseqs = {}
        indseq = {}
        for seqf in uniqseq:
            seqf2 = os.path.join(basedir2, seqf)
            if os.path.isfile(seqf2):
                indseqs[seqf] = get_seq(seqf2)

                
        # have to do this carefully bc
        # it turns out we don't have all the data
        for loc, seqf in mid.items():
            if mid[loc] in indseqs:
                indseq[loc] = indseqs[mid[loc]][loc]
        
    seqout = os.path.join(outdir, '%s.fasta' % sp)
    o = open(seqout, 'w')
    for s, id in indseq.items():
        o.write('>%s\n%s\n' % (s, id))
    o.close()
