import re
import gzip
import argparse
import os

parser = argparse.ArgumentParser(description = "Calculate dstat",
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--triad', type = str,
                    default = None, help = 'triad for which to run script')
parser.add_argument('--cov', type = int,
                    default = 5, help = 'minimum coverage required at a site')
parser.add_argument('--out', type = str,
                    default = None, help = 'output file (incl. dir)')
args = parser.parse_args()

mdir = '/home/babs/Desktop/birds'
vardir = os.path.join(mdir, 'variants')


def get_var(sp):
    vcf = os.path.join(vardir, '%s.qual_filtered20.cov_filtered5.vcf.gz' % sp)

    f = gzip.open(vcf, 'r')
    for l in f:
        l = l.decode('utf-8')
        if not re.search('^#', l):
            d = re.split('\t', l.rstrip())

            c = d[0]
            pos = d[1]
            a = [d[3]] + re.split(',', d[4])

            dpix = re.split(':', d[8]).index('DP')
            geno = re.split(':', d[9])

            dp = int(geno[dpix])
            geno = [a[int(x)] for x in re.split('/', geno[0])]

            # only want sites with high cov
            if dp >= args.cov:
                if c not in var:
                    var[c] = {}
                if pos not in var[c]:
                    var[c][pos] = []
                var[c][pos].append(geno)
    f.close()

    return var

def print_dstat(var):
    o = open(args.out, 'w')
    o.write('species,locus,aln_pos,abba,baba\n')
    for c in var:
        for pos in var[c]:
            # first check to see if len 4
            # bc that means fully defined across all inds
            if len(var[c][pos]) == 4:
                # now check that biallelic site
                nuc = [x for geno in var[c][pos] for x in geno]
                nuc = list(set(nuc))
                if len(nuc) == 2:

                    if var[c][pos][3].count(nuc[0]) < 2:
                        der = nuc[0]
                    else:
                        der = nuc[1]

                    var2 = []
                    for x in var[c][pos]:
                        af = x.count(der) /  2.0
                        var2.append(af)

                    abba = (1 - var2[0]) * var2[1] * var2[2] * (1 - var2[3])
                    baba = var2[0] * (1 - var2[1]) * var2[2] * (1 - var2[3])
                    
                    if abba > 0 or baba > 0:
                        triad = re.sub(',', ':', args.triad)
                        o.write('%s,%s,%s,%s,%s\n' % (triad, c, pos, abba, baba))
    o.close()

# get the species for which to run
sps = re.split(',', args.triad)
var = {}
for sp in sps:
    var = get_var(sp)
print_dstat(var)
