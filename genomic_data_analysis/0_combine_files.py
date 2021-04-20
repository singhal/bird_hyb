import re
import glob
import os
import subprocess

r = {}
f = open('/home/babs/Desktop/birds/metadata/file_info.csv', 'r')
for l in f:
    d = re.split(',', l.rstrip())
    if d[2] != 'SampleName':
        if d[2] not in r:
            r[d[2]] = []
        read = d[0]
        if re.search('_1_Rapid', read):
            r[d[2]].append(read)
f.close()

curdir = '/home/babs/Desktop/birds/Raw-Data/'
newdir = '/home/babs/Desktop/birds/Raw-Data-renamed/'

for ind in r:
    new1 = os.path.join(newdir, '%s_R1.fastq.gz' % ind)
    new2 = os.path.join(newdir, '%s_R2.fastq.gz' % ind)

    old1 = [os.path.join(curdir, read) for read in r[ind]]
    old2 = [re.sub('_1_Rapid', '_2_Rapid', read) for read in old1]

    old1ok = True
    old2ok = True
    for f in old1:
        if not os.path.isfile(f):
            old1ok = False
            print('%s\t%s' % (ind, f))
    for f in old2:
        if not os.path.isfile(f):
            old2ok = False
            print('%s\t%s' % (ind, f))

    if old1ok:
        subprocess.call("cat %s > %s" % (' '.join(old1), new1), shell = True)
        for f in old1:
            os.remove(f)
    if old2ok:
        subprocess.call("cat %s > %s" % (' '.join(old2), new2), shell = True)
        for f in old2:
            os.remove(f)
