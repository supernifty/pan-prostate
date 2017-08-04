#!/usr/bin/env python
# find donors
# takes list of CMHS filenames, find donors in metadata
# find /data/punim0261/data01/out/ -type l -name \*.fastq.gz | sed 's/.*\///g' > samples.txt
# usage: python sample_counts.py ../../cfg/sample-metadata.csv < samples.txt

import sys

#Sample UUID,Patient UUID,Lab ID,tissue_id,is_normal,Status,Pre-aligned,Validated
m = {} # sample -> is_normal
wbm = {}
for line in open(sys.argv[1], 'r'):
  f = line.strip('\n').split(',')
  if f[0] in m:
    print('warning: sample {} defined more than once in metadata'.format(f[0]))
  m[f[0]] = f[1] # sample -> donor
  wbm[f[0]] = f[4] # sample -> isnormal

d = set()
counts = {'N': 0, 'Y': 0, 'Samples': 0}
wbd = set()
files = set()
for line in sys.stdin: # available sample with filename
  x = line.strip('\n').replace('-', '_')
  #x = line.strip('\n')
  if '_R1' in x:
    name = x.split('_')
    if name[0] in files:
      print('warning: filename {} found more than once on filesystem'.format(name[0]))
    files.add(name[0])

    d.add(m[name[0]]) # donors
    counts[wbm[name[0]]] += 1
    if wbm[name[0]] == 'Y': # is normal
      if m[name[0]] in wbd:
        print("{} has multiple normal".format(name[0]))
      wbd.add(m[name[0]])
    counts['Samples'] += 1

counts['Donors'] = len(d)

# write donors found
#sys.stdout.write('\n'.join(list(d)))
#sys.stdout.write('\n'.join(list(wbd)))

# write wb vs tumous
print(counts)
