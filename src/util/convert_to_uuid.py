#!/usr/bin/env python

import os
import sys

uuid = {}
for line in sys.stdin: # map
    fields = line.strip('\n').split(',')
    uuid[fields[1]] = fields[0]

for line in open(sys.argv[1], 'r'):
    fn = line.strip()
    # e.g. a_wb_R1.fq.gz, a_wb_x_R1.fq.gz
    components = os.path.basename(fn).split('_')
    #sample = '{}_{}'.format( components[0], components[1] )
    sample = '_'.join(components[:-1])
    if sample in uuid:
        sys.stdout.write('ln -s {} {}\n'.format(fn, "{}_{}".format(uuid[sample], components[-1])))
    else:
        sys.stderr.write('ERROR: {} not found in UUID map\n'.format(fn))
