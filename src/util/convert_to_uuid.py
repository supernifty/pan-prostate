#!/usr/bin/env python
#
# python src/util/convert_to_uuid.py $SOURCES < $UUIDS > $LINKS
# sources = list of paths to fastq file
# uuids = cmhs,labid
# links = ln commands

import os
import sys

uuid = {}
duplicate_uuids = set()
for line in sys.stdin: # map cmhs to labid
    fields = line.strip('\n').split(',')
    if fields[1] in uuid:
        sys.stderr.write('WARNING: lab {} found more than once in uuid list: {}, {}\n'.format(fields[1], fields[0], uuid[fields[1]]))
        duplicate_uuids.add(fields[0])
    uuid[fields[1]] = fields[0]

filecount = 0
for line in open(sys.argv[1], 'r'):
    fn = line.strip()
    # e.g. labid_wb_R1.fq.gz, labid_wb_x_R1.fq.gz
    components = os.path.basename(fn).split('_')
    #sample = '{}_{}'.format( components[0], components[1] )
    sample = '_'.join(components[:-1])
    if sample in uuid:
        filecount += 1
        sys.stdout.write('ln -s {} {}\n'.format(fn, "{}_{}".format(uuid[sample], components[-1])))
        if uuid[sample] in duplicate_uuids:
            sys.stderr.write('WARNING: data files have multiple possible UUIDs: {} uuid: {} filename: {}\n'.format(sample, uuid[sample], fn))
    else:
        sys.stderr.write('ERROR: {} not found in UUID map\n'.format(fn))

sys.stderr.write('DEBUG: wrote {} symlinks\n'.format(filecount))
