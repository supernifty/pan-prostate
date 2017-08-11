#!/usr/bin/env python
'''
  find status of all samples for htsdb
'''

import glob
import os
import re
import sys

sys.stdout.write('#Sample,Status,Date\n')
for f in glob.glob('CMHS[0-9]*'): # pre-aligned bam is at least in progress
  # prealigned CMHSnnn.bam, CMHSnnn.validation
  m = re.search('^(CMHS[0-9]+).bam$', f)
  if m is not None:
    modified = os.stat(f).st_mtime
    sys.stdout.write('{},prealigned,{}\n'.format(m.group(1), modified))
    continue

  # aligned CMHSnnn.mapped.bam
  m = re.search('^(CMHS[0-9]+).mapped.bam$', f)
  if m is not None:
    modified = os.stat(f).st_mtime
    sys.stdout.write('{},aligned,{}\n'.format(m.group(1), modified))
    continue

  # aligned with qc - updated with populate_measurement

# variant calling - CMHSnnn.wgs/completed.finished
for f in glob.glob('CMHS[0-9]*.wgs/completed.finish'): # pre-aligned bam is at least in progress
  sample = f.split('.')[0]
  modified = os.stat(f).st_mtime
  sys.stdout.write('{},variants,{}\n'.format(sample, modified))

# delly
for f in glob.glob('CMHS[0-9]*.delly.completed'): # pre-aligned bam is at least in progress
  sample = f.split('.')[0]
  modified = os.stat(f).st_mtime
  sys.stdout.write('{},delly,{}\n'.format(sample, modified))
