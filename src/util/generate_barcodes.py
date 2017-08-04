#!/usr/bin/env python
'''
  find inferred flowcell and barcode for all samples
'''

import glob
import logging
import os
import re
import subprocess
import sys

logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)

sys.stdout.write('#Sample,Barcode,Flowcell\n')
found = 0
for f in glob.glob('CMHS[0-9]*.mapped.bam'): # pre-aligned bam is at least in progress
  logging.info('processing %s', f)
  for line in subprocess.check_output(['samtools', 'view', '-H', f]).split('\n'):
    if line.startswith('@RG'):
      fields = line.strip('\n').split('ID:')[1].split('\t')[0].split('-') # @RG     ID:CMHS100-H2CWCALXX-CGGCTATG.5 SM:b0be11c6-f452-4752-894a-868a8f45202b PL:ILLUMINA
      sys.stdout.write('{},{},{}\n'.format(f.split('.')[0], fields[2].split('.')[0], fields[1]))
      logging.info('processed %s', f)
      break
    else:
      pass
      #logging.debug('skipped %s', line)
