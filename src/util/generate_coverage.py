#!/usr/bin/env python
'''
  collect coverage details
'''

import glob
import logging
import os
import re
import subprocess
import sys

logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)

sys.stdout.write('#Sample,Name,Value\n')
found = 0
#CMHS70.genomecov.stats
# looks like:
#coverage_total: 404868170538
#coverage_bases: 3137454505
#coverage_mean: 129.043519163
#coverage_total_no_zero: 404868170538
#coverage_mean_no_zero: 139.722244685

#CMHS70.picard.stats
# looks like:
#GENOME_TERRITORY        MEAN_COVERAGE   SD_COVERAGE     MEDIAN_COVERAGE MAD_COVERAGE    PCT_EXC_MAPQ    PCT_EXC_DUPE    PCT_EXC_UNPAIRED        PCT_EXC_BASEQ   PCT_EXC_OVERLAP PCT_EXC_CAPPED  PCT_EXC_TOTAL   PCT_1X  PCT_5X  PCT_10X PCT_15X PCT_20X PCT_25X PCT_30X PCT_40X PCT_50X PCT_60X PCT_70X PCT_80X PCT_90X PCT_100X        HET_SNP_SENSITIVITY     HET_SNP_Q
#2900340137      124.041702      305.33184       122     12      0       0.096314        0.000593        0.000205        0.014019        0.000108        0.111239        0.999071        0.998584        0.997964        0.997164        0.996558        0.99594 0.995112        0.992166        0.983791        0.959589        0.922291        0.888627        0.862206        0.82599 0.99801 27

for f in glob.glob('CMHS[0-9]*.picard.stats'): # pre-aligned bam is at least in progress
  logging.info('processing %s', f)
  stage = 0
  for line in open(f, 'r'):
    if stage == 0:
      if line.startswith('## METRICS CLASS'):
        stage = 1
    elif stage == 1:
      #logging.debug('stage 1')
      stage = 2
      fieldnames = line.strip('\n').split('\t')
    elif stage == 2:
      #logging.debug('stage 2')
      fields = line.strip('\n').split('\t')
      if len(fields) != len(fieldnames):
        logging.warn('Fieldname count mismatch: %i vs %i', len(fields), len(fieldnames))
        break
      for idx, fieldname in enumerate(fieldnames):
        sys.stdout.write('{},picard.{},{}\n'.format(f.split('.')[0], fieldname, fields[idx]))
      break

