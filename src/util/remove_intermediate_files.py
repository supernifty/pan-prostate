#!/usr/bin/env python

# -rw-r--r-- 1 ubuntu ubuntu  5205054996 Mar 31 04:01 CMHS51_6_R1.fastq.gz
#-rw-r--r-- 1 ubuntu ubuntu  6055631947 Mar 31 04:25 CMHS51_6_R2.fastq.gz
#-rw-r--r-- 1 ubuntu ubuntu  5064887984 Mar 31 04:01 CMHS51_7_R1.fastq.gz
#-rw-r--r-- 1 ubuntu ubuntu  5893668696 Mar 31 04:25 CMHS51_7_R2.fastq.gz
#-rw-r--r-- 1 ubuntu ubuntu  4888585183 Mar 31 04:00 CMHS51_8_R1.fastq.gz
#-rw-r--r-- 1 ubuntu ubuntu  5687766745 Mar 31 04:25 CMHS51_8_R2.fastq.gz
#-rw-r--r-- 3 ubuntu ubuntu 31734848002 Mar 31 10:10 CMHS51.bam
#-rw-r--r-- 1 ubuntu ubuntu       78523 Mar 31 10:10 CMHS51.bam.log.err
#-rw-r--r-- 1 ubuntu ubuntu         174 Mar 31 10:10 CMHS51.bam.log.out
#-rw-rw-r-- 1 ubuntu ubuntu        1110 Apr  3 03:52 CMHS51.dockstore
#-rw-r--r-- 1 ubuntu ubuntu 10592468572 Mar 31 07:33 CMHS51_HCC3YALXX_CGCTCATT_L006.bam
#-rw-r--r-- 1 ubuntu ubuntu 10311402232 Mar 31 06:26 CMHS51_HCC3YALXX_CGCTCATT_L007.bam
#-rw-r--r-- 1 ubuntu ubuntu  9962544131 Mar 31 05:22 CMHS51_HCC3YALXX_CGCTCATT_L008.bam
#-rw-r--r-- 1 ubuntu ubuntu 27236020701 Apr  3 14:52 CMHS51.mapped.bam
#-rw-r--r-- 1 ubuntu ubuntu     9584272 Apr  3 14:48 CMHS51.mapped.bam.bai
#-rw-r--r-- 1 ubuntu ubuntu        1401 Apr  3 14:52 CMHS51.mapped.bam.bas
#-rw-r--r-- 1 ubuntu ubuntu        6766 Apr  3 14:48 CMHS51.mapped.bam.log.err
#-rw-r--r-- 1 ubuntu ubuntu       62419 Apr  3 14:52 CMHS51.mapped.bam.log.out
#-rw-r--r-- 1 ubuntu ubuntu         348 Apr  3 14:48 CMHS51.mapped.bam.maptime
#-rw-r--r-- 1 ubuntu ubuntu          32 Apr  3 14:52 CMHS51.mapped.bam.md5
#-rw-r--r-- 1 ubuntu ubuntu        1576 Apr  3 14:48 CMHS51.mapped.bam.met
#-rw-r--r-- 1 ubuntu ubuntu         306 Apr  2 18:46 CMHS51.validation
#-rw-r--r-- 1 ubuntu ubuntu        1122 Apr  2 18:46 CMHS51.validation.err
#-rw-r--r-- 1 ubuntu ubuntu           0 Apr  2 17:58 CMHS51.validation.out
#-rw-rw-r-- 1 ubuntu ubuntu         102 Apr  2 17:58 CMHS51.validation_src

import os
import glob

total_size = 0

# look for split fastqs
for f in glob.glob('CMHS[0-9]*.bam'): # pre-aligned bam is at least in progress
    components = f.split('.')
    if len(components) != 2:
        continue
    split_fastq = '{}_[0-9]_*.fastq.gz'.format(components[0])
    #split_fastq = '{}_*.fastq.gz'.format(components[0])
    for c in glob.glob(split_fastq):
        this_size = os.stat(c).st_size
        total_size += this_size
        print("can remove since {} exists saves {} of {}: {}".format(f, this_size, total_size, c))

# look for intermediate bams
for f in glob.glob('CMHS[0-9]*.mapped.bam'): # pre-aligned bam is at least in progress
    components = f.split('.')
    if len(components) != 3:
        continue
    lane_bams = '{}_*.bam'.format(components[0])
    #print("looking for {}".format(lane_bams))
    for c in glob.glob(lane_bams):
        this_size = os.stat(c).st_size
        total_size += this_size
        print("can remove since {} exists saves {} of {}: {}".format(f, this_size, total_size, c))

print("can save {} bytes".format(total_size))
