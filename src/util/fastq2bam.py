#!/usr/bin/env python

import collections
import gzip
import logging
import os
import re
import sys
import threading

ROOT='/data/projects/punim0095/pan-prostate'

def read_group_info(filename):
    logging.info('processing {}...'.format(filename))
    lanes = collections.defaultdict(dict)
    flowcells = collections.defaultdict(dict)
    idx = 0
    for idx, line in enumerate(gzip.open(filename, 'r')):
        if idx % 4 == 0:
            # e.g. @HWI-ST960:132:D1MAJACXX:1:1101:1237:2215 1:N:0:GCCAAT
            # based on https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers
            # 0: HWI-ST960: instrument name
            # 1: 132: run ID
            # 2: D1MAJACXX: flowcell ID
            # 3: 1: flowcell lane
            # 4: 1101: number in lane
            # 5: 1237: x
            # 6: 2215: y

            # 7: 1: paired end
            # 8: N: filtered (Y) or not (N)
            # 9: 0: control bits
            # 10: GCCAAT: index sequence
            fields = re.split(':| ', line.strip('\n'))
            if len(fields) != 11:
                logging.info('WARN: {} had unexpected number of fields {}'.format(filename, fields))
            else:
                flowcell = fields[2]
                lane = fields[3]
                barcode = fields[10]
                if barcode not in lanes[lane]:
                    lanes[lane][barcode] = 0
                lanes[lane][barcode] += 1
                if flowcell not in flowcells[lane]:
                    flowcells[lane][flowcell] = 0
                flowcells[lane][flowcell] += 1
        if idx % 10000000 == 0:
            logging.info('processed {} million lines. {} lanes'.format(idx / 1000000, len(lanes)))
    # finished with file
    result = {}
    for lane in lanes:
        total = 0
        best = None
        best_flowcell = None
        for barcode in lanes[lane]:
            total += lanes[lane][barcode]
            if best is None or lanes[lane][barcode] > lanes[lane][best]:
                best = barcode
        for flowcell in flowcells[lane]:
            if best_flowcell is None or flowcells[lane][flowcell] > flowcells[lane][best_flowcell]:
                best_flowcell = flowcell
        result[lane] = { 'barcode': best, 'flowcell': best_flowcell }

    return result

def fastqsplit(lanes, inputs):
    command = '{}/src/util/fastqsplit {}'.format(ROOT, ' '.join(inputs))
    result = os.system(command)
    if result != 0:
        raise "command {} failed: {}".format(command, result)

def fastq_to_merged_bam(fastq_read1_in, fastq_read2_in, output_dir, bam_out, use_same_library=True):
    '''
      processes the specific fastq read pair and looks for additional sequencing
      of the form sample-n_R[12].fastq.gz
    '''
    # 1. find lanes and barcodes
    # 2. split into lanes
    # 3. fastq2bam
    sample = os.path.basename(fastq_read1_in).split('_')[0]
    bams = fastq_to_bams(fastq_read1_in, fastq_read2_in, output_dir, sample, use_same_library=use_same_library)

    # is there more sequencing of the form sample-n?
    index = 1
    while True:
        candidate_fastq_read1 = fastq_read1_in.replace(sample, '{}-{}'.format(sample, index))
        if os.path.isfile(candidate_fastq_read1):
            candidate_fastq_read2 = fastq_read2_in.replace(sample, '{}-{}'.format(sample, index))
            bams.extend(fastq_to_bams(candidate_fastq_read1, candidate_fastq_read2, output_dir, sample, use_same_library=use_same_library))
            index += 1
        else:
            break
    
    # 4. merge bams
    command = 'java -jar {}/tools/picard-2.8.2.jar MergeSamFiles {inputs} OUTPUT={bam_out} ASSUME_SORTED=true'.format(ROOT, inputs=' '.join([ 'INPUT={}'.format(bam) for bam in bams]), bam_out=bam_out)
    logging.info("merging bams: {}".format(command))
    result = os.system(command)
    if result != 0:
        raise "command {} failed: {}".format(command, result)

def fastq_to_bams(fastq_read1_in, fastq_read2_in, output_dir, parent_sample, use_same_library=True):
    '''
      processes a specific fastq read pair, splits into lanes, runs fastq2bam
      returns list of generated bams
    '''
    # 1. find lanes and barcodes
    lanes = read_group_info(fastq_read1_in)

    # 2. split into lanes
    logging.info("splitting into lanes: {}".format((fastq_read1_in, fastq_read2_in)))
    #extract_all_lanes(lanes, (fastq_read1_in, fastq_read2_in), output_dir)
    fastqsplit(lanes, (fastq_read1_in, fastq_read2_in))

    # 3. fastq2bam
    logging.info("converting to bam: {}".format((fastq_read1_in, fastq_read2_in)))
    bams = []

    # cmhs sample format
    sample = os.path.basename(fastq_read1_in).split('_')[0]

    for lane in lanes.keys():
        # filename has to match util.extract_all_lanes
        #r1 = os.path.join(output_dir, '{sample}_{flowcell}_{barcode}_L00{lane}_R1.fastq.gz'.format(sample=sample, flowcell=lanes[lane]['flowcell'], barcode=lanes[lane]['barcode'], lane=lane))
        #r2 = os.path.join(output_dir, '{sample}_{flowcell}_{barcode}_L00{lane}_R2.fastq.gz'.format(sample=sample, flowcell=lanes[lane]['flowcell'], barcode=lanes[lane]['barcode'], lane=lane))
        #bam = os.path.join(output_dir, '{sample}_{flowcell}_{barcode}_L00{lane}.bam'.format(sample=sample, flowcell=lanes[lane]['flowcell'], barcode=lanes[lane]['barcode'], lane=lane))
        os.rename(os.path.join(os.path.dirname(fastq_read1_in), '{sample}_{lane}_R1.fastq.gz'.format(sample=sample, lane=lane)), os.path.join(output_dir, '{sample}_{lane}_R1.fastq.gz'.format(sample=sample, lane=lane)))
        os.rename(os.path.join(os.path.dirname(fastq_read2_in), '{sample}_{lane}_R2.fastq.gz'.format(sample=sample, lane=lane)), os.path.join(output_dir, '{sample}_{lane}_R2.fastq.gz'.format(sample=sample, lane=lane)))
        r1 = os.path.join(output_dir, '{sample}_{lane}_R1.fastq.gz'.format(sample=sample, lane=lane))
        r2 = os.path.join(output_dir, '{sample}_{lane}_R2.fastq.gz'.format(sample=sample, lane=lane))
        bam = os.path.join(output_dir, '{sample}_{flowcell}_{barcode}_L00{lane}.bam'.format(sample=sample, flowcell=lanes[lane]['flowcell'], barcode=lanes[lane]['barcode'], lane=lane))

        if use_same_library:
            library = parent_sample
        else:
            library = sample

        command = '{command_dir}/fastqtobam gz=1 namescheme=pairedfiles qualitymax=42 I={fastq_r1} I={fastq_r2} RGID={rgid} RGPU={rgpu} RGSM={rgsm} RGPL={rgpl} RGLB={rglb} verbose=1 > {bam}'.format(
            command_dir="{}/tools/biobambam2-2.0.65-release-20161130121735-x86_64-etch-linux-gnu/bin".format(ROOT), 
            fastq_r1=r1, 
            fastq_r2=r2, 
            bam=bam,
            rgid='{sample}-{flowcell}-{barcode}.{lane}'.format(sample=sample, flowcell=lanes[lane]['flowcell'], barcode=lanes[lane]['barcode'], lane=lane), 
            rgpu='{flowcell}-{barcode}.{lane}'.format(flowcell=lanes[lane]['flowcell'], barcode=lanes[lane]['barcode'], lane=lane),
            rgsm=parent_sample, # same across each input
            rgpl='ILLUMINA', 
            rglb='lib-{}'.format(library)) # could potentially be different
        logging.info("fastq2bam: {}".format(command))
        result = os.system(command)
        if result != 0:
            raise "command {} failed: {}".format(command, result)
        bams.append(bam)

    return bams

if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG)
    logging.info('starting fastq2bam with parameters {}...'.format(sys.argv))
    import argparse
    parser = argparse.ArgumentParser(description='Prealign')
    parser.add_argument('--r1', required=True, help='read 1')
    parser.add_argument('--r2', required=True, help='read 2')
    parser.add_argument('--output_dir', required=True, help='output directory')
    parser.add_argument('--bam', required=True, help='output bam')
    args = parser.parse_args()
    fastq_to_merged_bam(args.r1, args.r2, args.output_dir, args.bam, use_same_library=False)
    logging.info('fastq2bam with parameters {}: finished'.format(sys.argv))
