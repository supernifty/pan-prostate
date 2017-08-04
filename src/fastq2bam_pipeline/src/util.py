
import collections
import gzip
import logging
import os
import re
import sys
import threading

def find_normal(tumour_id, metadata_fh):
    normals = {}
    tumour_donor = None
    # Sample UUID,Patient UUID,Lab ID,tissue_id,is_normal,Status,Pre-aligned,Validated
    for line in metadata_fh:
        fields = line.strip('\n').split(',')
        sample_id = fields[0]
        donor_id = fields[1]
        is_normal = fields[4]
        if tumour_id == sample_id:
            if is_normal == 'Y':
                return None # not a tumour sample
            tumour_donor = fields[1]
        elif is_normal == 'Y':
            normals[fields[1]] = fields[0]

    # have now read the whole file
    if tumour_donor is None: # tumour not found
        raise Exception('failed to find sample {} in sample metadata file'.format(tumour_id))

    if tumour_donor in normals:
        return normals[tumour_donor]

    raise Exception('no normal found for sample {} with donor {}'.format(tumour_id, tumour_donor))


def read_group_info(filename):
    sys.stderr.write('processing {}...\n'.format(filename))
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
                sys.stderr.write('WARN: {} had unexpected number of fields {}\n'.format(filename, fields))
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
            sys.stderr.write('processed {} million lines. {} lanes\n'.format(idx / 1000000, len(lanes)))
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

def extract_one_lane(required_lane, lane_info, inputs, output_dir):
  outputs = {}
  buffer = [ '', '', '', '' ]

  for inp in inputs:
    logging.info('processing lane {} of {}'.format(required_lane, inp))
    sample = '{}_{}'.format( os.path.basename(inp).split('_')[0], os.path.basename(inp).split('_')[1] )
    tail = inp.split('_')[-1]
    written = 0
    with gzip.open(inp, 'r') as fh_in:
      for idx, line in enumerate(fh_in):
        buffer[idx % 4] = line
        if idx % 4 == 3: # 4 lines written
          # @ST-E00104:566:HLL3HCCXX:1:1101:28168:2047 1:N:0:NTTACTCG
          fields = buffer[0].strip('\n').split(':')
          lane = fields[3]
          if lane == required_lane:
            flowcell = fields[2]
            barcode = lane_info['barcode']
            filename = os.path.join(output_dir, '{sample}_{flowcell}_{barcode}_L00{lane}_{tail}'.format(sample=sample, flowcell=flowcell, barcode=barcode, lane=lane, tail=tail))
            if filename not in outputs:
              logging.info('writing to {}'.format(filename))
              outputs[filename] = gzip.open(filename, 'w')
            outputs[filename].write(buffer[0])
            outputs[filename].write(buffer[1])
            outputs[filename].write(buffer[2])
            outputs[filename].write(buffer[3])
            written += 4

        if idx % 10000000 == 0:
          logging.info('processed {} million lines, wrote {} looking for lane {} in {}'.format(idx / 1000000, written, required_lane, inp))
    logging.info('finished reading {}. wrote {} lines'.format(inp, written))

  logging.info('finished writing to {}'.format(outputs.keys()))

def extract_all_lanes(lanes, inputs, output_dir):
    threads = []
    for lane in lanes.keys():
        logging.info('starting thread for lane {}'.format(lane))
        lane_thread = threading.Thread(target=extract_one_lane, args=(str(lane), lanes[lane], inputs, output_dir))
        threads.append(lane_thread)
        lane_thread.start()

    for thread in threads:
        thread.join()

