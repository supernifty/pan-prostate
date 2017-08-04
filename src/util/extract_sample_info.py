#!/usr/bin/env python
'''
  given a list of sample IDs, write out the fields
  - donor ID
  - local normal sample ID
  - local tumour sample ID
  - normal uuid
  - tumour uuid
  - country (australia)
  for the train 1 list

  e.g.
  python extract_sample_info.py < train.in > train.out
'''

import collections
import logging
import os.path
import sys

VALIDATION_RESULT='/data/projects/punim0095/pan-prostate/out/{}.validation'

def find_uuid(sample, donor):
  #Donor_ID       Tissue_ID       is_normal (Yes/No,Y/N)  Sample_ID       relative_file_path      bam_state       Donor_UUID      Tissue_UUID     Sample_UUID     bam_md5sum
  filename = VALIDATION_RESULT.format(sample)
  if not os.path.isfile(filename):
    logging.warn('Failed to find validation file {} for sample {} (affects {})'.format(filename, sample, donor))
    return None

  for line in open(filename, 'r'):
    fields = line.strip('\n').split('\t')
    if fields[3] == sample:
      return fields[8]

  logging.warn('Failed to find sample {} in {}'.format(sample, filename))
  return None

def process(in_fh, out_fh, sample_fh):
  logging.info('reading sample metadata...')
  samples = {}

  for line in sample_fh:
    # Sample UUID,Patient UUID,Lab ID,tissue_id,is_normal,Status,Pre-aligned,Validated
    fields = line.strip('\n').split(',')
    # samples[uuid] = (patient uuid, is normal)
    samples[fields[0]] = (fields[1], fields[4])

  logging.info('reading sample metadata: done. {} samples read'.format(len(samples)))
  
  result = {}
  current_donor = None
  for line in sys.stdin:
    if line.startswith('#'):
      continue
    fields = line.strip('\n').split(',')
    donor = fields[0]
    if donor == '':
      donor = current_donor
    if donor not in result:
      result[donor] = {}

    current_donor = donor

    sample = fields[1]

    if sample not in samples:
      logging.warn('Provided sample {} not in sample metadata.'.format(sample))
    elif donor != samples[sample][0]:
      logging.warn('Provided donor {} for sample {} does not match metadata donor {}.'.format(donor, sample, samples[sample][0]))
    elif samples[sample][1] in result[donor]:
      #logging.warn('Added additional sample {} to donor {} is_normal {}'.format(sample, donor, samples[sample][1]))
      result[donor][samples[sample][1]].append(sample)
    else:
      result[donor][samples[sample][1]] = [sample]

  logging.info('finished reading input. {} donors'.format(len(result)))

  # look up sample details
  out_fh.write('{},{},{},{},{},{}\n'.format('Donor ID', 'Normal Sample ID', 'Tumour Sample ID', 'Normal UUID', 'Tumour UUID', 'Country'))
  donor_out = 0
  for donor in result:
    if len(result[donor]) != 2: # Y, N
      logging.warn('Donor {} has {} samples. Expected 2 samples.'.format(donor, len(result[donor])))
      continue

    # look up uuids
    if len(result[donor]['Y']) != 1:
      logging.warn('Donor {} has {} wb samples. Expected 1 sample.'.format(donor, len(result[donor]['Y'][0])))

    normal_uuid = find_uuid(result[donor]['Y'][0], donor) # assume only 1 wb

    # each tumor for donor
    for tumor_sample in result[donor]['N']:
      tumour_uuid = find_uuid(tumor_sample, donor)

      if normal_uuid is not None and tumour_uuid is not None:
        out_fh.write('{},{},{},{},{},{}\n'.format(donor, result[donor]['Y'][0], tumor_sample, normal_uuid, tumour_uuid, 'Australia'))
        donor_out += 1
      else:
        logging.warn('Failed to find at least one uuid for donor {} - normal {} uuid: {}. tumour {} uuid: {}'.format(donor, result[donor]['N'], normal_uuid, tumor_sample, tumour_uuid))

  logging.info('Wrote {} samples from {} donors'.format(donor_out, len(result)))

if __name__ == '__main__':
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  if len(sys.argv) != 2:
    sys.stderr.write('Usage: basename sample-metadata < list\n')
  else:
    process(sys.stdin, sys.stdout, open(sys.argv[1], 'r'))
