'''
Individual stages of the pipeline implemented as functions from
input files to output files.

The run_stage function knows everything about submitting jobs and, given
the state parameter, has full access to the state of the pipeline, such
as config, options, DRMAA and the logger.
'''

from pipeline_base.utils import safe_make_dir
from pipeline_base.runner import run_stage

import os
import re
import sys
import uuid

ROOT='/data/projects/punim0095/pan-prostate'
OUT='/data/punim0261/data01/out'
IN='/data/punim0261/data02/original_data'
REFERENCE='/data/projects/punim0095/pan-prostate/reference'
TMP='/data/projects/punim0095/pan-prostate/tmp'

class Stages(object):
    def __init__(self, state):
        self.state = state

    def get_stage_options(self, stage, *options):
        return self.state.config.get_stage_options(stage, *options)

    def get_options(self, *options):
        return self.state.config.get_options(*options)

    def original_files(self, output):
        '''Original files'''
        pass

    def fastqc(self, fastq_in, dir_out):
        '''Quality check fastq file using fastqc'''
        safe_make_dir(dir_out)
        command = 'fastqc --extract -o {dir} {fastq}'.format(dir=dir_out, fastq=fastq_in)
        run_stage(self.state, 'fastqc', command)

    def fastq2bam(self, inputs, bam_out, sample):
        '''
          Convert fastq to a prealigned bam. 
          stages:
          1 infer lanes and indexes
          2 split into lanes
          3 fastq2bam
          4 merge
        '''

        # input filenames
        fastq_read1_in, fastq_read2_in = inputs
        output_dir = os.path.dirname(bam_out)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        log_out = os.path.join(output_dir, '{}.log.out'.format(bam_out))
        log_err = os.path.join(output_dir, '{}.log.err'.format(bam_out))

        command = "python {}/src/util/fastq2bam.py --r1 {} --r2 {} --output_dir {} --bam {} 1>{} 2>{}".format(ROOT, fastq_read1_in, fastq_read2_in, output_dir, bam_out, log_out, log_err)
        run_stage(self.state, 'fastq2bam', command)

    def validate_prealigned_bam(self, input, validation_out):
        '''
            run validation script
            @input: the pre-aligned bam
            @validation_out: tsv file with validation details
        '''
        prefix = re.sub('.bam$', '', input)
        sample = re.sub('.bam$', '', os.path.basename(input))

        validation_in = '{}.validation_src'.format(prefix)
        # read in additional metadata
        found = False
        for line in open("{}/cfg/sample-metadata.csv".format(ROOT), 'r'):
            # Sample UUID,Patient UUID,Lab ID,tissue_id,is_normal
            fields = line.strip('\n').split(',')
            if fields[0] == sample:
                donor_id = fields[1]
                tissue_id = fields[3]
                is_normal = fields[4]
                found = True
                break

        if not found:
            raise Exception("Sample '{}' not found in metadata file".format(sample))

        # generate input to the validation script
        with open(validation_in, 'w') as validation_src:
            validation_src.write('#Donor_ID\tTissue_ID\tis_normal (Yes/No,Y/N)\tSample_ID\trelative_file_path\n')
            validation_src.write('{donor_id}\t{tissue_id}\t{is_normal}\t{sample_id}\t{sample}.bam\n'.format(
                donor_id=donor_id, 
                tissue_id=tissue_id, 
                is_normal=is_normal, 
                sample_id=sample, 
                sample=sample))

        # make our own align script
        tmp_id = '{}-{}'.format(sample, str(uuid.uuid4()))
        tmp_dir = '{root}/tmp/{tmp_id}'.format(root=ROOT, tmp_id=tmp_id)
        safe_make_dir(tmp_dir)
        with open('{tmp_dir}/validate.sh'.format(tmp_dir=tmp_dir), 'w') as align_fh:
            for line in open('{root}/src/util/validate.sh.template'.format(root=ROOT), 'r'):
                new_line = re.sub('TMP_ID', tmp_id, line)
                new_line = re.sub('SAMPLE', sample, new_line)
                align_fh.write(new_line)

        # run the validation script and generate output
        #command = ". {root}/src/util/profile; validate_sample_meta.pl -in {validation_in} -out {validation_out} -f tsv 1>{prefix}.validation.out 2>{prefix}.validation.err".format(root=ROOT, validation_in=validation_in, validation_out=validation_out, prefix=prefix)
        command = 'singularity exec -i --bind {in_dir}:/mnt/in,{out}:/mnt/out,{reference}:/mnt/reference,{tmp}:/mnt/tmp --workdir {tmp_dir} --contain {root}/img/cgpqc.img bash /mnt/tmp/{tmp_id}/validate.sh'.format(root=ROOT, in_dir=IN, out=OUT, reference=REFERENCE, tmp=TMP, tmp_dir=tmp_dir, tmp_id=tmp_id)
        run_stage(self.state, 'validate_prealigned_bam', command)

        # check that it worked - but run_stage doesn't block
        #lines = open(validation_out, 'r').readlines()
        #if len(lines) != 2:
        #    raise Exception('{} contained {} lines. Expected 2 lines.'.format(validation_out, len(lines)))
        #fields = lines[1].strip('\n').split('\t')
        #if len(fields) != 10:
        #    raise Exception('{} contained {} fields. Expected 10.'.format(validation_out, len(fields)))

    def align(self, inputs, bam_out):
        '''
          run the alignment dockstore image
          @input: the pre-aligned bam
          @bam_out: aligned bam
        '''
        # generate dockstore file as sample.dockstore
        validation, bam = inputs
        prefix = re.sub('.bam$', '', bam) # full path without .bam
        sample_filename = prefix.split('/')[-1] # e.g. CMHS1
        dockstore_out = re.sub('.bam$', '.dockstore', bam)

        # determine sample from validation file
        for line in open(validation, 'r'):
            if line.startswith('#'):
                continue
            fields = line.strip('\n').split('\t')
            sample = fields[8]

        if input == dockstore_out:
            raise Exception("Unexpected input file {}".format(bam))

        #log_out = '{}.log.out'.format(bam_out)
        #log_err = '{}.log.err'.format(bam_out)

        # make our own align script
        tmp_id = '{}-{}'.format(sample, str(uuid.uuid4()))
        tmp_dir = '{root}/tmp/{tmp_id}'.format(root=ROOT, tmp_id=tmp_id)
        safe_make_dir(tmp_dir)
        with open('{tmp_dir}/align.sh'.format(tmp_dir=tmp_dir), 'w') as align_fh:
            for line in open('{root}/src/util/align.sh.template'.format(root=ROOT), 'r'):
                new_line = re.sub('TMP_ID', tmp_id, line)
                new_line = re.sub('SAMPLE_FILENAME', sample_filename, new_line)
                new_line = re.sub('SAMPLE_ID', sample, new_line)
                align_fh.write(new_line)

        command = 'singularity exec -i --bind {in_dir}:/mnt/in,{out}:/mnt/out,{reference}:/mnt/reference,{tmp}:/mnt/tmp --workdir {tmp_dir} --contain {root}/img/cgpmap.img bash /mnt/tmp/{tmp_id}/align.sh 1>{prefix}.mapped.log.out 2>{prefix}.mapped.log.err && rm -rf "{tmp_dir}"'.format(root=ROOT, in_dir=IN, out=OUT, reference=REFERENCE, tmp=TMP, tmp_dir=tmp_dir, tmp_id=tmp_id, prefix=prefix)
        run_stage(self.state, 'align', command)

    def align_stats_bedtools(self, inputs, stats_out):
        '''
          generate coverage stats from bam
        '''
        mapped_bam = inputs
        command = 'bedtools genomecov -ibam {mapped_bam} | python {root}/src/util/coverage_histogram.py {stats_out}.histogram.html 1>{stats_out} 2>{stats_out}.err'.format(root=ROOT, mapped_bam=mapped_bam, stats_out=stats_out)
        run_stage(self.state, 'align_stats_bedtools', command)

    def align_stats_picard(self, inputs, stats_out):
        '''
          generate coverage stats from bam
        '''
        mapped_bam = inputs
        command = 'java -jar {root}/tools/picard-2.8.2.jar CollectRawWgsMetrics INPUT={input} OUTPUT={output} REFERENCE_SEQUENCE={root}/reference/core_ref_GRCh37d5/genome.fa INCLUDE_BQ_HISTOGRAM=true 1>{output}.log.out 2>{output}.log.err'.format(root=ROOT, input=mapped_bam, output=stats_out)
        run_stage(self.state, 'align_stats_picard', command)

