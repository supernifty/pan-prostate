
This is the pan prostate pipeline installed for spartan.

# Folders

The code and working directory is /data/projects/punim0095/pan-prostate

* deploy: contains the deployed ruffus pipeline
* img: contains generated singularity images
* in: input data
* out: generated data
* reference: common reference data
* scripts: helper scripts used for testing
* src: source code for ruffus pipeline, other helper scripts and templates
* tmp: temporary files created during analysis
* tools: third party tools used in the pipeline

The data is mounted at:
* /data/punim0261/data01 - generated results
* /data/punim0261/data02 - input fastq files

The account with access is:
* punim0261

# Installation

## Installing dependencies

* bedtools must be on the path
* biobambam2: wget https://github.com/gt1/biobambam2/releases/download/2.0.65-release-20161130121735/biobambam2-2.0.65-release-20161130121735-x86_64-etch-linux-gnu.tar.gz in tools
* fastqc must be on the path
* picard: wget https://github.com/broadinstitute/picard/releases/download/2.8.2/picard-2.8.2.jar in tools

## Installing the pipeline

```
cd deploy
virtualenv venv
source ./venv/bin/activate
pip install -U pip
pip install ../tools/ruffus
pip install ../tools/pipeline_base
pip install -U ../src/fastq2bam_pipeline/
```

```
cp ../src/util/pipeline.config .
```

# Inputs

* ROOT/cfg/sample-metadata.csv contains sample metadata of the form: 
```
Sample UUID,Patient UUID,Lab ID,tissue_id,is_normal,Status,Pre-aligned,Validated
CMHS1,CMHP1,299bT,BT,N,,,
```

## Symlinks

* pipeline.config expects the fastq inputs to be in the out directory.
* generate symlinks to the actual input files.

```
./src/util/make_symlinks.sh
```

This generates a *sources* file containing all fastq files, then uses ./cfg/uuid_map to convert these files to UUIDs.

A *links* file is created, which can be executed in the out directory to generate all symlinks.
```
cd out
bash ../links-YYYYMMDD.sh

# Running the pipeline

* root is set to ROOT='/data/projects/punim0095/pan-prostate' in stages.py and fastq2bam.py
* copy src/util/pipeline.config deploy/
* edit deploy/pipeline.config to include all the samples. Use links-*.sh to get the list.

```
cd deploy
. ../src/util/env.sh
fastq2bam_pipeline --config pipeline.config --verbose 3 --jobs 10 --use_threads --log_file ./test`date +%Y%m%d`.log
```

# Components

* img/cgpqc.img: runs a validation script
* img/cgpmap.img: aligns the pre-aligned bam

# TODO

# Notes

## Converting the docker container to run bwa alignment to a singularity container

Steps to building the singularity image:
Locally:
```
git clone https://github.com/cancerit/dockstore-cgpmap
git checkout master
docker build -t cgpmap .
docker run -it cgpmap:latest bash
docker ps
docker export adoring_blackwell > cgpmap.export.tar
singularity create cgpmap.img
singularity import cgpmap.img cgpmap.export.tar
```

Running the singularity image:
```
module load Singularity
```

A docker command from nectar:
docker run -i --volume=/mnt/vicnode_nfs/jobs/fastq2bam/./datastore/launcher-9b14e477-c989-4ee8-9685-c83057055ad5/inputs/35bba1e6-0503-441b-9805-ea787d38708b/bwa_idx_GRCh37d5.tar.gz:/var/lib/cwl/stg0185577a-0938-4c49-9a54-d68e7869850b/bwa_idx_GRCh37d5.tar.gz:ro --volume=/mnt/vicnode_nfs/jobs/fastq2bam/./datastore/launcher-9b14e477-c989-4ee8-9685-c83057055ad5/inputs/93494ef6-d851-405a-9ccc-624aa721385b/CMHS231.bam:/var/lib/cwl/stg495d3698-e127-4365-989b-16794e2a0060/CMHS231.bam:ro --volume=/mnt/vicnode_nfs/jobs/fastq2bam/./datastore/launcher-9b14e477-c989-4ee8-9685-c83057055ad5/inputs/06404b25-e971-4a5b-8ac3-fd7d6aa6ddb6/core_ref_GRCh37d5.tar.gz:/var/lib/cwl/stgcf238f68-053e-46e0-854a-8dbdc0f7cc00/core_ref_GRCh37d5.tar.gz:ro --volume=/mnt/vicnode_nfs/dockstore-tmp/ba232dd7-037c-48fa-9dbd-55d2061ddce6-8b4e45d8-3e6b-47cf-a7cf-e41d198938ce/tmpWk0CbJ:/var/spool/cwl:rw --volume=/mnt/vicnode_nfs/jobs/fastq2bam/datastore/launcher-9b14e477-c989-4ee8-9685-c83057055ad5/working_kf81v:/tmp:rw --workdir=/var/spool/cwl --read-only=true --user=1000 --rm --env=TMPDIR=/tmp --env=HOME=/var/spool/cwl quay.io/wtsicgp/dockstore-cgpmap:2.0.0 /opt/wtsi-cgp/bin/ds-wrapper.pl -reference /var/lib/cwl/stgcf238f68-053e-46e0-854a-8dbdc0f7cc00/core_ref_GRCh37d5.tar.gz -bwa_idx /var/lib/cwl/stg0185577a-0938-4c49-9a54-d68e7869850b/bwa_idx_GRCh37d5.tar.gz -sample ba232dd7-037c-48fa-9dbd-55d2061ddce6 -scramble  -bwa  -Y -K 100000000 -t 4 /var/lib/cwl/stg495d3698-e127-4365-989b-16794e2a0060/CMHS231.bam


## ruffus pipeline to run pre-alignment and alignment