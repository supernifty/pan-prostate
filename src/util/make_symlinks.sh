#!/bin/bash

IN=/data/projects/punim0095/pan-prostate/in/
UUIDS=/data/projects/punim0095/pan-prostate/cfg/uuids_map

# target files
SOURCES=./sources-$(date +%Y%m%d)
LINKS=./links-$(date +%Y%m%d).sh

# generates a script to generate links,
# then generates the links

# generate a list of inputs
# all fastq inputs
echo "finding all fastq files"
find -L $IN -name \*.fastq.gz > $SOURCES

# priority samples
#find /mnt/vicnode_nfs/original_data/1611KHX-0054_hdd1 -name \*.fastq.gz > sources-170321
#find /mnt/vicnode_nfs/original_data/1611KHX-0054_hdd2 -name \*.fastq.gz >> sources-170321

# make a script to generate links
echo "converting to UUIDs"
python src/util/convert_to_uuid.py $SOURCES < $UUIDS > $LINKS

# generate the links
echo "create links by running $LINKS in target dir"
#/bin/bash $LINKS

echo "done"
