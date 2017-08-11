#!/usr/bin/env python

import datetime
import glob
import sys
import dateutil.parser


for f in glob.glob('CMHS[0-9]*.bam.log.err'):
  if '.mapped.' in f:
    fields = f.split('.')
    sample = fields[0]
    job = 'Align'
    started = None
    finished = None
    for line in open(f, 'r'):
      if 'Start workflow: ' in line:
        date = line.split('Start workflow: ')[1]
        #sys.stderr.write('{}: {}'.format(f, date))
        started = dateutil.parser.parse(date)
      if 'Workflow end: ' in line:
        date = line.split('Workflow end: ')[1]
        #sys.stderr.write('{}: {}'.format(f, date))
        finished = dateutil.parser.parse(date)
    if started is not None and finished is not None:
      sys.stdout.write('{},{},{},{},{}\n'.format(sample, job, started.strftime('%Y-%m-%d %H:%M:%S'), finished.strftime('%Y-%m-%d %H:%M:%S'), (finished - started).total_seconds()))
  else:
    fields = f.split('.')
    sample = fields[0]
    job = 'prealign'
    started = None
    finished = None
    try:
      for line in open(f, 'r'):
        if not started:
          date = ' '.join(line.strip().split(' ')[0:2])
          #sys.stderr.write('{}: {}'.format(f, date))
          started = dateutil.parser.parse(date)
      
      date = ' '.join(line.strip().split(' ')[0:2])
      #sys.stderr.write('{}: {}'.format(f, date))
      finished = dateutil.parser.parse(date)
      if started is not None and finished is not None:
        sys.stdout.write('{},{},{},{},{}\n'.format(sample, job, started.strftime('%Y-%m-%d %H:%M:%S'), finished.strftime('%Y-%m-%d %H:%M:%S'), (finished - started).total_seconds()))
    except ValueError:
      sys.stderr.write('skipped {}\n'.format(f))

# wgs
for f in glob.glob('CMHS[0-9]*.wgs.*.log.out'):
  fields = f.split('.')
  sample = fields[0]
  job = fields[2]
  started = None
  finished = None
  for line in open(f, 'r'):
    # look for starting
    if 'Starting {}'.format(job) in line:
      date = line.split('Starting {}'.format(job))[1].split(' and ')[0]
      started = dateutil.parser.parse(date)
    if 'Finished {}'.format(job) in line:
      date = line.split('Finished {}'.format(job))[1].split(' with ')[0]
      finished = dateutil.parser.parse(date)
  if started is not None and finished is not None:
    sys.stdout.write('{},{},{},{},{}\n'.format(sample, job, started.strftime('%Y-%m-%d %H:%M:%S'), finished.strftime('%Y-%m-%d %H:%M:%S'), (finished - started).total_seconds()))

# delly
for f in glob.glob('CMHS[0-9]*.delly.log.err'):
  fields = f.split('.')
  sample = fields[0]
  job = 'delly'
  started = None
  finished = None
  for line in open(f, 'r'):
    # look for dates of the form INFO : 2017-08-10 22:06:55,910 : Poll job KBPKQ941GT
    if started is None and line.startswith('INFO'):
      date = line.split('INFO : '.format(job))[1].split(',')[0]
      started = dateutil.parser.parse(date)
    # last line is INFO : 2017-08-10 22:24:38,072 : Archive bam_qc directory /tmp/workdir/bam_qc
    if finished is None and 'Archive bam_qc directory' in line:
      date = line.split('INFO : '.format(job))[1].split(',')[0]
      finished = dateutil.parser.parse(date)
  if started is not None and finished is not None:
    sys.stdout.write('{},{},{},{},{}\n'.format(sample, job, started.strftime('%Y-%m-%d %H:%M:%S'), finished.strftime('%Y-%m-%d %H:%M:%S'), (finished - started).total_seconds()))
