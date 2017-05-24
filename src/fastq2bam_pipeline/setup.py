#!/usr/bin/env python

from setuptools import setup

setup(
    name='fastq2bam_pipeline',
    version='0.0.6',
    author='Peter Georgeson',
    author_email='peter.georgeson@unimelb.edu.au',
    packages=['src'],
    entry_points={
        'console_scripts': ['fastq2bam_pipeline = src.main:main']
    },
    url='https://github.com/supernifty/fastq2bam_pipeline',
    license='LICENSE',
    description='fastq2bam_pipeline is a demonstration pipeline based on Ruffus', 
    long_description=open('README.md').read(),
    install_requires=[
        "ruffus == 2.6.3",
        "pipeline_base == 1.0.0",
        "plotly == 2.0.7"
    ],
)
