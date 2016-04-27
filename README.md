# OSG-GEM
Open Science Grid Workflow That Creates Gene Expression Matrices (GEMs) from SRA/FASTQ NGS Files 

## Introduction

(...)

## Configuration

The workflow is configured in a file named _osg-gem.conf_. A template is provided with the
workflow, so unless you already have a _osg-gem.conf_ file, copy the template:

    $ cp osg-gem.conf.template osg-gem.conf

Copy your reference genomes into the _references/_ directory. Note that the it is expected
that the file has already been indexed and that the index files are copied as well.

In _osg-gem.conf_, specify the prefix of the reference genome. Also specify the paths/URLs
for the input genome.

## Submitting / Monitoring

When the workflow has been configured, please start it by running:

    $ ./submit

In the output, a path for the submitted workflow will be provided together with an example
on how to monitor progress with _pegasus\_status_


