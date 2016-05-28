# OSG-GEM
Pegasus workflow that utilizes Open Science Grid(OSG) resources to produce a Gene Expression Matrix(GEM) from DNA sequencing FASTQ files.

## Introduction

This workflow processes raw or compressed paired end FASTQ files to produce a two column matrix of normalized RNA molecule counts(FPKM).  An indexed reference genome along with gene model annotation files must be obtained prior to configuring
and running the workflow.
The following tasks are directed by the Pegasus workflow manager:

-Splitting input FASTQ files into files containing 20,000 sequences each.

-Trimming raw sequences with Trimmomatic

-Aligning reads to the reference genome using Hisat2 or Tophat2

-Merging BAM alignment files into a single sorted BAM file using Samtools

-Quantifying RNA transcript levels using StringTie or Cufflinks

It is suggested that the user become familiar with the documentation associated with the following software packages:

* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [Hisat2](https://ccb.jhu.edu/software/hisat2/manual.shtml)
* [Tophat2](https://ccb.jhu.edu/software/tophat/manual.shtml)
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
* [Samtools](http://www.htslib.org/doc/samtools.html)
* [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)
* [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/manual/)

## Open Science Grid / Execution Environment

The OSG-GEM workflow is designed to execute on the [Open Science Grid](https://www.opensciencegrid.org/) via the
[OSG Connect](https://osgconnect.net/) infrastructure. Access to the system can be requested on the 
[sign up page](https://osgconnect.net/signup). 

Once you have an account and have joined a project, a workflow specific ssh key has to be created. This key
is used for some of the data staging steps of the workflow. 

        $ mkdir -p ~/.ssh
        $ ssh-keygen -t rsa -b 2048 -f ~/.ssh/workflow
          (just hit enter when asked for a passphrase)
        $ cat ~/.ssh/workflow.pub >>~/.ssh/authorized_keys


## Example Workflow Setup

The worklow cloned from github contains an example config file as well as example input files from the 21st chromosome of Gencode Release 24 of the GRCh38 build of the human reference genome.  Two small FASTQ files containing
200,000 sequences from NCBI dataset SRR1825962 lie within the _Test_data_ directory of the workflow.  To run the test workflow, the user must copy the _osg-gem.config.template_ file:

        $ cp osg-gem.conf.template osg-gem.conf

The workflow, configured to run Hisat2 and Stringtie, can then be launched by running:

        $ ./submit

From here, the user may follow our documentation to modify the software options as well as point to their own input dataset.



## Pre-Workflow User Input

The user must provide indexed reference genome files as well as gene model annotation information prior to submitting the workflow.  The user must select a reference prefix($REF_PREFIX) that will be recognized by
Pegasus as well as by Hisat2 or Tophat2. In addition, information about splice sites or a reference transcriptome must be provided in order to guide accurate mapping of split input files.  Once the user has downloaded a reference genome fasta
file and gene annotation in GTF/GFF3 format, the following commands can be used to produce the necessary input files, using GRCh38 as an example $REF_PREFIX for Gencode Release 24 of the human reference genome:

### If the user would like to use Hisat2:

####Index the reference genome

       $ hisat2-build -f GRCh38.fa GRCh38

####Generate Tab delimited list of splice sites using gene model GTF file as input (Python DefaultDictionary Module necessary)

       $ python hisat2_extract_splice_sites.py GRCh38-gencode.v24.annotation.gtf > GRCh38.Splice_Sites.txt

### If the user would like to use Tophat2:

####Index the reference genome

       $ bowtie2-build GRCh38.fa GRCh38

####Generate and Index Reference Transcriptome

       $ tophat2 -G GRCh38.gencode.v24.annotation.gff3 --transcriptome-index=transcriptome_data/GRCh38 GRCh38
       $ tar czf GRCh38.transcriptome_data.tar.gz transcriptome_data/




## Workflow Configuration
Once the user has obtained necessary input files, the _osg-gem.config_ file must be appropriately modified and reference files must be placed into the _reference_ directory with appropriate filenames.

### Place Files in _reference_ directory

####If the user would like to use Hisat2, the following files must be present in the _reference_ directory:
$REF_PREFIX.fa

$REF_PREFIX.1.ht2 … $REF_PREFIX.N.ht2

$REF_PREFIX.Splice_Sites.txt

$REF_PREFIX.gff3

####If the user would like to use Tophat2, the following files must be present in the _reference_ directory:

$REF_PREFIX.fa

$REF_PREFIX.1.bt2 … $REF_PREFIX.N.bt2

$REF_PREFIX.rev.1.bt2

$REF_PREFIX.rev.2.bt2

$REF_PREFIX.transcriptome_data.tar.gz

$REF_PREFIX.gff3



### Modify _osg-gem.config_ file

#### Specify reference prefix that matches filenames in the _reference_  directory

[reference]

reference_prefix = $REF_PREFIX

#### Specify file path to forward and reverse FASTQ file for a given dataset($DATASET)

[inputs]

forward = /path/to/$DATASET_1.fastq.gz

reverse = /path/to/$DATASET_2.fastq.gz


#### Select software options

[config]

tophat2 = 'True' or 'False'

hisat2 = 'True' or 'False'

cufflinks = 'True' or 'False'

stringtie = 'True' or 'False'


#### Example _osg-gem.config_ file:

If a user cloned OSG-GEM into '/stash2/user/username/GEM_test', and placed input FASTQ files for dataset 'TEST' in '/stash2/user/username/Data'. To process this dataset using Hisat2 and StringTie with the GRCh38 build of the human reference genome, the osg-gem.config file would be modified as follows:

[reference]

reference_prefix = GRCh38

[inputs]

forward = /stash2/user/username/Data/TEST_1.fastq.gz

reverse = /stash2/user/username/Data/TEST_2.fastq.gz

[config]

tophat2 = False

hisat2 = True

cufflinks = False

stringtie = True


## Monitoring Workflow

Pegasus provides a set of commands to monitor workflow progress.  The path to the worklow files as well as commands to monitor the workflow will print to screen upon submitting the workflow.  For example:

        2016.05.26 23:31:03.859 CDT:   Your workflow has been started and is running in
        2016.05.26 23:31:03.869 CDT:     /stash2/user/username/workflows/osg-gem-x/workflow/osg-gem-x
        2016.05.26 23:31:03.880 CDT:   *** To monitor the workflow you can run ***
        
        2016.05.26 23:31:03.891 CDT:     pegasus-status -l /stash2/user/username/workflows/osg-gem-x/workflow/osg-gem-x
        
        2016.05.26 23:31:03.901 CDT:   *** To remove your workflow run ***
        
        2016.05.26 23:31:03.912 CDT:     pegasus-remove /stash2/user/username/workflows/osg-gem-x/workflow/osg-gem-x


Output will be transferred at the base of this directory upon completion.  For example:

        $ cd /stash2/user/username/workflows/osg-gem-x
        $ ls
        data  outputs  scratch  workflow

        $ cd outputs
        $ head -1 TEST-merged_counts.fpkm
        TranscriptID    FPKM


## User Modifications to Workflow

To customize OSG-GEM parameters, basic understanding of the directory structure of the workflow is necessary.  

### Workflow Directory Structure

####Test_data

This contains small FASTQ files for testing.  The user may place their own data in this directory or elsehwere on the OSG filesystems.  

####reference

Contains all reference genome and annotation files, as described previously.

####Tools

This directory contains job wrappers for each step of the workflow.  It is suggested that the user becomes familiar with the parameters set for each software to determine if they would like to make changes.  If the user would like to change software parameters, they may modify the commands in the files here.  Note that any changes to input filenames in the commands must match the files that are catalogued in the _task-files_ directory (explained below)

####task-files

This directory contains subdirectories for each job that utilizes specific files(eg., python script to parse StringTie output, fasta_adapters.txt file for trimmomatic).

Any files placed in these directories will be transferred to OSG compute nodes for the corresponding jobs.  For example, if the user would like to use a different fasta adapters file 'NewAdapters.txt' for read trimming for the hisat2 job, they would copy this file to the _hisat2_ directory.  Note that the job wrapper in the _tools_ directory must now be modified to match this filename.  

#### Base directory

The base directory of the worfklow contains the _osg-gem.config_ file, the _submit_ script, and a pegasus configuration file. 

The execution environment is catalogued in the _submit_ script, allowing the user to alter the resources requested by the workflow.  

For example:

        <site  handle="condorpool" arch="x86_64" os="LINUX">
                <profile namespace="pegasus" key="style" >condor</profile>
                <profile namespace="condor" key="universe" >vanilla</profile>
                <profile namespace="condor" key="requirements" >OSGVO_OS_STRING == "RHEL 6" &amp;&amp; HAS_MODULES == True &amp;&amp; HAS_SCP == True &amp;&amp; GLIDEIN_ResourceName != "Hyak" &amp;&amp; TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName1 &amp;&amp; TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName2 &amp;&amp; TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName3 &amp;&amp; TARGET.GLIDEIN_ResourceName =!= MY.MachineAttrGLIDEIN_ResourceName4</profile>
                <profile namespace="condor" key="request_memory" >5 GB</profile>
                <profile namespace="condor" key="request_disk" >30 GB</profile>
                <profile namespace="condor" key="+WantsStashCache" >True</profile>
        </site>

By default, the workflow is cloned with requests for at least 5 GB of memory and 30 GB of disk space on OSG compute nodes.  If the user is working with an organism with a large reference genome and finds that 5 GB is insufficient, they may change:

        <profile namespace="condor" key="request_memory" >5 GB</profile>
        
to request six gigabytes of RAM:

        <profile namespace="condor" key="request_memory" >6 GB</profile>
        

If the user finds that 5 gigabytes of RAM per job is unnecessary and would like to speed up queue times, they may change:

        <profile namespace="condor" key="request_memory" >5 GB</profile>

to request only 3 gigabytes of RAM per job:

        <profile namespace="condor" key="request_memory" >3 GB</profile>
        
###Interchanging Software

This workflow utilizes OASIS software modules that OSG compute nodes can access.  Job wrappers in this workflow load these modules to utilize specific versions of software.  For example, the following software modules are loaded for all _tophat_ jobs using the 'module load' command:

        module load tophat/2.1.1
        module load samtools/1.3.1
        module load bowtie/2.2.9
        module load java/8u25
        
If the user would like to plug in alternate software, or would like to use a different version of the available software, an osgconnect user support ticket may be submited to have their software of choice installed as an OASIS module.  

We have also found that precompiled software packages for linux x86_64 architecture have been stable on OSG compute nodes.  The user may utilize these software packages by adding a tar archive of the package to the appropriate task-files directory of the workflow.  This will then be transferred as input to the job, which can be unpacked and utilized for the user's task.  
        
        

        
        

        
        



