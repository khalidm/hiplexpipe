# hiplexpipe
***
## A bioinformatics pipeline for variant calling for [Hi-Plex](http://hiplex.org/) sequencing.

Author: Khalid Mahmood (kmahmood@unimelb.edu.au)

hiplexpipe is based on the [Ruffus](http://www.ruffus.org.uk/) library for writing bioinformatics pipelines. Its features include:

 * Job submission on a cluster using DRMAA (currently only tested with SLURM).
 * Job dependency calculation and checkpointing.
 * Pipeline can be displayed as a flowchart.
 * Re-running a pipeline will start from the most up-to-date stage. It will not redo previously completed tasks.

## License

See LICENSE.txt in source repository.

## Installation dependencies

#### External tools dependencies

`hiplexpipe` depends on the following programs and libraries:

 * [python](https://www.python.org/download/releases/2.7.5/) (version 2.7.5)
 * [java](https://java.com/en/download/) (version 1.8)
 * [DRMAA](http://www.drmaa.org/) for submitting jobs to the cluster (it uses the Python wrapper to do this).
   You need to install your own `libdrama.so` for your local job submission system. There are versions
   available for common schedulers such as Torque/PBS, [SLURM](http://apps.man.poznan.pl/trac/slurm-drmaa) and so on.
 <!-- * [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (version 0.10.1) -->
 * [SAMtools](http://www.htslib.org/doc/samtools-1.1.html) (version 1.3.1)
 * [bwa](http://bio-bwa.sourceforge.net/) for aligning reads to the reference genome (version 0.7.15)  
 * [GATK](https://software.broadinstitute.org/gatk/) for calling variants and genotyping (version 3.6)
 * [BEDTools](http://bedtools.readthedocs.io/) for calculating sequencing coverage statistics (version 2.26.0)

`hiplexpipe` assumes the tools above are installed by the users themselves.

#### Python dependencies

`hiplexpipe` depends on the following python libraries, tools and wrappers.

* Python 2.7.5
* [PyVCF](https://pypi.python.org/pypi/PyVCF)  
* [Biopython](https://pypi.python.org/pypi/biopython)
* [pybedtools](https://daler.github.io/pybedtools/)
* [cyvcf2](http://brentp.github.io/cyvcf2/)

We recommend using a python virtual environment. Following is an examples of how to setup a `hiplexpipe` virtual environment ready for analysis:

## Setup: New environment

The following instruction are based on the Melbourne Bioinformatics computing infrastructure.

    module load Python/2.7.12-GCC-4.9.3
    export DRMAA_LIBRARY_PATH=/usr/local/slurm_drmaa/1.0.7-GCC/lib/libdrmaa.so
    virtualenv --system-site-packages hiplexpipe
    source hiplexpipe/bin/activate
    pip install jupyter
    pip install plotly
    pip install pybedtools
    pip install -U https://github.com/khalidm/undr_rover/archive/master.zip
    pip install -U https://github.com/khalidm/hiplexpipe/archive/simple.zip
    pip install -U https://github.com/khalidm/offtarget/archive/master.zip
    mkdir references
    mkdir coverage
##### Test pipeline works with:

    hiplexpipe --config pipeline.config --use_threads --log_file pipeline.log --jobs 10 --verbose 3 --just_print

## Setup: New gene panel

#### Hi-Plex primer files

You should have two target interval files for every Hi-Plex experiment.

* rover.txt - this contains the amplicon regions and primer sequences.
* idt.txt - this file contains the primer sequences and their names matching the names in the above rover.txt file.

_Make sure heel sequences are removed from rover.txt file_

#### Additional interval files

Follow instructions below to prepare the intervals files for the pipeline. (We are working on a tool to automate this task).

##### Main rover bed file. (rover.bed)
Each interval in this bed file is the midpoint of each amplicon. This file is used to calculate alignment and coverage statistics.

    cut -f1,2,3,4,5 rover.txt > rover.bed
or

    awk ' BEGIN{FS="\t";OFS="\t"}; { print $1,int($2+($3-$2)/2),int($3-($3-$2)/2),$4,$5} ' rover.txt > rover.bed

##### Primer coordinates file. (primer.bedpe)
This file is used to clip primer sequences from the alignments.

    awk ' BEGIN{FS="\t";OFS="\t"}; { print $1,$7,$8,$1,$12,$11} ' rover.txt > primer.bedpe

##### Create intervals for GATK variant calling (gatk.bed)
This creates a bed file of intervals for GATK variant calling. Note this is different from rover.bed as this merges overlapping targets and mainly functions to provide a target for variant calling.

    cut -f1,2,3 rover.txt | bedtools slop -i - -b 10 -g hg19.genome | bedtools merge -i - > rover.gatk.bed

## New analysis  
#### Step 1. Load software requirements
    module load Python/2.7.12-GCC-4.9.3
    export DRMAA_LIBRARY_PATH=/usr/local/slurm_drmaa/1.0.7-GCC/lib/libdrmaa.so
    module load BEDTools/2.26.0-vlsci_intel-2015.08.25
    module load SAMtools
    module load VCFtools

#### Step 2. Preparing pipeline config files
I have created a template config file (pipeline_template.config) for all these analysis.

1. Create a new directory for the analysis
1. Make a copy of pipeline_template.config in the new analysis directory.
2. Make relevant changes to the new config file.
  - change pipeline_id
  - add fastq file paths
  - under the comment "hiplex files" - amend paths to files relevant to the design

#### Step 3: Create new screen and load modules
Log into snowy (HPC)

Run following commands:

    module load Python/2.7.12-GCC-4.9.3
    screen -S new_analysis
    module purge
    module load vlsci
    module load Python/2.7.12-GCC-4.9.3
    module load SAMtools
    module load VCFtools
    source hiplexpipe/bin/activate

#### Step 4: Run hiplexpipe

    hiplexpipe --config pipeline.config --use_threads --log_file pipeline.log --jobs 50 --verbose 2

## Generate statistics

### Alignment statistics

From within the virtualenv, run the following command:

    python alignment_stats.py > stats.txt

This will generate a table containing various alignment statistics for each sample.

### Heatmaps for alignment coverage

    jupyter nbconvert --ExecutePreprocessor.timeout=6000 --to html --execute coverage_analysis_main.ipynb

This will output `coverage_analysis_main.html` file.

### Offtarget
Generates statistics on which amplicons are mapping to incorrect regions of the genome, or not mapping at all.

Run for a few samples picked at random.

    offtarget --primers <rover.txt> --fastq1 <fastq_read_1> --fastq2 <fastq_read_2> --bam <sorted bam file> --log offtarget.log > output.txt
