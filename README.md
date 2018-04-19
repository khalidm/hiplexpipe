# hiplexpipe

## A bioinformatics pipeline for variant calling for [Hi-Plex](http://hiplex.org/) sequencing.

Author: Khalid Mahmood (kmahmood@unimelb.edu.au)

hiplexpipe is based on the [Ruffus](http://www.ruffus.org.uk/) library for writing bioinformatics pipelines. Its features include:

 * Job submission on a cluster using DRMAA (currently only tested with SLURM).
 * Job dependency calculation and checkpointing.
 * Pipeline can be displayed as a flowchart.
 * Re-running a pipeline will start from the most up-to-date stage. It will not redo previously completed tasks.

## License

See LICENSE.txt in source repository.

## Installation

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

## Getting started

The following instruction are based on the Melbourne Bioinformatics computing infrastructure.

### Step 1. Setup Python
    `module load Python/2.7.12-GCC-4.9.3`
### Step 2. Load software requirements
    `export DRMAA_LIBRARY_PATH=/usr/local/slurm_drmaa/1.0.7-GCC/lib/libdrmaa.so`
    Other Software dependancies    
    `module load BEDTools/2.26.0-vlsci_intel-2015.08.25`
    `bgzip, samtools, picard`

### Step 3. Installing hiplexpipe.
    ```
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
    ```

    Test:

    ```
    hiplexpipe --config pipeline.config --use_threads --log_file pipeline.log --jobs 10 --verbose 3
     --just_print
     ```

### Step 4. Preparing pipeline config files

I have created a template config file (pipeline_template.config) for all these analysis.

    For a new analysis - create a new directory ( I have created 4gp_analysis for the 4gene panel).
    Make a copy of pipeline_template.config in the new analysis directory. Make relevant changes to the new config file.
        - change pipeline_id and add fastq file paths
        - under the comment "hiplex files" - amend paths to files relevant to the design
        - see 4gp_analysis/pipeline.config as an example

### Step 5. Run the pipeline - this is an example for the '4gp' analysis

a. log in to snowy
b. module load Python/2.7.12-GCC-4.9.3
c. screen -S new_analysis
d. module purge
e. module load vlsci
f. module load Python/2.7.12-GCC-4.9.3
g. module load SAMtools
h. module load VCFtools
i. source hiplexpipe/bin/activate (assuming above Step 3 is performed)
j. command: `hiplexpipe --config pipeline.config --use_threads --log_file pipeline.log --jobs 50 --verbose 2`

### Step 6. Generate alignment statistics

a. alignment statistics: from within the virtualenv run the follwing command:
    `python alignment_stats.py > stats.txt`

### Step 7. Generate heatmaps for alignment coverage
    `jupyter nbconvert --ExecutePreprocessor.timeout=6000 --to html --execute coverage_analysis_main.ipynb`
    This will output 'coverage_analysis_main.html` file.

## Preparing target region files

You should have two target interval files for every Hi-Plex experiment.

rover.txt - this contains the amplicon regions and primer sequences.
idt.txt - this file contains the primer sequences and their names matching the names in the above rover.txt file.

Follow instructions below to prepare the intervals files for the pipeline. (We are working on a tool to automate this task).

* Main rover bed file. (rover.bed) This file is used to calculate alignment and coverage statistics. cut -f1,2,3,4,5 rover.txt > rover.bed or awk ' BEGIN{FS="\t";OFS="\t"}; { print $1,int($2+($3-$2)/2),int($3-($3-$2)/2),$4,$5} ' rover.txt > rover.bed
* Primer coordinates file. (primer.bedpe) This file is used to clip primer sequences from the alignments. awk ' BEGIN{FS="\t";OFS="\t"}; { print $1,$7,$8,$1,$12,$11} ' rover.txt > primer.bedpe
* Create intervals for GATK variant calling
`cut -f1,2,3 rover.txt | bedtools slop -i - -b 10 -g hg19.genome | bedtools merge -i - > rover.gatk.bed`
