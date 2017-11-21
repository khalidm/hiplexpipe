# hiplexpipe

## A bioinformatics pipeline for variant calling for [Hi-Plex] (http://hiplex.org/)

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

We recommend using a python virtual environment. Following is an examples of how to setup a `hiplexpipe` virtual environment ready for analysis:

## Requirements
  * Python 2.7
  * [PyVCF](https://pypi.python.org/pypi/PyVCF)  
  * [Biopython](https://pypi.python.org/pypi/biopython)
  * [pybedtools](https://daler.github.io/pybedtools/)
  * [cyvcf2](http://brentp.github.io/cyvcf2/)

## Usage example on Melbourne Bioinformatics (VLSCI)
```
module load Python/2.7.10-vlsci_intel-2015.08.25
export DRMAA_LIBRARY_PATH=/usr/local/slurm_drmaa/1.0.7-GCC/lib/libdrmaa.so
virtualenv --system-site-packages venv
source venv/bin/activate
pip install -U https://github.com/khalidm/undr_rover/archive/master.zip
pip install -U https://github.com/khalidm/hiplexpipe/archive/benchmark.zip
hiplexpipe --config pipeline.config --use_threads --log_file pipeline.log --jobs 10 --verbose 3 --just_print
```
