# hiplexpipe


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
