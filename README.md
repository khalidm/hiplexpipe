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


## make interval files for the targets:

1. Main rover file. (rover.txt)
2. Primer sequence file (idt.csv)
3. Main rover bed file. (rover.bed)
    cut -f1,2,3,4,5 CRC_10g_23May16.final.rover.txt > CRC_10g_23May16.final.rover.bed
    or
    awk ' BEGIN{FS="\t";OFS="\t"}; { print $1,int($2+($3-$2)/2),int($3-($3-$2)/2),$4,$5} ' CRC_10g_23May16.final.rover.txt > CRC_10g_23May16.final.rover.bed
4. !Interval file. (rover.interval_list) - not required as the input bam is now clipped.
    java -jar picard.jar BedToIntervalList I=rover.bed SD=<hg19.dict> -O=rover.interval_list
5. Primer coordinates file. (primer.bedpe)
    awk ' BEGIN{FS="\t";OFS="\t"}; { print $1,$7,$8,$1,$12,$11} ' CRC_10g_23May16.final.rover.txt > CRC_10g_23May16.final.rover.primer.bedpe



## config file
# annotation tools
snpeff_path: /vlsci/VR0002/kmahmood/Programs/snpeff/snpEff4.3/snpEff/snpEff.jar
snpeff_conf: /vlsci/VR0002/kmahmood/reference/snpeff/4.1/snpEff_37.config
vep_path: /vlsci/VR0002/kmahmood/Programs/vep/ensembl-tools-release-87/scripts/variant_effect_predictor/
other_vep: --dir_cache /vlsci/VR0002/kmahmood/Programs/vep/barcoo/ensembl-tools-release-87/scripts/variant_effect_predictor/data/
anno: /vlsci/VR0002/kmahmood/Programs/gemini/gemini/gemini_data/gem.conf
annolua: /vlsci/VR0002/kmahmood/Programs/gemini/data/gemini_data/custom2.lua

#other tools
vt_path: /vlsci/VR0002/kmahmood/Programs/vt/vt/vt
bamclipper: /vlsci/VR0002/kmahmood/Programs/bamclipper/bamclipper-master/bamclipper.sh

# hiplex files
coord_file: /vlsci/VR0182/shared/km/hiplex/CRC_10g_23May16.final.rover.txt # primer coordinates - undr_rover
primer_file: /vlsci/VR0182/shared/km/hiplex/CRC_10g_23May16.final.idt.csv # primer sequences file - undr_rover
target_bed: /vlsci/VR0182/shared/km/hiplex/CRC_10g_23May16.final.rover.bed # bedtools multicov
interval_file: /vlsci/VR0182/shared/km/hiplex/CRC_10g_23May16.final.rover.interval_list # haplotypecaller gatk
primer_bedpe_file: /vlsci/UOM0040/shared/km/hiplex/GSN/CRC_10g_23May16.final.rover.primer.bedpe # bamclipper
#fragment_bed: goc.bed
#target_bed_merged: /vlsci/VR0182/shared/km/hiplex/CRC_merged.bed

# gatk
mills_hg19: /vlsci/VR0002/kmahmood/reference/Homo_sapiens/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
one_k_g_hg19_indels: /vlsci/VR0002/kmahmood/reference/Homo_sapiens/b37/1000G_phase1.indels.b37.vcf
one_k_g_snps: /vlsci/VR0002/kmahmood/reference/Homo_sapiens/b37/1000G_omni2.5.b37.vcf
one_k_g_highconf_snps: /vlsci/VR0002/kmahmood/reference/Homo_sapiens/b37/1000G_phase1.snps.high_confidence.b37.vcf
one_k_g_indels: /vlsci/VR0002/kmahmood/reference/Homo_sapiens/b37/1000G_phase1.indels.b37.vcf
hapmap: /vlsci/VR0002/kmahmood/reference/Homo_sapiens/b37/hapmap_3.3.b37.vcf
dbsnp_hg19: /vlsci/VR0002/kmahmood/reference/Homo_sapiens/b37/dbsnp_138.b37.vcf.gz

# undr_rover
proportionthresh: 0.3
absthresh: 10
maxvariants: 4 # max variants allowed per read

#coverdirs: /vlsci/VR0182/shared/km/hiplex/variants/undr_rover/coverdir/

# The Human Genome in FASTA format.
ref_grch37: /vlsci/UOM0040/shared/reference/human_g1k_v37_decoy.fasta
hrfile: /vlsci/VR0182/shared/hiplex_test/4genes/N7_plus_hg19_210817_sorted_slop.bed.gz # Homopolymer run of length
