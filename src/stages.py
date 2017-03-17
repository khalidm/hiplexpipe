'''
Individual stages of the pipeline implemented as functions from
input files to output files.

The run_stage function knows everything about submitting jobs and, given
the state parameter, has full access to the state of the pipeline, such
as config, options, DRMAA and the logger.
'''

from utils import safe_make_dir
from runner import run_stage
import os

# PICARD_JAR = '$PICARD_HOME/lib/picard-1.69.jar'
PICARD_JAR = '/vlsci/VR0002/kmahmood/Programs/picard/picard-tools-2.0.1/picard.jar'
SNPEFF_JAR = '/usr/local/easybuild/software/snpEff/4.1d-Java-1.7.0_80/snpEff.jar'

GATK_JAR = '$GATK_HOME/GenomeAnalysisTK.jar'

def java_command(jar_path, mem_in_gb, command_args):
    '''Build a string for running a java command'''
    # Bit of room between Java's max heap memory and what was requested.
    # Allows for other Java memory usage, such as stack.
    java_mem = mem_in_gb - 2
    return 'java -Xmx{mem}g -jar {jar_path} {command_args}'.format(
        jar_path=jar_path, mem=java_mem, command_args=command_args)


def run_java(state, stage, jar_path, mem, args):
    command = java_command(jar_path, mem, args)
    run_stage(state, stage, command)


class Stages(object):
    def __init__(self, state):
        self.state = state
        self.reference = self.get_options('ref_grch37')
        # self.dbsnp_hg19 = self.get_options('dbsnp_hg19')
        # self.mills_hg19 = self.get_options('mills_hg19')
        # self.one_k_g_snps = self.get_options('one_k_g_snps')
        # self.one_k_g_indels = self.get_options('one_k_g_indels')
        # self.one_k_g_highconf_snps = self.get_options('one_k_g_highconf_snps')
        # self.hapmap = self.get_options('hapmap')
        # self.interval_hg19 = self.get_options('exome_bed_hg19')
        # self.CEU_mergeGvcf = self.get_options('CEU_mergeGvcf')
        self.snpeff_conf = self.get_options('snpeff_conf')
        self.vep_path = self.get_options('vep_path')
        self.vt_path = self.get_options('vt_path')
        self.coord_file = self.get_options('coord_file')
        self.primer_file = self.get_options('primer_file')
        self.proportionthresh = self.get_options('proportionthresh')
        self.absthresh = self.get_options('absthresh')
        self.coverdir = self.get_options('coverdirs')
        # self.GBR_mergeGvcf = self.get_options('GBR_mergeGvcf')
        # self.FIN_mergeGvcf = self.get_options('FIN_mergeGvcf')

    def run_picard(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
        return run_java(self.state, stage, PICARD_JAR, mem, args)

    def run_snpeff(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
        return run_java(self.state, stage, SNPEFF_JAR, mem, args)

    def run_gatk(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
        return run_java(self.state, stage, GATK_JAR, mem, args)

    def get_stage_options(self, stage, *options):
        return self.state.config.get_stage_options(stage, *options)

    def get_options(self, *options):
        return self.state.config.get_options(*options)

    def original_fastqs(self, output):
        '''Original fastq files'''
        # print output
        pass

    def align_bwa(self, inputs, bam_out, sample_id, read_id, lane, lib):
        # def align_bwa(self, inputs, bam_out, sample_id):
        '''Align the paired end fastq files to the reference genome using bwa'''
        fastq_read1_in, fastq_read2_in = inputs
        cores = self.get_stage_options('align_bwa', 'cores')
        safe_make_dir('alignments/{sample}'.format(sample=sample_id))
        read_group = '"@RG\\tID:{readid}\\tSM:{sample}_{read_id}\\tPU:lib1\\tLN:{lane}\\tPL:Illumina"' \
            .format(readid=read_id, lib=lib, lane=lane, sample=sample_id)
        command = 'bwa mem -M -t {cores} -R {read_group} {reference} {fastq_read1} {fastq_read2} ' \
                  '| samtools view -b -h -o {bam} -' \
                  .format(cores=cores,
                          read_group=read_group,
                          fastq_read1=fastq_read1_in,
                          fastq_read2=fastq_read2_in,
                          reference=self.reference,
                          bam=bam_out)
        run_stage(self.state, 'align_bwa', command)

    def apply_undr_rover(self, inputs, vcf_output, sample_id, read_id, lane, lib):
        # def align_bwa(self, inputs, bam_out, sample_id):
        '''Align the paired end fastq files to the reference genome using bwa'''
        fastq_read1_in, fastq_read2_in = inputs
        cores = self.get_stage_options('align_bwa', 'cores')
        safe_make_dir('variants/undr_rover')
        read_group = '"@RG\\tID:{readid}\\tSM:{sample}_{read_id}\\tPU:lib1\\tLN:{lane}\\tPL:Illumina"' \
            .format(readid=read_id, lib=lib, lane=lane, sample=sample_id)
        command = 'undr_rover --primer_coords {coord_file} ' \
                  '--primer_sequences {primer_file} ' \
                  '--reference {reference} ' \
                  '--out {vcf_output} ' \
                  '--genotype ' \
                  '--coverdir {coverdir}' \
                  '--proportionthresh {propmoortionthresh} ' \
                  '--absthresh {absthresh} ' \
                  '{fastq_read1} {fastq_read2}'.format(
                        coord_file=self.coord_file, primer_file=self.primer_file,
                        reference=self.reference,
                        vcf_output=vcf_output,
                        proportionthresh=self.proportionthresh,
                        absthresh=self.absthresh,
                        coverdir=self.coverdir,
                        fastq_read1=fastq_read1_in,
                        fastq_read2=fastq_read2_in)
        run_stage(self.state, 'apply_undr_rover', command)

    def sort_bam_picard(self, bam_in, sorted_bam_out):
        '''Sort the BAM file using Picard'''
        picard_args = 'SortSam INPUT={bam_in} OUTPUT={sorted_bam_out} ' \
                      'VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate ' \
                      'MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=True'.format(
                          bam_in=bam_in, sorted_bam_out=sorted_bam_out)
        self.run_picard('sort_bam_picard', picard_args)

    #samtools
    def apply_samtools_mpileup(self, bam_in, mpileup_out_bcf):
        '''Samtools mpileup'''
        # bam_in = bam_in
        bams = ' '.join([bam for bam in bam_in])
        safe_make_dir('variants')
        command = 'samtools mpileup -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -go {mpileup_out_bcf} ' \
                  '-f {reference} {bams}'.format(
                          mpileup_out_bcf=mpileup_out_bcf,reference=self.reference,bams=bams)
        run_stage(self.state, 'apply_samtools_mpileup', command)

    #bcftools
    def apply_bcftools(self, mpileup_in, vcf_out):
        '''Bcftools call variants'''
        mpileup_in = mpileup_in
        # mpileup_in = ' '.join([vcf for vcf in vcf_files_in])
        command = 'bcftools call -vmO v -o {vcf_out} {mpileup_in}'.format(
                          vcf_out=vcf_out,mpileup_in=mpileup_in)
        run_stage(self.state, 'apply_bcftools', command)

    def apply_vt(self, inputs, vcf_out):
        '''Apply NORM'''
        vcf_in = inputs
        cores = self.get_stage_options('apply_vt', 'cores')
        vt_command = "{vt_path} decompose -s {vcf_in} - | {vt_path2} normalize -r {reference} - | " \
                    "{vt_path3} uniq - -o {vcf_out}".format(
                    vt_path=self.vt_path, vcf_in=vcf_in, vt_path2=self.vt_path, reference=self.reference,
                    vt_path3=self.vt_path, vcf_out=vcf_out)
        run_stage(self.state, 'apply_vt', vt_command)

    def apply_vep(self, inputs, vcf_out):
        '''Apply VEP'''
        vcf_in = inputs
        cores = self.get_stage_options('apply_vep', 'cores')
        vep_command = "{vep_path}/variant_effect_predictor.pl --cache -i {vcf_in} --format vcf -o {vcf_vep} --force_overwrite --vcf " \
                    "--sift b --polyphen b --symbol --numbers --biotype --total_length --offline --fields " \
                    "Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE".format(
                    vep_path=self.vep_path, vcf_in=vcf_in, vcf_vep=vcf_out)
        run_stage(self.state, 'apply_vep', vep_command)
        #--cache -i VTVCF --cache --sift b --polyphen b --symbol --numbers --biotype --total_length
        #--force_overwrite --fork THREADS -o VEPVCF --vcf --offline --fields Consequence,Codons,Amino_acids,Gene,
        #SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE

    def apply_bcf(self, inputs, vcf_out):
        '''Apply BCF'''
        vcf_in = inputs
        cores = self.get_stage_options('apply_bcf', 'cores')
        command = "bcftools filter -e \"ALT='*'\" {vcf_in} > {vcf_out}".format(cores=cores,
                            vcf_in=vcf_in, vcf_out=vcf_out)
        run_stage(self.state, 'apply_bcf', command)

    def apply_snpeff(self, inputs, vcf_out):
        '''Apply SnpEFF'''
        vcf_in = inputs
        #cores = self.get_stage_options('apply_snpeff', 'cores')
        snpeff_command = "eff -c {snpeff_conf} -canon GRCh37.75 {vcf_in} > {vcf_out}".format(
                    snpeff_conf=self.snpeff_conf, vcf_in=vcf_in, vcf_out=vcf_out)
        self.run_snpeff('apply_snpeff', snpeff_command)
        #run_snpeff(self.state, 'apply_snpeff', snpeff_command)

    def apply_vcfanno(self, inputs, vcf_out):
        '''Apply anno'''
        vcf_in = inputs
        #cores = self.get_stage_options('apply_snpeff', 'cores')
        anno_command = "./vcfanno_linux64 -lau {annolau} {anno} {vcf_in} > {vcf_out}".format(
                    anno=self.anno, vcf_in=vcf_in, vcf_out=vcf_out)
        self.run_snpeff('apply_vcfanno', anno_command)
        #run_snpeff(self.state, 'apply_snpeff', snpeff_command)
