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
# PICARD_JAR = '/vlsci/VR0002/kmahmood/Programs/Picard/picard-tools-2.8.3/picard.jar'
PICARD_JAR = '/usr/local/easybuild/software/picard/2.3.0/picard.jar'
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
        self.dbsnp_hg19 = self.get_options('dbsnp_hg19')
        self.mills_hg19 = self.get_options('mills_hg19')
        self.one_k_g_snps = self.get_options('one_k_g_snps')
        self.one_k_g_indels = self.get_options('one_k_g_indels')
        self.one_k_g_highconf_snps = self.get_options('one_k_g_highconf_snps')
        self.hapmap = self.get_options('hapmap')
        # self.interval_hg19 = self.get_options('exome_bed_hg19')
        # self.CEU_mergeGvcf = self.get_options('CEU_mergeGvcf')
        self.snpeff_conf = self.get_options('snpeff_conf')
        self.bamclipper = self.get_options('bamclipper')
        self.vep_path = self.get_options('vep_path')
        self.vt_path = self.get_options('vt_path')
        self.coord_file = self.get_options('coord_file')
        self.target_bed = self.get_options('target_bed')
        # self.interval_file = self.get_options('interval_file')
        self.primer_file = self.get_options('primer_file')
        self.primer_bedpe_file = self.get_options('primer_bedpe_file')
        self.proportionthresh = self.get_options('proportionthresh')
        self.absthresh = self.get_options('absthresh')
        self.maxvariants = self.get_options('maxvariants')
        # self.fragment_bed = self.get_options('fragment_bed')
        self.annolua = self.get_options('annolua')
        self.anno = self.get_options('anno')
        self.hrfile = self.get_options('hrfile')
        self.other_vep = self.get_options('other_vep')
        self.snpeff_path = self.get_options('snpeff_path')
        self.gatk_bed = self.get_options('gatk_bed')

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
        safe_make_dir('alignments/{sample}_{readid}'.format(sample=sample_id, readid=read_id))
        read_group = '"@RG\\tID:{readid}\\tSM:{sample}_{readid}\\tPU:lib1\\tLN:{lane}\\tPL:Illumina"' \
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

    def apply_undr_rover(self, inputs, vcf_output, sample_id, readid):
        # def align_bwa(self, inputs, bam_out, sample_id):
        '''Apply undr_rover to call variants from paired end fastq files'''
        fastq_read1_in, fastq_read2_in = inputs
        cores = self.get_stage_options('apply_undr_rover', 'cores')
        safe_make_dir('variants/undr_rover')
        safe_make_dir('variants/undr_rover/coverdir')
        coverfile = "variants/undr_rover/coverdir/" + sample_id + "_" + readid + ".coverage"
        # read_group = '"@RG\\tID:{readid}\\tSM:{sample}_{readid}\\tPU:lib1\\tLN:{lane}\\tPL:Illumina"' \
            # .format(readid=read_id, lib=lib, lane=lane, sample=sample_id)

        command = 'undr_rover --primer_coords {coord_file} ' \
                  '--primer_sequences {primer_file} ' \
                  '--reference {reference} ' \
                  '--out {vcf_output} ' \
                  '--coverfile {coverfile} ' \
                  '--proportionthresh {proportionthresh} ' \
                  '--absthresh {absthresh} ' \
                  '--max_variants {maxvariants} ' \
                  '{fastq_read1} {fastq_read2}'.format(
                        coord_file=self.coord_file, primer_file=self.primer_file,
                        reference=self.reference,
                        vcf_output=vcf_output,
                        #coverdir=self.coverdir,
                        proportionthresh=self.proportionthresh,
                        absthresh=self.absthresh,
                        maxvariants=self.maxvariants,
                        coverfile=coverfile,
                        fastq_read1=fastq_read1_in,
                        fastq_read2=fastq_read2_in)
        run_stage(self.state, 'apply_undr_rover', command)

    def clip_bam(self, bam_in, sorted_bam_out):
        '''Clip the BAM file using Bamclipper'''
        bamclipper_args = '{bamclipper} -b {bam_in} -p {primer_bedpe_file} -n 1'.format(
                          bamclipper=self.bamclipper, bam_in=bam_in, primer_bedpe_file=self.primer_bedpe_file)
        run_stage(self.state, 'clip_bam', bamclipper_args)

    def sort_bam_picard(self, bam_in, sorted_bam_out):
        '''Sort the BAM file using Picard'''
        picard_args = 'SortSam INPUT={bam_in} OUTPUT={sorted_bam_out} ' \
                      'VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate ' \
                      'MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=True'.format(
                          bam_in=bam_in, sorted_bam_out=sorted_bam_out)
        self.run_picard('sort_bam_picard', picard_args)

    def primary_bam(self, bam_in, sbam_out):
        '''On keep primary alignments in the BAM file using samtools'''
        command = 'samtools view -h -q 1 -f 2 -F 4 -F 8 -F 256 -b ' \
                    '-o {sbam_out} {bam_in}'.format(
                        bam_in=bam_in, sbam_out=sbam_out)
        run_stage(self.state, 'primary_bam', command)

    # index sorted bam file
    def index_sort_bam_picard(self, bam_in, bam_index):
        '''Index sorted bam using samtools'''
        command = 'samtools index {bam_in} {bam_index}'.format(
                          bam_in=bam_in, bam_index=bam_index)
        run_stage(self.state, 'index_sort_bam_picard', command)

    ##########
    def call_haplotypecaller_gatk(self, bam_in, vcf_out):
        '''Call variants using GATK'''
        safe_make_dir('variants/gatk')
        # safe_make_dir('variants}'.format(sample=sample_id))
        gatk_args = "-T HaplotypeCaller -R {reference} --min_base_quality_score 20 " \
                    "--emitRefConfidence GVCF " \
                    "-A AlleleBalance -A AlleleBalanceBySample " \
                    "-A ChromosomeCounts -A ClippingRankSumTest " \
                    "-A Coverage -A DepthPerAlleleBySample " \
                    "-A DepthPerSampleHC -A FisherStrand " \
                    "-A GCContent -A GenotypeSummaries " \
                    "-A HardyWeinberg -A HomopolymerRun " \
                    "-A LikelihoodRankSumTest -A LowMQ " \
                    "-A MappingQualityRankSumTest -A MappingQualityZero " \
                    "-A QualByDepth " \
                    "-A RMSMappingQuality -A ReadPosRankSumTest " \
                    "-A SampleList -A SpanningDeletions " \
                    "-A StrandBiasBySample -A StrandOddsRatio " \
                    "-A TandemRepeatAnnotator -A VariantType " \
                    "--dontUseSoftClippedBases " \
                    "-I {bam} -o {out}".format(reference=self.reference,
                                                                  bam=bam_in, out=vcf_out)
        self.run_gatk('call_haplotypecaller_gatk', gatk_args)

    def combine_gvcf_gatk(self, vcf_files_in, vcf_out):
        '''Combine G.VCF files for all samples using GATK'''
        g_vcf_files = ' '.join(['--variant ' + vcf for vcf in vcf_files_in])
        gatk_args = "-T CombineGVCFs -R {reference} -L {gatk_bed} " \
                    "--disable_auto_index_creation_and_locking_when_reading_rods " \
                    "{g_vcf_files} -o {vcf_out}".format(reference=self.reference, gatk_bed=self.gatk_bed,
                                                        g_vcf_files=g_vcf_files, vcf_out=vcf_out)
        self.run_gatk('combine_gvcf_gatk', gatk_args)

    def genotype_gvcf_gatk(self, combined_vcf_in, vcf_out):
        '''Genotype G.VCF files using GATK'''
        cores = self.get_stage_options('genotype_gvcf_gatk', 'cores')
        gatk_args = "-T GenotypeGVCFs -R {reference} " \
                    "--disable_auto_index_creation_and_locking_when_reading_rods " \
                    "--dbsnp {dbsnp} -L {gatk_bed} " \
                    "--num_threads {cores} --variant {combined_vcf} --out {vcf_out}" \
                    .format(reference=self.reference, dbsnp=self.dbsnp_hg19, gatk_bed=self.gatk_bed,
                            cores=cores, combined_vcf=combined_vcf_in, vcf_out=vcf_out)
        self.run_gatk('genotype_gvcf_gatk', gatk_args)
