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

    def align_bwa(self, inputs, bam_out, read_id, lib, lane, sample_id):
        # def align_bwa(self, inputs, bam_out, sample_id):
        '''Align the paired end fastq files to the reference genome using bwa'''
        fastq_read1_in, fastq_read2_in = inputs
        cores = self.get_stage_options('align_bwa', 'cores')
        safe_make_dir('alignments/{sample}'.format(sample=sample_id))
        read_group = '"@RG\\tID:{readid}\\tSM:{sample}\\tPU:lib1\\tLN:{lane}\\tPL:Illumina"' \
            .format(readid=read_id, lib=lib, lane=lane, sample=sample_id)
        command = 'bwa mem -t {cores} -R {read_group} {reference} {fastq_read1} {fastq_read2} ' \
                  '| samtools view -b -h -o {bam} -' \
                  .format(cores=cores,
                          read_group=read_group,
                          fastq_read1=fastq_read1_in,
                          fastq_read2=fastq_read2_in,
                          reference=self.reference,
                          bam=bam_out)
        run_stage(self.state, 'align_bwa', command)

    def sort_bam_picard(self, bam_in, sorted_bam_out):
        '''Sort the BAM file using Picard'''
        picard_args = 'SortSam INPUT={bam_in} OUTPUT={sorted_bam_out} ' \
                      'VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate ' \
                      'MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=True'.format(
                          bam_in=bam_in, sorted_bam_out=sorted_bam_out)
        self.run_picard('sort_bam_picard', picard_args)

    # def mark_duplicates_picard(self, bam_in, outputs):
    #     '''Mark duplicate reads using Picard'''
    #     dedup_bam_out, metrics_out = outputs
    #     picard_args = 'MarkDuplicates INPUT={bam_in} OUTPUT={dedup_bam_out} ' \
    #                   'METRICS_FILE={metrics_out} VALIDATION_STRINGENCY=LENIENT ' \
    #                   'MAX_RECORDS_IN_RAM=5000000 ASSUME_SORTED=True ' \
    #                   'CREATE_INDEX=True'.format(bam_in=bam_in, dedup_bam_out=dedup_bam_out,
    #                                              metrics_out=metrics_out)
    #     self.run_picard('mark_duplicates_picard', picard_args)

    # def realigner_target_creator(self, inputs, intervals_out):
    #     '''Generate chromosome intervals using GATK'''
    #     bam_in, _metrics_dup = inputs
    #     cores = self.get_stage_options('chrom_intervals_gatk', 'cores')
    #     gatk_args = '-T RealignerTargetCreator -R {reference} -I {bam} ' \
    #                 '--num_threads {threads} --known {mills_hg19} ' \
    #                 '--known {one_k_g_indels} ' \
    #                 '--known {one_k_g_indels} ' \
    #                 '-o {out}'.format(reference=self.reference, bam=bam_in,
    #                                   threads=cores, mills_hg19=self.mills_hg19,
    #                                   one_k_g_indels=self.one_k_g_indels,
    #                                   out=intervals_out)
    #     self.run_gatk('chrom_intervals_gatk', gatk_args)
    #
    # def local_realignment_gatk(self, inputs, bam_out):
    #     '''Local realign reads using GATK'''
    #     target_intervals_in, bam_in = inputs
    #     gatk_args = "-T IndelRealigner -R {reference} -I {bam} " \
    #                 "-targetIntervals {target_intervals} -known {mills_hg19} " \
    #                 "-known {one_k_g_indels} " \
    #                 "-o {out}".format(reference=self.reference, bam=bam_in,
    #                                   mills_hg19=self.mills_hg19,
    #                                   one_k_g_indels=self.one_k_g_indels,
    #                                   target_intervals=target_intervals_in,
    #                                   out=bam_out)
    #     self.run_gatk('local_realignment_gatk', gatk_args)

    # XXX I'm not sure that --num_cpu_threads_per_data_thread has any benefit
    # here
    # def base_recalibration_gatk(self, bam_in, outputs):
    #     '''Base recalibration using GATK'''
    #     csv_out, log_out = outputs
    #     gatk_args = "-T BaseRecalibrator -R {reference} -I {bam} " \
    #                 "--num_cpu_threads_per_data_thread 4 --knownSites {dbsnp_hg19} " \
    #                 "--knownSites {mills_hg19} --knownSites {one_k_g_indels} " \
    #                 "-log {log} -o {out}".format(reference=self.reference, bam=bam_in,
    #                                              mills_hg19=self.mills_hg19, dbsnp_hg19=self.dbsnp_hg19,
    #                                              one_k_g_indels=self.one_k_g_indels,
    #                                              log=log_out, out=csv_out)
    #     self.run_gatk('base_recalibration_gatk', gatk_args)

    # XXX I'm not sure that --num_cpu_threads_per_data_thread has any benefit
    # here
    # def print_reads_gatk(self, inputs, bam_out):
    #     '''Print reads using GATK'''
    #     [csv_in, _log], bam_in = inputs
    #     gatk_args = "-T PrintReads -R {reference} -I {bam} --BQSR {recal_csv} " \
    #                 "-o {out} --num_cpu_threads_per_data_thread 4".format(reference=self.reference,
    #                                              bam=bam_in, recal_csv=csv_in, out=bam_out)
    #     self.run_gatk('print_reads_gatk', gatk_args)

    # Merge per lane bam into a single bam per sample
    # def merge_sample_bams(self, bam_files_in, bam_out):
    #     '''Merge per lane bam into a merged bam file'''
    #     bam_files = ' '.join(['INPUT=' + bam for bam in bam_files_in])
    #     picard_args = 'MergeSamFiles {bams_in} OUTPUT={merged_bam_out} ' \
    #                   'VALIDATION_STRINGENCY=LENIENT ' \
    #                   'MAX_RECORDS_IN_RAM=5000000 ASSUME_SORTED=True ' \
    #                   'CREATE_INDEX=True'.format(
    #                       bams_in=bam_files, merged_bam_out=bam_out)
    #     self.run_picard('merge_sample_bams', picard_args)

    # def call_haplotypecaller_gatk(self, bam_in, vcf_out):
    #     '''Call variants using GATK'''
    #     # safe_make_dir('variants}'.format(sample=sample_id))
    #     gatk_args = "-T HaplotypeCaller -R {reference} --min_base_quality_score 20 " \
    #                 "--emitRefConfidence GVCF " \
    #                 "-A AlleleBalance -A AlleleBalanceBySample " \
    #                 "-A ChromosomeCounts -A ClippingRankSumTest " \
    #                 "-A Coverage -A DepthPerAlleleBySample " \
    #                 "-A DepthPerSampleHC -A FisherStrand " \
    #                 "-A GCContent -A GenotypeSummaries " \
    #                 "-A HardyWeinberg -A HomopolymerRun " \
    #                 "-A LikelihoodRankSumTest -A LowMQ " \
    #                 "-A MappingQualityRankSumTest -A MappingQualityZero " \
    #                 "-A QualByDepth " \
    #                 "-A RMSMappingQuality -A ReadPosRankSumTest " \
    #                 "-A SampleList -A SpanningDeletions " \
    #                 "-A StrandBiasBySample -A StrandOddsRatio " \
    #                 "-A TandemRepeatAnnotator -A VariantType " \
    #                 "-I {bam} -L {interval_list} -o {out}".format(reference=self.reference,
    #                                                               bam=bam_in, interval_list=self.interval_hg19, out=vcf_out)
    #     self.run_gatk('call_haplotypecaller_gatk', gatk_args)
    #
    # def call_haplotypecaller_gatk_nct(self, bam_in, vcf_out):
    #     '''Call variants using GATK'''
    #     # safe_make_dir('variants}'.format(sample=sample_id))
    #     gatk_args = "-T HaplotypeCaller -R {reference} --min_base_quality_score 20 " \
    #                 "--standard_min_confidence_threshold_for_calling 30.0 " \
    #                 "--num_cpu_threads_per_data_thread 4 " \
    #                 "--variant_index_type LINEAR " \
    #                 "--standard_min_confidence_threshold_for_emitting 30.0 " \
    #                 "-I {bam} -L {interval_list} -o {out}".format(reference=self.reference,
    #                                                               bam=bam_in, interval_list=self.interval_hg19, out=vcf_out)
    #     self.run_gatk('call_haplotypecaller_gatk', gatk_args)
    #
    # def combine_gvcf_gatk(self, vcf_files_in, vcf_out):
    #     '''Combine G.VCF files for all samples using GATK'''
    #     g_vcf_files = ' '.join(['--variant ' + vcf for vcf in vcf_files_in])
    #     gatk_args = "-T CombineGVCFs -R {reference} " \
    #                 "--disable_auto_index_creation_and_locking_when_reading_rods " \
    #                 "{g_vcf_files} -o {vcf_out}".format(reference=self.reference,
    #                                                     g_vcf_files=g_vcf_files, vcf_out=vcf_out)
    #     # "{g_vcf_files} -o {vcf_out} --variant {CEU}".format(reference=self.reference,
    #     # g_vcf_files=g_vcf_files, vcf_out=vcf_out, CEU=self.CEU_mergeGvcf)
    #     self.run_gatk('combine_gvcf_gatk', gatk_args)
    #
    # def genotype_gvcf_gatk(self, merged_vcf_in, vcf_out):
    #     '''Genotype G.VCF files using GATK'''
    #     cores = self.get_stage_options('genotype_gvcf_gatk', 'cores')
    #     gatk_args = "-T GenotypeGVCFs -R {reference} " \
    #                 "--disable_auto_index_creation_and_locking_when_reading_rods " \
    #                 "-A AlleleBalance -A AlleleBalanceBySample " \
    #                 "-A ChromosomeCounts -A ClippingRankSumTest " \
    #                 "-A Coverage -A DepthPerAlleleBySample " \
    #                 "-A DepthPerSampleHC -A FisherStrand " \
    #                 "-A GCContent -A GenotypeSummaries " \
    #                 "-A HardyWeinberg -A HomopolymerRun " \
    #                 "-A LikelihoodRankSumTest " \
    #                 "-A MappingQualityRankSumTest -A MappingQualityZero " \
    #                 "-A QualByDepth " \
    #                 "-A RMSMappingQuality -A ReadPosRankSumTest " \
    #                 "-A SampleList -A SpanningDeletions " \
    #                 "-A StrandBiasBySample -A StrandOddsRatio " \
    #                 "-A TandemRepeatAnnotator -A VariantType " \
    #                 "--dbsnp {dbsnp} " \
    #                 "--num_threads {cores} --variant {merged_vcf} --out {vcf_out}" \
    #                 .format(reference=self.reference, dbsnp=self.dbsnp_hg19,
    #                         cores=cores, merged_vcf=merged_vcf_in, vcf_out=vcf_out)
    #     self.run_gatk('genotype_gvcf_gatk', gatk_args)

    # def genotype_gvcf_gatk(self, merged_vcf_in, vcf_out):
    #     '''Genotype G.VCF files using GATK'''
    #     cores = self.get_stage_options('genotype_gvcf_gatk', 'cores')
    #     gatk_args = "-T GenotypeGVCFs -R {reference} " \
    #                 "--disable_auto_index_creation_and_locking_when_reading_rods " \
    #                 "--num_threads {cores} --variant {merged_vcf} --out {vcf_out} " \
    #                 "--variant {CEU_mergeGvcf} --variant {GBR_mergeGvcf} " \
    #                 "--variant {FIN_mergeGvcf}".format(reference=self.reference,
    #                         cores=cores, merged_vcf=merged_vcf_in, vcf_out=vcf_out,
    #                         CEU_mergeGvcf=self.CEU_mergeGvcf, GBR_mergeGvcf=self.GBR_mergeGvcf,
    #                         FIN_mergeGvcf=self.FIN_mergeGvcf)
    #     self.run_gatk('genotype_gvcf_gatk', gatk_args)

    # def snp_recalibrate_gatk(self, genotype_vcf_in, outputs):
    #     '''SNP recalibration using GATK'''
    #     recal_snp_out, tranches_snp_out, snp_plots_r_out = outputs
    #     cores = self.get_stage_options('snp_recalibrate_gatk', 'cores')
    #     gatk_args = "-T VariantRecalibrator --disable_auto_index_creation_and_locking_when_reading_rods " \
    #                 "-R {reference} --minNumBadVariants 5000 --num_threads {cores} " \
    #                 "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} " \
    #                 "-resource:omni,known=false,training=true,truth=true,prior=12.0 {one_k_g_snps} " \
    #                 "-resource:1000G,known=false,training=true,truth=false,prior=10.0 {one_k_g_highconf_snps} " \
    #                 "-an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR " \
    #                 "-input {genotype_vcf} --recal_file {recal_snp} --tranches_file {tranches_snp} " \
    #                 "-rscriptFile {snp_plots} -mode SNP".format(reference=self.reference,
    #                                                             cores=cores, hapmap=self.hapmap, one_k_g_snps=self.one_k_g_snps,
    #                                                             one_k_g_highconf_snps=self.one_k_g_highconf_snps, genotype_vcf=genotype_vcf_in,
    #                                                             recal_snp=recal_snp_out, tranches_snp=tranches_snp_out, snp_plots=snp_plots_r_out)
    #     self.run_gatk('snp_recalibrate_gatk', gatk_args)
    #
    # def indel_recalibrate_gatk(self, genotype_vcf_in, outputs):
    #     '''INDEL recalibration using GATK'''
    #     recal_indel_out, tranches_indel_out, indel_plots_r_out = outputs
    #     cores = self.get_stage_options('indel_recalibrate_gatk', 'cores')
    #     gatk_args = "-T VariantRecalibrator --disable_auto_index_creation_and_locking_when_reading_rods " \
    #                 "-R {reference} --minNumBadVariants 5000 --num_threads {cores} " \
    #                 "-resource:mills,known=false,training=true,truth=true,prior=12.0 {mills_hg19} " \
    #                 "-resource:1000G,known=false,training=true,truth=true,prior=10.0 {one_k_g_indels} " \
    #                 "-an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR " \
    #                 "-input {genotype_vcf} -recalFile {recal_indel} " \
    #                 "-tranchesFile {tranches_indel} -rscriptFile {indel_plots} " \
    #                 " -mode INDEL --maxGaussians 4".format(reference=self.reference,
    #                                       cores=cores, mills_hg19=self.mills_hg19, one_k_g_indels=self.one_k_g_indels,
    #                                       genotype_vcf=genotype_vcf_in, recal_indel=recal_indel_out,
    #                                       tranches_indel=tranches_indel_out, indel_plots=indel_plots_r_out)
    #     self.run_gatk('indel_recalibrate_gatk', gatk_args)
    #
    # def apply_snp_recalibrate_gatk(self, inputs, vcf_out):
    #     '''Apply SNP recalibration using GATK'''
    #     genotype_vcf_in, [recal_snp, tranches_snp] = inputs
    #     cores = self.get_stage_options('apply_snp_recalibrate_gatk', 'cores')
    #     gatk_args = "-T ApplyRecalibration --disable_auto_index_creation_and_locking_when_reading_rods " \
    #                 "-R {reference} --ts_filter_level 99.5 --excludeFiltered --num_threads {cores} " \
    #                 "-input {genotype_vcf} -recalFile {recal_snp} -tranchesFile {tranches_snp} " \
    #                 "-mode SNP -o {vcf_out}".format(reference=self.reference,
    #                                                 cores=cores, genotype_vcf=genotype_vcf_in, recal_snp=recal_snp,
    #                                                 tranches_snp=tranches_snp, vcf_out=vcf_out)
    #     self.run_gatk('apply_snp_recalibrate_gatk', gatk_args)
    #
    # def apply_indel_recalibrate_gatk(self, inputs, vcf_out):
    #     '''Apply INDEL recalibration using GATK'''
    #     genotype_vcf_in, [recal_indel, tranches_indel] = inputs
    #     cores = self.get_stage_options('apply_indel_recalibrate_gatk', 'cores')
    #     gatk_args = "-T ApplyRecalibration --disable_auto_index_creation_and_locking_when_reading_rods " \
    #                 "-R {reference} --ts_filter_level 99.0 --excludeFiltered --num_threads {cores} " \
    #                 "-input {genotype_vcf} -recalFile {recal_indel} -tranchesFile {tranches_indel} " \
    #                 "-mode INDEL -o {vcf_out}".format(reference=self.reference,
    #                                                   cores=cores, genotype_vcf=genotype_vcf_in, recal_indel=recal_indel,
    #                                                   tranches_indel=tranches_indel, vcf_out=vcf_out)
    #     self.run_gatk('apply_indel_recalibrate_gatk', gatk_args)

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
                    "--fork {threads} --everything --offline --coding_only --no_intergenic".format(
                    vep_path=self.vep_path, vcf_in=vcf_in, vcf_vep=vcf_out, threads=cores)
        run_stage(self.state, 'apply_vep', vep_command)

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

    def combine_variants_gatk(self, inputs, vcf_out):
        '''Combine variants using GATK'''
        recal_snp, [recal_indel] = inputs
        cores = self.get_stage_options('combine_variants_gatk', 'cores')
        gatk_args = "-T CombineVariants -R {reference} --disable_auto_index_creation_and_locking_when_reading_rods " \
                    "--num_threads {cores} --genotypemergeoption UNSORTED --variant {recal_snp} " \
                    "--variant {recal_indel} -o {vcf_out}".format(reference=self.reference,
                                                                  cores=cores, recal_snp=recal_snp, recal_indel=recal_indel,
                                                                  vcf_out=vcf_out)
        self.run_gatk('combine_variants_gatk', gatk_args)

    def select_variants_gatk(self, combined_vcf, vcf_out):
        '''Select variants using GATK'''
        gatk_args = "-T SelectVariants -R {reference} --disable_auto_index_creation_and_locking_when_reading_rods " \
                    "--variant {combined_vcf} -select 'DP > 100' -o {vcf_out}".format(reference=self.reference,
                                                                                      combined_vcf=combined_vcf, vcf_out=vcf_out)
        self.run_gatk('select_variants_gatk', gatk_args)
