'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from
from stages import Stages


def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name='hiplexpipe')
    # Get a list of paths to all the FASTQ files
    fastq_files = state.config.get_option('fastqs')
    # Stages are dependent on the state
    stages = Stages(state)

    # The original FASTQ files
    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.original_fastqs,
        name='original_fastqs',
        output=fastq_files)

    # Align paired end reads in FASTQ to the reference producing a BAM file
    pipeline.transform(
        task_func=stages.align_bwa,
        name='align_bwa',
        input=output_from('original_fastqs'),
        # Match the R1 (read 1) FASTQ file and grab the path and sample name.
        # This will be the first input to the stage.
        # We assume the sample name may consist of only alphanumeric
        # characters.
        # filter=formatter('(?P<path>.+)/
        # (?P<readid>[a-zA-Z0-9-\.]+)_(?P<lib>[a-zA-Z0-9-]+)_(?P<lane>[a-zA-Z0-9]+)_(?P<sample>[a-zA-Z0-9]+)_1.fastq.gz'),
        # 00-002-0640_S192_L001_R1_001.fastq

        filter=formatter(
            '.+/(?P<readid>[a-zA-Z0-9-]+)_(?P<sample>[a-zA-Z0-9-]+)_(?P<lane>[a-zA-Z0-9]+)_R1_(?P<lib>[a-zA-Z0-9-:]+).fastq.gz'),
        # Add one more inputs to the stage:
        #    1. The corresponding R2 FASTQ file
        # e.g. C2WPF.5_Solexa-201237_5_X4311_1.fastq.gz
        add_inputs=add_inputs(
            '{path[0]}/{readid[0]}_{sample[0]}_{lane[0]}_R2_{lib[0]}.fastq.gz'),
        # Add an "extra" argument to the state (beyond the inputs and outputs)
        # which is the sample name. This is needed within the stage for finding out
        # sample specific configuration options
        extras=['{readid[0]}', '{lib[0]}', '{lane[0]}', '{sample[0]}'],
        # extras=['{sample[0]}'],
        # The output file name is the sample name with a .bam extension.
        output='alignments/{sample[0]}/{sample[0]}.bam')

    # Sort the BAM file using Picard
    pipeline.transform(
        task_func=stages.sort_bam_picard,
        name='sort_bam_picard',
        input=output_from('align_bwa'),
        filter=suffix('.bam'),
        output='.sort.bam')

    # Apply samtools
    (pipeline.merge(
        task_func=stages.apply_samtools_mpileup,
        name='apply_samtools_mpileup',
        input=output_from('sort_bam_picard'),
        filter=suffix('.sort.bam'),
        # add_inputs=add_inputs(['variants/ALL.indel_recal', 'variants/ALL.indel_tranches']),
        output='.bcf')
        .follows('sort_bam_picard'))

    # Apply bcftools
    (pipeline.merge(
        task_func=stages.apply_samtools_mpileup,
        name='apply_bcftools',
        input=output_from('apply_samtools_mpileup'),
        filter=suffix('.bcf'),
        # add_inputs=add_inputs(['variants/ALL.indel_recal', 'variants/ALL.indel_tranches']),
        output='.raw.vcf')
        .follows('apply_samtools_mpileup'))

    # Apply NORM
    (pipeline.transform(
        task_func=stages.apply_vt,
        name='apply_vt',
        input=output_from('apply_bcftools'),
        filter=suffix('.raw.vcf'),
        # add_inputs=add_inputs(['variants/ALL.indel_recal', 'variants/ALL.indel_tranches']),
        output='.raw.vt.vcf')
        .follows('sort_bam_picard'))

    # Apply VEP
    (pipeline.transform(
        task_func=stages.apply_vep,
        name='apply_vep',
        input=output_from('apply_vt'),
        filter=suffix('.sort.bam'),
        # add_inputs=add_inputs(['variants/ALL.indel_recal', 'variants/ALL.indel_tranches']),
        output='.raw.vt.vep.vcf')
        .follows('apply_vt'))

    # Apply BCF
    (pipeline.transform(
        task_func=stages.apply_bcf,
        name='apply_bcf',
        input=output_from('apply_vep'),
        filter=suffix('.raw.vqsr.vt.vep.vcf'),
        # add_inputs=add_inputs(['variants/ALL.indel_recal', 'variants/ALL.indel_tranches']),
        output='.raw.vt.vep.bcf.vcf')
        .follows('apply_vep'))

    # Apply SnpEff
    (pipeline.transform(
        task_func=stages.apply_snpeff,
        name='apply_snpeff',
        input=output_from('apply_bcf'),
        filter=suffix('.raw.vqsr.vt.vep.bcf.vcf'),
        # add_inputs=add_inputs(['variants/ALL.indel_recal', 'variants/ALL.indel_tranches']),
        output='.annotated.vcf')
        .follows('apply_bcf'))

    # Combine variants using GATK
    # (pipeline.transform(
    #     task_func=stages.combine_variants_gatk,
    #     name='combine_variants_gatk',
    #     input=output_from('apply_snp_recalibrate_gatk'),
    #     filter=suffix('.recal_SNP.vcf'),
    #     add_inputs=add_inputs(['ALL.recal_INDEL.vcf']),
    #     # output='.combined.vcf')
    #     output='ALL.raw.vqsr.vcf')
    #     .follows('apply_indel_recalibrate_gatk'))
    #
    # # Select variants using GATK
    # pipeline.transform(
    #     task_func=stages.select_variants_gatk,
    #     name='select_variants_gatk',
    #     input=output_from('combine_variants_gatk'),
    #     filter=suffix('.combined.vcf'),
    #     output='.selected.vcf')


    return pipeline
