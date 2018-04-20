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
        # Hi-Plex example: OHI031002-P02F04_S318_L001_R1_001.fastq
        # new sample name = OHI031002-P02F04
        filter=formatter(
            '.+/(?P<sample>[a-zA-Z0-9-]+)_(?P<readid>[a-zA-Z0-9-]+)_(?P<lane>[a-zA-Z0-9]+)_R1_(?P<lib>[a-zA-Z0-9-:]+).fastq'),

        # Add one more inputs to the stage:
        #    1. The corresponding R2 FASTQ file
        # Hi-Plex example: OHI031002-P02F04_S318_L001_R2_001.fastq
        add_inputs=add_inputs(
            '{path[0]}/{sample[0]}_{readid[0]}_{lane[0]}_R2_{lib[0]}.fastq'),

        # Add an "extra" argument to the state (beyond the inputs and outputs)
        # which is the sample name. This is needed within the stage for finding out
        # sample specific configuration options
        extras=['{sample[0]}', '{readid[0]}', '{lane[0]}', '{lib[0]}'],
        # The output file name is the sample name with a .bam extension.
        output='alignments/{sample[0]}_{readid[0]}/{sample[0]}_{readid[0]}.bam')

    # Call variants using undr_rover
    pipeline.transform(
        task_func=stages.apply_undr_rover,
        name='apply_undr_rover',
        input=output_from('original_fastqs'),
        # Match the R1 (read 1) FASTQ file and grab the path and sample name.
        # This will be the first input to the stage.
        filter=formatter(
            '.+/(?P<sample>[a-zA-Z0-9-]+)_(?P<readid>[a-zA-Z0-9-]+)_(?P<lane>[a-zA-Z0-9]+)_R1_(?P<lib>[a-zA-Z0-9-:]+).fastq'),

        add_inputs=add_inputs(
            '{path[0]}/{sample[0]}_{readid[0]}_{lane[0]}_R2_{lib[0]}.fastq'),
        # extras=['{sample[0]}', '{readid[0]}', '{lane[0]}', '{lib[0]}'],
        extras=['{sample[0]}', '{readid[0]}'],

        # The output file name is the sample name with a .bam extension.
        output='variants/undr_rover/{sample[0]}_{readid[0]}.vcf')

    # Sort the BAM file using Picard
    pipeline.transform(
        task_func=stages.sort_bam_picard,
        name='sort_bam_picard',
        input=output_from('align_bwa'),
        filter=suffix('.bam'),
        output='.sort.bam')

    # High quality and primary alignments
    pipeline.transform(
        task_func=stages.primary_bam,
        name='primary_bam',
        input=output_from('sort_bam_picard'),
        filter=suffix('.sort.bam'),
        output='.primary.bam')

    # index bam file
    pipeline.transform(
        task_func=stages.index_sort_bam_picard,
        name='index_bam',
        input=output_from('primary_bam'),
        filter=suffix('.primary.bam'),
        output='.primary.bam.bai')

    # Clip the primer_seq from BAM File
    (pipeline.transform(
        task_func=stages.clip_bam,
        name='clip_bam',
        input=output_from('primary_bam'),
        filter=suffix('.primary.bam'),
        output='.primary.primerclipped.bam')
        .follows('index_bam'))

    ###### GATK VARIANT CALLING ######
    # Call variants using GATK
    pipeline.transform(
        task_func=stages.call_haplotypecaller_gatk,
        name='call_haplotypecaller_gatk',
        input=output_from('clip_bam'),
        # filter=suffix('.merged.dedup.realn.bam'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9-_]+).primary.primerclipped.bam'),
        output='variants/gatk/{sample[0]}.g.vcf')
        # .follows('index_sort_bam_picard'))

    # Combine G.VCF files for all samples using GATK
    pipeline.merge(
        task_func=stages.combine_gvcf_gatk,
        name='combine_gvcf_gatk',
        input=output_from('call_haplotypecaller_gatk'),
        output='variants/gatk/ALL.combined.vcf')

    # Genotype G.VCF files using GATK
    pipeline.transform(
        task_func=stages.genotype_gvcf_gatk,
        name='genotype_gvcf_gatk',
        input=output_from('combine_gvcf_gatk'),
        filter=suffix('.combined.vcf'),
        output='.raw.vcf')

    # Annotate VCF file using GATK
    pipeline.transform(
       task_func=stages.variant_annotator_gatk,
       name='variant_annotator_gatk',
       input=output_from('genotype_gvcf_gatk'),
       filter=suffix('.raw.vcf'),
       output='.raw.annotate.vcf')

    # Apply VariantFiltration using GATK
    pipeline.transform(
        task_func=stages.apply_variant_filtration_gatk_lenient,
        name='apply_variant_filtration_gatk_lenient',
        input=output_from('variant_annotator_gatk'),
        filter=suffix('.raw.annotate.vcf'),
        output='.raw.annotate.filtered_lenient.vcf')

    return pipeline
