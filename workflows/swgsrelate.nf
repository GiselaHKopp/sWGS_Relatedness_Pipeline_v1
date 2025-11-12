/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_swgsrelate_pipeline'

include { BASE_QUALITY_SCORE_RECALIBRATION } from '../subworkflows/local/base_quality_score_recalibration'
include { BOOTSTRAP_VARIANT_SET as BOOTSTRAP_VARIANT_SET_1 } from '../subworkflows/local/bootstrap_variant_set'
include { BOOTSTRAP_VARIANT_SET as BOOTSTRAP_VARIANT_SET_2 } from '../subworkflows/local/bootstrap_variant_set'
include { PREPARE_GENOME                   } from '../subworkflows/local/prepare_genome'
include { PREPARE_INTERVALS                } from '../subworkflows/local/prepare_intervals'
include { PREPROCESS                       } from '../subworkflows/local/preprocess'
include { CALL_VARIANTS_BCFTOOLS           } from '../subworkflows/local/call_variants_bcftools'
include { CALL_VARIANTS_GATK               } from '../subworkflows/local/call_variants_gatk'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SWGSRELATE {

    take:
    samplesheet // channel: [ meta, list(fastq) ]

    main:
    // Channel for collecting software versions
    ch_versions = channel.empty()

    // Channel for collecting MultiQC files
    ch_multiqc_files = channel.empty()

    // Define reference genome and index
    ch_fasta = params.fasta ?
    channel.fromPath(params.fasta)
        .map { f -> [ [id: f.baseName], f ] }
        .collect()
    : channel.empty()

    //
    // SUBWORKFLOW: PREPARE_GENOME
    //
    PREPARE_GENOME(ch_fasta)
    ch_fai = PREPARE_GENOME.out.fai
    ch_dict = PREPARE_GENOME.out.dict
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW: PREPARE_INTERVALS
    //
    PREPARE_INTERVALS(ch_fai)
    ch_intervals_split = PREPARE_INTERVALS.out.intervals_split
    ch_versions = ch_versions.mix(PREPARE_INTERVALS.out.versions)

    if(params.stages.contains('preprocess')) {
        //
        // SUBWORKFLOW: PREPROCESS
        //
        ch_preprocessed = PREPROCESS(samplesheet, ch_fasta, ch_fai)
        ch_cram = ch_preprocessed.cram
        ch_crai = ch_preprocessed.crai
        ch_versions = ch_versions.mix(ch_preprocessed.versions)
        ch_multiqc_files = ch_multiqc_files.mix(ch_preprocessed.multiqc_files)
    } else {
        // TODO: Load BAMs from samplesheet without preprocessing
        //ch_cram = channel.empty()
        //ch_crai = channel.empty()
        //ch_fai = channel.empty()
        //ch_dict = channel.empty()
    }

    if(params.stages.contains('bootstrap_variant_set')) {
        //
        // SUBWORKFLOW: BOOTSTRAP_VARIANT_SET - ROUND 1
        //
        ch_cram.map { meta, cram_file ->
            tuple( meta + ['bootstrapping_round': 1], cram_file ) }
            .set { ch_cram }
        ch_crai.map { meta, crai_file ->
            tuple( meta + ['bootstrapping_round': 1], crai_file ) }
            .set { ch_crai }
        BOOTSTRAP_VARIANT_SET_1(
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_intervals_split,
            ch_cram,
            ch_crai
        )
        ch_versions = ch_versions.mix(BOOTSTRAP_VARIANT_SET_1.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(BOOTSTRAP_VARIANT_SET_1.out.multiqc_files)
        ch_cram = BOOTSTRAP_VARIANT_SET_1.out.cram
        ch_crai = BOOTSTRAP_VARIANT_SET_1.out.crai
        ch_vcf  = BOOTSTRAP_VARIANT_SET_1.out.vcf
        ch_tbi  = BOOTSTRAP_VARIANT_SET_1.out.tbi

        //
        // SUBWORKFLOW: BOOTSTRAP_VARIANT_SET - ROUND 2
        //
        ch_cram.map { meta, cram_file ->
            tuple( meta + ['bootstrapping_round': 2], cram_file ) }
            .set { ch_cram }
        ch_crai.map { meta, crai_file ->
            tuple( meta + ['bootstrapping_round': 2], crai_file ) }
            .set { ch_crai }
        if (params.bqsr_rounds > 1) {
            BOOTSTRAP_VARIANT_SET_2(
                ch_fasta,
                ch_fai,
                ch_dict,
                ch_intervals_split,
                ch_cram,
                ch_crai
            )
            ch_cram = BOOTSTRAP_VARIANT_SET_2.out.cram
            ch_crai = BOOTSTRAP_VARIANT_SET_2.out.crai
            ch_vcf  = BOOTSTRAP_VARIANT_SET_2.out.vcf
            ch_tbi  = BOOTSTRAP_VARIANT_SET_2.out.tbi
        }

        // Remove bootstrapping metadata from CRAM channel
        ch_cram.map { meta, cram_file ->
            tuple( meta - meta.subMap('bootstrapping_round'), cram_file ) }
            .set { ch_cram }
        ch_crai.map { meta, crai_file ->
            tuple( meta - meta.subMap('bootstrapping_round'), crai_file ) }
            .set { ch_crai }
    } else {
        // TODO: Load existing variant set
        //ch_cram = channel.empty()
        //ch_crai = channel.empty()
        //ch_vcf = channel.empty()
        //ch_tbi = channel.empty()
    }

    if(params.stages.contains('base_quality_score_recalibration')) {
        //
        // SUBWORKFLOW: BASE_QUALITY_SCORE_RECALIBRATION
        //
        BASE_QUALITY_SCORE_RECALIBRATION(
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_intervals_split,
            ch_cram,
            ch_crai,
            ch_vcf,
            ch_tbi
        )
        ch_versions = ch_versions.mix(BASE_QUALITY_SCORE_RECALIBRATION.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(BASE_QUALITY_SCORE_RECALIBRATION.out.multiqc_files)
        ch_cram = BASE_QUALITY_SCORE_RECALIBRATION.out.cram
        ch_crai = BASE_QUALITY_SCORE_RECALIBRATION.out.crai
    }

    if(params.stages.contains('variant_calling')) {
        //
        // SUBWORKFLOW: CALL_VARIANTS_GATK
        //
        CALL_VARIANTS_GATK(
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_intervals_split,
            ch_cram,
            ch_crai
        )
        ch_versions = ch_versions.mix(CALL_VARIANTS_GATK.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(CALL_VARIANTS_GATK.out.multiqc_files)

        //
        // SUBWORKFLOW: CALL_VARIANTS_BCFTOOLS
        //
        CALL_VARIANTS_BCFTOOLS(
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_intervals_split,
            ch_cram,
            ch_crai
        )
        ch_versions = ch_versions.mix(CALL_VARIANTS_BCFTOOLS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(CALL_VARIANTS_BCFTOOLS.out.multiqc_files)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'swgsrelate_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
