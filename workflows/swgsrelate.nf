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

include { PREPROCESS             } from '../subworkflows/local/preprocess'
include { CALL_VARIANTS          } from '../subworkflows/local/call_variants'
include { FILTER_VARIANTS        } from '../subworkflows/local/filter_variants'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SWGSRELATE {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:
    // Channel for collecting software versions
    ch_versions = channel.empty()

    // Channel for collecting MultiQC files
    ch_multiqc_files = channel.empty()

    // Define reference genome and index
    ch_reference_fasta = params.fasta ?
    channel.fromPath(params.fasta)
        .map { f -> [ [id: f.baseName], f ] }
        .collect()
    : channel.empty()


    if(params.stages.contains('preprocess')) {
        //
        // SUBWORKFLOW: PREPROCESS
        //
        ch_preprocessed = PREPROCESS(samplesheet, ch_reference_fasta)
        ch_versions = ch_versions.mix(ch_preprocessed.versions)
        ch_multiqc_files = ch_multiqc_files.mix(ch_preprocessed.multiqc_files)
    } else {
        // TODO: Load BAMs from samplesheet without preprocessing
        ch_preprocessed = [
            bam: channel.empty(),
            bam_index: channel.empty(),
            fai: channel.empty(),
            dict: channel.empty(),
            versions: channel.empty(),
            multiqc_files: channel.empty()
        ]
    }

    if(params.stages.contains('prepare_variant_set')) {
        //
        // SUBWORKFLOW: CALL_VARIANTS
        //
        CALL_VARIANTS(
            ch_preprocessed.bam,
            ch_preprocessed.bam_index,
            ch_preprocessed.fai,
            ch_preprocessed.dict,
            ch_reference_fasta
        )
        ch_versions = ch_versions.mix(CALL_VARIANTS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(CALL_VARIANTS.out.multiqc_files)

        //
        // SUBWORKFLOW: FILTER_VARIANTS
        //
        FILTER_VARIANTS(
            CALL_VARIANTS.out.vcf,
            CALL_VARIANTS.out.tbi,
            ch_preprocessed.fai,
            ch_preprocessed.dict,
            ch_reference_fasta
        )
        ch_versions = ch_versions.mix(FILTER_VARIANTS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(FILTER_VARIANTS.out.multiqc_files)
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
