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

include { BCFTOOLS_INDEX                                   } from '../modules/nf-core/bcftools/index/main'

include { BASE_QUALITY_SCORE_RECALIBRATION                 } from '../subworkflows/local/base_quality_score_recalibration'
include { BOOTSTRAP_VARIANT_SET as BOOTSTRAP_VARIANT_SET_1 } from '../subworkflows/local/bootstrap_variant_set'
include { BOOTSTRAP_VARIANT_SET as BOOTSTRAP_VARIANT_SET_2 } from '../subworkflows/local/bootstrap_variant_set'
include { BOOTSTRAP_VARIANT_SET as BOOTSTRAP_VARIANT_SET_3 } from '../subworkflows/local/bootstrap_variant_set'
include { CALL_VARIANTS_BCFTOOLS                           } from '../subworkflows/local/call_variants_bcftools'
include { CALL_VARIANTS_GATK                               } from '../subworkflows/local/call_variants_gatk'
include { PREPARE_GENOME                                   } from '../subworkflows/local/prepare_genome'
include { PREPARE_INTERVALS                                } from '../subworkflows/local/prepare_intervals'
include { PREPROCESS                                       } from '../subworkflows/local/preprocess'
include { RELATEDNESS_READ                                 } from '../subworkflows/local/relatedness_read'
include { VCF_INTERSECTION_THINNING                        } from '../subworkflows/local/vcf_intersection_thinning'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SWGSRELATE {

    take:
    samplesheet // channel: [ meta, list(fastq) ]

    main:
    ch_versions = channel.empty()
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
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // Gather built indices or get them from the params
    // Built from the fasta file:
    ch_dict = params.dict
        ? channel.fromPath(params.dict).map { it -> [[id: it.baseName], it] }.collect()
        : PREPARE_GENOME.out.dict
    ch_fasta_fai = params.fasta_fai
        ? channel.fromPath(params.fasta_fai).map { it -> [[id: it.baseName], it] }.collect()
        : PREPARE_GENOME.out.fasta_fai
    ch_bwamem2 = params.bwamem2_index
        ? channel.fromPath(params.bwamem2_index).map { it -> [[id: it.baseName], it] }.collect()
        : PREPARE_GENOME.out.bwamem2_index

    //
    // SUBWORKFLOW: PREPARE_INTERVALS
    //
    PREPARE_INTERVALS(ch_fasta_fai)
    ch_intervals_split = PREPARE_INTERVALS.out.intervals_split
    ch_versions = ch_versions.mix(PREPARE_INTERVALS.out.versions)

    //
    // SUBWORKFLOW: PREPROCESS
    //
    ch_preprocessed = PREPROCESS(samplesheet, ch_fasta, ch_bwamem2, ch_fasta_fai)
    ch_cram = ch_preprocessed.cram
    ch_crai = ch_preprocessed.crai
    ch_versions = ch_versions.mix(ch_preprocessed.versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_preprocessed.multiqc_files)

    //
    // SUBWORKFLOW: BOOTSTRAP_VARIANT_SET - ROUND 1
    //
    BOOTSTRAP_VARIANT_SET_1(
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_intervals_split,
        ch_cram,
        ch_crai,
        1
    )
    ch_versions = ch_versions.mix(BOOTSTRAP_VARIANT_SET_1.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(BOOTSTRAP_VARIANT_SET_1.out.multiqc_files)

    //
    // SUBWORKFLOW: BOOTSTRAP_VARIANT_SET - ROUND 2
    //
    BOOTSTRAP_VARIANT_SET_2(
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_intervals_split,
        BOOTSTRAP_VARIANT_SET_1.out.cram,
        BOOTSTRAP_VARIANT_SET_1.out.crai,
        2
    )

    //
    // SUBWORKFLOW: BOOTSTRAP_VARIANT_SET - ROUND 3
    //
    BOOTSTRAP_VARIANT_SET_3(
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_intervals_split,
        BOOTSTRAP_VARIANT_SET_2.out.cram,
        BOOTSTRAP_VARIANT_SET_2.out.crai,
        3
    )

    // Select which set of CRAM/VCF to use based on params.bootstrapping_rounds
    def cram_channels = [
        1: BOOTSTRAP_VARIANT_SET_1.out.cram,
        2: BOOTSTRAP_VARIANT_SET_2.out.cram,
        3: BOOTSTRAP_VARIANT_SET_3.out.cram
    ]
    def crai_channels = [
        1: BOOTSTRAP_VARIANT_SET_1.out.crai,
        2: BOOTSTRAP_VARIANT_SET_2.out.crai,
        3: BOOTSTRAP_VARIANT_SET_3.out.crai
    ]
    def vcf_channels = [
        1: BOOTSTRAP_VARIANT_SET_1.out.vcf,
        2: BOOTSTRAP_VARIANT_SET_2.out.vcf,
        3: BOOTSTRAP_VARIANT_SET_3.out.vcf
    ]
    def tbi_channels = [
        1: BOOTSTRAP_VARIANT_SET_1.out.tbi,
        2: BOOTSTRAP_VARIANT_SET_2.out.tbi,
        3: BOOTSTRAP_VARIANT_SET_3.out.tbi
    ]

    ch_cram = params.known_variants_vcf ? ch_cram : cram_channels[ params.bootstrapping_rounds ]
    ch_crai = params.known_variants_vcf ? ch_crai : crai_channels[ params.bootstrapping_rounds ]
    ch_vcf  = params.known_variants_vcf ? channel.fromPath(params.known_variants_vcf)
        .map { it -> [[id: 'known_variants_vcf'], it] }
        .collect()
        :
        vcf_channels[ params.bootstrapping_rounds ]
    ch_tbi = params.known_variants_vcf ? (params.known_variants_tbi ? channel.fromPath(params.known_variants_tbi)
        .map { tbi -> [ [id: 'known_variants_tbi'], tbi ] }.collect()
        : BCFTOOLS_INDEX(ch_vcf).out.tbi.collect())
        :
        tbi_channels[ params.bootstrapping_rounds ]

    //
    // SUBWORKFLOW: BASE_QUALITY_SCORE_RECALIBRATION
    //
    BASE_QUALITY_SCORE_RECALIBRATION(
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_intervals_split,
        ch_cram,
        ch_crai,
        ch_vcf,
        ch_tbi
    )
    ch_versions = ch_versions.mix(BASE_QUALITY_SCORE_RECALIBRATION.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(BASE_QUALITY_SCORE_RECALIBRATION.out.multiqc_files)
    ch_cram = BASE_QUALITY_SCORE_RECALIBRATION.out.recalibrated_cram
    ch_crai = BASE_QUALITY_SCORE_RECALIBRATION.out.recalibrated_crai

    //
    // SUBWORKFLOW: CALL_VARIANTS_GATK
    //
    CALL_VARIANTS_GATK(
        ch_fasta,
        ch_fasta_fai,
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
        ch_fasta_fai,
        ch_dict,
        ch_intervals_split,
        ch_cram,
        ch_crai
    )
    ch_versions = ch_versions.mix(CALL_VARIANTS_BCFTOOLS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(CALL_VARIANTS_BCFTOOLS.out.multiqc_files)

    //
    // SUBWORKFLOW: VCF_INTERSECTION
    //
    VCF_INTERSECTION_THINNING(
        CALL_VARIANTS_GATK.out.vcf,
        CALL_VARIANTS_GATK.out.tbi,
        CALL_VARIANTS_BCFTOOLS.out.vcf,
        CALL_VARIANTS_BCFTOOLS.out.tbi,
        PREPARE_INTERVALS.out.intervals_combined
    )
    ch_versions = ch_versions.mix(VCF_INTERSECTION_THINNING.out.versions)

    //
    // SUBWORKFLOW: RELATEDNESS_READ
    //
    RELATEDNESS_READ(VCF_INTERSECTION_THINNING.out.intersection)
    ch_versions = ch_versions.mix(RELATEDNESS_READ.out.versions)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
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
