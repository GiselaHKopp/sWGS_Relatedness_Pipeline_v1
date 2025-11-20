/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BASE_QUALITY_SCORE_RECALIBRATION as BQSR_BOOTSTRAP } from '../base_quality_score_recalibration'
include { CALL_VARIANTS_GATK as CALL_VARIANTS_GATK_BOOTSTRAP } from '../call_variants_gatk'
include { FILTER_VARIANTS as FILTER_VARIANTS_BOOTSTRAP       } from '../filter_variants'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow BOOTSTRAP_VARIANT_SET {
    take:
    fasta       // channel: [ meta, fasta]
    fai         // channel: [ meta, fai]
    dict        // channel: [ meta, dict]
    intervals   // channel: [ meta, intervals, number_of_intervals]
    cram        // channel: [ meta, cram]
    crai        // channel: [ meta, crai]

    main:
    versions = channel.empty()
    multiqc_files = channel.empty()

    //
    // SUBWORKFLOW: CALL_VARIANTS_GATK_BOOTSTRAP
    //
    CALL_VARIANTS_GATK_BOOTSTRAP(
        fasta,
        fai,
        dict,
        intervals,
        cram,
        crai
    )
    versions = versions.mix(CALL_VARIANTS_GATK_BOOTSTRAP.out.versions)
    multiqc_files = multiqc_files.mix(CALL_VARIANTS_GATK_BOOTSTRAP.out.multiqc_files)

    //
    // SUBWORKFLOW: FILTER_VARIANTS
    //
    FILTER_VARIANTS_BOOTSTRAP(
        fasta,
        fai,
        dict,
        CALL_VARIANTS_GATK_BOOTSTRAP.out.vcf,
        CALL_VARIANTS_GATK_BOOTSTRAP.out.tbi
    )
    versions = versions.mix(FILTER_VARIANTS_BOOTSTRAP.out.versions)
    multiqc_files = multiqc_files.mix(FILTER_VARIANTS_BOOTSTRAP.out.multiqc_files)

    //
    // SUBWORKFLOW: BQSR_BOOTSTRAP
    //
    BQSR_BOOTSTRAP(
        fasta,
        fai,
        dict,
        intervals,
        cram,
        crai,
        FILTER_VARIANTS_BOOTSTRAP.out.vcf,
        FILTER_VARIANTS_BOOTSTRAP.out.tbi
    )
    versions = versions.mix(BQSR_BOOTSTRAP.out.versions)
    multiqc_files = multiqc_files.mix(BQSR_BOOTSTRAP.out.multiqc_files)

    emit:
    cram = BQSR_BOOTSTRAP.out.recalibrated_cram
    crai = BQSR_BOOTSTRAP.out.recalibrated_crai
    vcf = FILTER_VARIANTS_BOOTSTRAP.out.vcf
    tbi = FILTER_VARIANTS_BOOTSTRAP.out.tbi
    multiqc_files
    versions
}
