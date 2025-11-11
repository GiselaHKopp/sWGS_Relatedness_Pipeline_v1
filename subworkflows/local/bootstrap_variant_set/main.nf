/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BASE_QUALITY_SCORE_RECALIBRATION } from '../base_quality_score_recalibration'
include { CALL_VARIANTS_GATK               } from '../call_variants_gatk'
include { FILTER_VARIANTS                  } from '../filter_variants'

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
    // SUBWORKFLOW: CALL_VARIANTS_GATK
    //
    CALL_VARIANTS_GATK(
        fasta,
        fai,
        dict,
        intervals,
        cram,
        crai
    )
    versions = versions.mix(CALL_VARIANTS_GATK.out.versions)
    multiqc_files = multiqc_files.mix(CALL_VARIANTS_GATK.out.multiqc_files)

    //
    // SUBWORKFLOW: FILTER_VARIANTS
    //
    FILTER_VARIANTS(
        fasta,
        fai,
        dict,
        CALL_VARIANTS_GATK.out.vcf,
        CALL_VARIANTS_GATK.out.tbi
    )
    versions = versions.mix(FILTER_VARIANTS.out.versions)
    multiqc_files = multiqc_files.mix(FILTER_VARIANTS.out.multiqc_files)

    //
    // SUBWORKFLOW: BASE_QUALITY_SCORE_RECALIBRATION
    //
    BASE_QUALITY_SCORE_RECALIBRATION(
        fasta,
        fai,
        dict,
        intervals,
        cram,
        crai,
        FILTER_VARIANTS.out.vcf,
        FILTER_VARIANTS.out.tbi
    )
    versions = versions.mix(BASE_QUALITY_SCORE_RECALIBRATION.out.versions)
    multiqc_files = multiqc_files.mix(BASE_QUALITY_SCORE_RECALIBRATION.out.multiqc_files)

    emit:
    cram = BASE_QUALITY_SCORE_RECALIBRATION.out.recalibrated_cram
    crai = BASE_QUALITY_SCORE_RECALIBRATION.out.recalibrated_crai
    vcf = FILTER_VARIANTS.out.vcf
    tbi = FILTER_VARIANTS.out.tbi
    multiqc_files
    versions
}
