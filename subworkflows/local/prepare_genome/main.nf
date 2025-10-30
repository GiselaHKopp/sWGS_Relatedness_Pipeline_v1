/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { SAMTOOLS_FAIDX                 } from '../../../modules/nf-core/samtools/faidx'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PREPARE_GENOME {
    take:
    fasta

    main:
    versions = channel.empty()

    // Build the FASTA index (fai)
    SAMTOOLS_FAIDX(fasta, [[id: 'no_fai'], []], false)
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)

    // Build the sequence dictionary
    GATK4_CREATESEQUENCEDICTIONARY(fasta)
    versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

    emit:
    fai = SAMTOOLS_FAIDX.out.fai
    dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
    versions
}
