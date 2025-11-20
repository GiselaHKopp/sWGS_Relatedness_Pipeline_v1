/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BWAMEM2_INDEX                  } from '../../../modules/nf-core/bwamem2/index'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { SAMTOOLS_FAIDX                 } from '../../../modules/nf-core/samtools/faidx'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PREPARE_GENOME {
    take:
    fasta   // channel: [ meta, fasta]

    main:
    versions = channel.empty()

    // Build the BWA index from the provided FASTA
    BWAMEM2_INDEX(fasta)
    versions = versions.mix(BWAMEM2_INDEX.out.versions)

    // Build the sequence dictionary
    GATK4_CREATESEQUENCEDICTIONARY(fasta)
    versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

    // Build the FASTA index (fai)
    SAMTOOLS_FAIDX(fasta, [[id: 'no_fai'], []], false)
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)

    emit:
    bwamem2_index   = BWAMEM2_INDEX.out.index.collect()
    dict            = GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect()
    fasta_fai       = SAMTOOLS_FAIDX.out.fai.collect()
    versions
}
