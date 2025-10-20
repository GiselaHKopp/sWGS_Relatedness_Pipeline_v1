/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BWAMEM2_INDEX                } from '../../../modules/nf-core/bwamem2/index'
include { BWAMEM2_MEM                  } from '../../../modules/nf-core/bwamem2/mem'
include { FASTP                        } from '../../../modules/nf-core/fastp'
include { GATK4_ADDORREPLACEREADGROUPS } from '../../../modules/nf-core/gatk4/addorreplacereadgroups'
include { SAMTOOLS_FAIDX               } from '../../../modules/nf-core/samtools/faidx'
include { SPRING_DECOMPRESS            } from '../../../modules/nf-core/spring/decompress'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PREPROCESS {
    take:
    samplesheet
    reference_fasta

    main:
    // Collect software versions and QC reports
    ch_versions = Channel.empty()
    qc_reports = Channel.empty()

    // Split by file type (spring vs fastq)
    samplesheet.branch {
        spring: it[1].every { f -> f.getName().endsWith('.spring') }
        fastq : it[1].every { f -> f.getName().endsWith('.fastq') || f.getName().endsWith('.fastq.gz') || f.getName().endsWith('.fq.gz') }
    }.set{ input_branches }

    // Decompress SPRING → FASTQ pairs
    spring_pairs = SPRING_DECOMPRESS(input_branches.spring, false)
    ch_versions = ch_versions.mix(spring_pairs.versions)

    // Merge with normal FASTQs into one unified channel
    merged_fastqs = input_branches.fastq.mix(spring_pairs.fastq)

    // Trim & QC with FASTP
    fastp_input = merged_fastqs.map { meta, reads -> tuple(meta, reads, []) }
    fastp_results = FASTP(fastp_input, false, false, false)
    ch_versions = ch_versions.mix(fastp_results.versions)
    qc_reports = fastp_results.html
                  .map { meta, file -> file }
                  .mix(fastp_results.json.map { meta, file -> file })

    // Build the BWA index from the provided FASTA
    bwa_index = BWAMEM2_INDEX(reference_fasta)
    ch_versions = ch_versions.mix(bwa_index.versions)

    // Map to reference
    bwa_results = BWAMEM2_MEM(fastp_results.reads, bwa_index.index, reference_fasta, true)
    ch_versions = ch_versions.mix(bwa_results.versions)

    // Build the FASTA index (fai)
    faidx_result = SAMTOOLS_FAIDX(reference_fasta, [[id: 'no_fai'], []], false)
    ch_versions = ch_versions.mix(faidx_result.versions)

    // Add read groups
    bwa_results.bam.view { "BWA emits: ${it} (${it.getClass()})" }
    rg_bams = GATK4_ADDORREPLACEREADGROUPS(bwa_results.bam, bwa_index.index, reference_fasta)
    ch_versions = ch_versions.mix(rg_bams.versions)

    emit:
    bams = rg_bams.bam
    qc_reports
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
