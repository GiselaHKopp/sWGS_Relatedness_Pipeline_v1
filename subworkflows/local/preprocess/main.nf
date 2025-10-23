/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BWAMEM2_INDEX                  } from '../../../modules/nf-core/bwamem2/index'
include { BWAMEM2_MEM                    } from '../../../modules/nf-core/bwamem2/mem'
include { FASTP                          } from '../../../modules/nf-core/fastp'
include { GATK4_ADDORREPLACEREADGROUPS   } from '../../../modules/nf-core/gatk4/addorreplacereadgroups'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { GATK4SPARK_MARKDUPLICATES      } from '../../../modules/nf-core/gatk4spark/markduplicates'
include { MOSDEPTH                       } from '../../../modules/nf-core/mosdepth'
include { PRESEQ_CCURVE                  } from '../../../modules/nf-core/preseq/ccurve'
include { PRESEQ_LCEXTRAP                } from '../../../modules/nf-core/preseq/lcextrap'
include { SAMTOOLS_FAIDX                 } from '../../../modules/nf-core/samtools/faidx'
include { SAMTOOLS_STATS                 } from '../../../modules/nf-core/samtools/stats'
include { SPRING_DECOMPRESS              } from '../../../modules/nf-core/spring/decompress'

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
    ch_multiqc_files = Channel.empty()

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
    ch_multiqc_files = fastp_results.html.map { meta, file -> file }.mix(fastp_results.json.map { meta, file -> file })

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
    rg_bams = GATK4_ADDORREPLACEREADGROUPS(bwa_results.bam, reference_fasta, faidx_result.fai)
    ch_versions = ch_versions.mix(rg_bams.versions)

    // Mark duplicates
    seq_dict = GATK4_CREATESEQUENCEDICTIONARY(reference_fasta)
    ch_versions = ch_versions.mix(seq_dict.versions)

    // Mark duplicates
    markdup_results = GATK4SPARK_MARKDUPLICATES(rg_bams.bam, reference_fasta.map { it[1] }, faidx_result.fai.map{ it[1] }, seq_dict.dict.map{ it[1] })
    ch_versions = ch_versions.mix(markdup_results.versions)
    ch_multiqc_files = ch_multiqc_files.mix(markdup_results.metrics.map { it[1] })

    // Preseq analyses
    preseq_c_curve = PRESEQ_CCURVE(markdup_results.output)
    ch_versions = ch_versions.mix(preseq_c_curve.versions)
    ch_multiqc_files = ch_multiqc_files.mix(preseq_c_curve.c_curve.map { meta, file -> file }).mix(preseq_c_curve.log.map{ meta, file -> file })

    preseq_lc_extrap = PRESEQ_LCEXTRAP(markdup_results.output)
    ch_versions = ch_versions.mix(preseq_lc_extrap.versions)
    ch_multiqc_files = ch_multiqc_files.mix(preseq_lc_extrap.lc_extrap.map { meta, file -> file }).mix(preseq_lc_extrap.log.map{ meta, file -> file })

    // Samtools stats on final BAMs
    samstats_input = markdup_results.output.join(markdup_results.bam_index).map { meta, bam, bai -> tuple(meta, bam, bai) }
    samstats_results = SAMTOOLS_STATS(samstats_input, reference_fasta)
    ch_versions = ch_versions.mix(samstats_results.versions)
    ch_multiqc_files = ch_multiqc_files.mix(samstats_results.stats.map { it[1] })

    // Coverage calculation with mosdepth
    mosdepth_input = markdup_results.output.join(markdup_results.bam_index).map { meta, bam, bai -> tuple(meta, bam, bai, []) }
    mosdepth_results = MOSDEPTH(mosdepth_input, reference_fasta)
    ch_versions = ch_versions.mix(mosdepth_results.versions)
    ch_multiqc_files = ch_multiqc_files.mix(mosdepth_results.global_txt.map { meta, file -> file }).mix(mosdepth_results.summary_txt.map { meta, file -> file })

    emit:
    bams = markdup_results.output
    ch_multiqc_files
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
