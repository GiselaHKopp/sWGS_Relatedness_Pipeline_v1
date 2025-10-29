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
include { GATK4_MARKDUPLICATES           } from '../../../modules/nf-core/gatk4/markduplicates'
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
    versions = channel.empty()
    multiqc_files = channel.empty()

    // Split by file type (spring vs fastq)
    samplesheet.branch { row ->
        spring: row[1].every { file -> file.getName().endsWith('.spring') }
        fastq : row[1].every { file -> file.getName().endsWith('.fastq') || file.getName().endsWith('.fastq.gz') || file.getName().endsWith('.fq.gz') }
    }.set{ input_branches }

    // Decompress SPRING → FASTQ pairs
    spring_pairs = SPRING_DECOMPRESS(input_branches.spring, false)
    versions = versions.mix(spring_pairs.versions)

    // Merge with normal FASTQs into one unified channel
    merged_fastqs = input_branches.fastq.mix(spring_pairs.fastq)

    // Trim & QC with FASTP
    fastp_input = merged_fastqs.map { meta, reads -> tuple(meta, reads, []) }
    fastp_results = FASTP(fastp_input, false, false, false)
    versions = versions.mix(fastp_results.versions)
    multiqc_files = fastp_results.html.map { _meta, file -> file }.mix(fastp_results.json.map { _meta, file -> file })

    // Build the BWA index from the provided FASTA
    bwa_index = BWAMEM2_INDEX(reference_fasta)
    versions = versions.mix(bwa_index.versions)

    // Map to reference
    bwa_results = BWAMEM2_MEM(fastp_results.reads, bwa_index.index, reference_fasta, true)
    versions = versions.mix(bwa_results.versions)

    // Build the FASTA index (fai)
    faidx_result = SAMTOOLS_FAIDX(reference_fasta, [[id: 'no_fai'], []], false)
    versions = versions.mix(faidx_result.versions)

    // Add read groups
    rg_bams = GATK4_ADDORREPLACEREADGROUPS(bwa_results.bam, reference_fasta, faidx_result.fai)
    versions = versions.mix(rg_bams.versions)

    // Mark duplicates
    seq_dict = GATK4_CREATESEQUENCEDICTIONARY(reference_fasta)
    versions = versions.mix(seq_dict.versions)

    // Mark duplicates
    if(params.use_gatk_spark) {
        markduplicates_results = GATK4SPARK_MARKDUPLICATES(rg_bams.bam, reference_fasta.map { tuple -> tuple[1] }, faidx_result.fai.map{ tuple -> tuple[1] }, seq_dict.dict.map{ tuple -> tuple[1] })
        ch_markduplicates_bam = markduplicates_results.output
        ch_markduplicates_bai = markduplicates_results.bam_index
    } else {
        markduplicates_results = GATK4_MARKDUPLICATES(rg_bams.bam, reference_fasta.map { tuple -> tuple[1] }, faidx_result.fai.map{ tuple -> tuple[1] })
        ch_markduplicates_bam = markduplicates_results.bam
        ch_markduplicates_bai = markduplicates_results.bai
    }
    versions = versions.mix(markduplicates_results.versions)
    multiqc_files = multiqc_files.mix(markduplicates_results.metrics.map { tuple -> tuple[1] })

    // Preseq analyses
    preseq_c_curve = PRESEQ_CCURVE(ch_markduplicates_bam)
    versions = versions.mix(preseq_c_curve.versions)
    multiqc_files = multiqc_files.mix(preseq_c_curve.c_curve.map { _meta, file -> file }).mix(preseq_c_curve.log.map{ _meta, file -> file })

    preseq_lc_extrap = PRESEQ_LCEXTRAP(ch_markduplicates_bam)
    versions = versions.mix(preseq_lc_extrap.versions)
    multiqc_files = multiqc_files.mix(preseq_lc_extrap.lc_extrap.map { _meta, file -> file }).mix(preseq_lc_extrap.log.map{ _meta, file -> file })

    // Samtools stats on final BAMs
    samstats_input = ch_markduplicates_bam.join(ch_markduplicates_bai).map { meta, bam, bai -> tuple(meta, bam, bai) }
    samstats_results = SAMTOOLS_STATS(samstats_input, reference_fasta)
    versions = versions.mix(samstats_results.versions)
    multiqc_files = multiqc_files.mix(samstats_results.stats.map { tuple -> tuple[1] })

    // Coverage calculation with mosdepth
    mosdepth_input = ch_markduplicates_bam.join(ch_markduplicates_bai).map { meta, bam, bai -> tuple(meta, bam, bai, []) }
    mosdepth_results = MOSDEPTH(mosdepth_input, reference_fasta)
    versions = versions.mix(mosdepth_results.versions)
    multiqc_files = multiqc_files.mix(mosdepth_results.global_txt.map { _meta, file -> file }).mix(mosdepth_results.summary_txt.map { _meta, file -> file })

    emit:
    bam = ch_markduplicates_bam        // val(meta), path(bam)
    bam_index = ch_markduplicates_bai  // val(meta), path(bai)
    fai = faidx_result.fai             // val(meta), path(fai)
    dict = seq_dict.dict               // val(meta), path(dict)
    multiqc_files
    versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
