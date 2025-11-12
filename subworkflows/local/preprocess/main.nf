/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BWAMEM2_INDEX                  } from '../../../modules/nf-core/bwamem2/index'
include { BWAMEM2_MEM                    } from '../../../modules/nf-core/bwamem2/mem'
include { FASTP                          } from '../../../modules/nf-core/fastp'
include { GATK4_ADDORREPLACEREADGROUPS   } from '../../../modules/nf-core/gatk4/addorreplacereadgroups'
include { GATK4_MARKDUPLICATES           } from '../../../modules/nf-core/gatk4/markduplicates'
include { MOSDEPTH                       } from '../../../modules/nf-core/mosdepth'
include { PRESEQ_CCURVE                  } from '../../../modules/nf-core/preseq/ccurve'
include { PRESEQ_LCEXTRAP                } from '../../../modules/nf-core/preseq/lcextrap'
include { SAMTOOLS_STATS                 } from '../../../modules/nf-core/samtools/stats'
include { SPRING_DECOMPRESS              } from '../../../modules/nf-core/spring/decompress'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PREPROCESS {
    take:
    samplesheet // channel: [ meta, list(fastq) ]
    fasta       // channel: [ meta, fasta]
    fai         // channel: [ meta, fai]

    main:
    versions = channel.empty()
    multiqc_files = channel.empty()

    // Split by file type (spring vs fastq)
    samplesheet.branch { row ->
        spring: row[1].every { file -> file.getName().endsWith('.spring') }
        fastq : row[1].every { file -> file.getName().endsWith('.fastq') || file.getName().endsWith('.fastq.gz') || file.getName().endsWith('.fq.gz') }
    }.set{ ch_input_branches }

    // Decompress SPRING to FASTQ pairs
    SPRING_DECOMPRESS(ch_input_branches.spring, false)
    versions = versions.mix(SPRING_DECOMPRESS.out.versions)

    // Merge with normal FASTQs into one unified channel
    merged_fastqs = ch_input_branches.fastq.mix(SPRING_DECOMPRESS.out.fastq)

    // Trim and QC with FASTP
    ch_fastp_input = merged_fastqs.map { meta, reads -> tuple(meta, reads, []) }
    FASTP(ch_fastp_input, false, false, false)
    versions = versions.mix(FASTP.out.versions)
    multiqc_files = FASTP.out.html.map { _meta, file -> file }.mix(FASTP.out.json.map { _meta, file -> file })

    // Build the BWA index from the provided FASTA
    BWAMEM2_INDEX(fasta)
    versions = versions.mix(BWAMEM2_INDEX.out.versions)

    // Map to reference
    BWAMEM2_MEM(FASTP.out.reads, BWAMEM2_INDEX.out.index, fasta, true)
    versions = versions.mix(BWAMEM2_MEM.out.versions)

    // Add read groups
    GATK4_ADDORREPLACEREADGROUPS(BWAMEM2_MEM.out.bam, fasta, fai)
    versions = versions.mix(GATK4_ADDORREPLACEREADGROUPS.out.versions)

    // Mark duplicates
    markduplicates_results = GATK4_MARKDUPLICATES(GATK4_ADDORREPLACEREADGROUPS.out.bam, fasta.map { tuple -> tuple[1] }, fai.map{ tuple -> tuple[1] })
    versions = versions.mix(markduplicates_results.versions)
    multiqc_files = multiqc_files.mix(GATK4_MARKDUPLICATES.out.metrics.map { tuple -> tuple[1] })
/*
    // Preseq analyses
    PRESEQ_CCURVE(GATK4_MARKDUPLICATES.out.cram)
    versions = versions.mix(PRESEQ_CCURVE.out.versions)
    multiqc_files = multiqc_files.mix(PRESEQ_CCURVE.out.c_curve.map { _meta, file -> file }).mix(PRESEQ_CCURVE.out.log.map{ _meta, file -> file })

    PRESEQ_LCEXTRAP(GATK4_MARKDUPLICATES.out.cram)
    versions = versions.mix(PRESEQ_LCEXTRAP.out.versions)
    multiqc_files = multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap.map { _meta, file -> file }).mix(PRESEQ_LCEXTRAP.out.log.map{ _meta, file -> file })
*/
    // Samtools stats on final CRAMs
    ch_samstats_input = GATK4_MARKDUPLICATES.out.cram.join(GATK4_MARKDUPLICATES.out.crai).map { meta, cram, crai -> tuple(meta, cram, crai) }
    SAMTOOLS_STATS(ch_samstats_input, fasta)
    versions = versions.mix(SAMTOOLS_STATS.out.versions)
    multiqc_files = multiqc_files.mix(SAMTOOLS_STATS.out.stats.map { tuple -> tuple[1] })

    // Coverage calculation with mosdepth
    ch_mosdepth_input = GATK4_MARKDUPLICATES.out.cram.join(GATK4_MARKDUPLICATES.out.crai).map { meta, cram, crai -> tuple(meta, cram, crai, []) }
    MOSDEPTH(ch_mosdepth_input, fasta)
    versions = versions.mix(MOSDEPTH.out.versions)
    multiqc_files = multiqc_files.mix(MOSDEPTH.out.global_txt.map { _meta, file -> file }).mix(MOSDEPTH.out.summary_txt.map { _meta, file -> file })

    emit:
    cram = GATK4_MARKDUPLICATES.out.cram
    crai = GATK4_MARKDUPLICATES.out.crai
    multiqc_files
    versions
}
