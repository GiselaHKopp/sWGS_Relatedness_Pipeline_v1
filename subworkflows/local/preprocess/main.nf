/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BWAMEM2_MEM                    } from '../../../modules/nf-core/bwamem2/mem'
include { FASTP                          } from '../../../modules/nf-core/fastp'
include { GATK4_ADDORREPLACEREADGROUPS   } from '../../../modules/nf-core/gatk4/addorreplacereadgroups'
include { GATK4_MARKDUPLICATES           } from '../../../modules/nf-core/gatk4/markduplicates'
include { MOSDEPTH                       } from '../../../modules/nf-core/mosdepth'
include { PRESEQ_CCURVE                  } from '../../../modules/nf-core/preseq/ccurve'
include { PRESEQ_LCEXTRAP                } from '../../../modules/nf-core/preseq/lcextrap'
include { SAMTOOLS_INDEX                 } from '../../../modules/nf-core/samtools/index'
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
    bwamem2     // channel: [ meta, bwamem2 ]
    fai         // channel: [ meta, fai]

    main:
    versions = channel.empty()
    multiqc_files = channel.empty()

    // Split by file type (spring vs fastq)
    samplesheet.branch { row ->
        spring: row[1].every { file -> file.getName().endsWith('.spring') }
        fastq : row[1].every { file -> file.getName().endsWith('.fastq') || file.getName().endsWith('.fastq.gz') || file.getName().endsWith('.fq.gz') }
        bam   : row[1].every { file -> file.getName().endsWith('.bam') }
        cram  : row[1].every { file -> file.getName().endsWith('.cram') }
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

    // Map to reference
    BWAMEM2_MEM(FASTP.out.reads, bwamem2, fasta, true)
    versions = versions.mix(BWAMEM2_MEM.out.versions)
    bam = BWAMEM2_MEM.out.bam.mix(ch_input_branches.bam)

    // Add read groups
    GATK4_ADDORREPLACEREADGROUPS(bam, fasta, fai)
    versions = versions.mix(GATK4_ADDORREPLACEREADGROUPS.out.versions)

    // Mark duplicates
    GATK4_MARKDUPLICATES(GATK4_ADDORREPLACEREADGROUPS.out.bam, fasta.map { tuple -> tuple[1] }, fai.map{ tuple -> tuple[1] })
    versions = versions.mix(GATK4_MARKDUPLICATES.out.versions)
    multiqc_files = multiqc_files.mix(GATK4_MARKDUPLICATES.out.metrics.map { tuple -> tuple[1] })
    ch_cram = GATK4_MARKDUPLICATES.out.cram.mix(ch_input_branches.cram)

    // Compute index
    SAMTOOLS_INDEX(ch_input_branches.cram)
    versions = versions.mix(SAMTOOLS_INDEX.out.versions)
    ch_crai = GATK4_MARKDUPLICATES.out.crai.mix(SAMTOOLS_INDEX.out.crai)
/*
    // Preseq analyses
    PRESEQ_CCURVE(ch_cram)
    versions = versions.mix(PRESEQ_CCURVE.out.versions)
    multiqc_files = multiqc_files.mix(PRESEQ_CCURVE.out.c_curve.map { _meta, file -> file }).mix(PRESEQ_CCURVE.out.log.map{ _meta, file -> file })

    PRESEQ_LCEXTRAP(ch_cram)
    versions = versions.mix(PRESEQ_LCEXTRAP.out.versions)
    multiqc_files = multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap.map { _meta, file -> file }).mix(PRESEQ_LCEXTRAP.out.log.map{ _meta, file -> file })
*/
    // Samtools stats on final CRAMs
    ch_samstats_input = ch_cram.join(ch_crai).map { meta, cram, crai -> tuple(meta, cram, crai) }
    SAMTOOLS_STATS(ch_samstats_input, fasta)
    multiqc_files = multiqc_files.mix(SAMTOOLS_STATS.out.stats.map { tuple -> tuple[1] })

    // Coverage calculation with mosdepth
    ch_mosdepth_input = ch_cram.join(ch_crai).map { meta, cram, crai -> tuple(meta, cram, crai, []) }
    MOSDEPTH(ch_mosdepth_input, fasta)
    versions = versions.mix(MOSDEPTH.out.versions)
    multiqc_files = multiqc_files.mix(MOSDEPTH.out.global_txt.map { _meta, file -> file }).mix(MOSDEPTH.out.summary_txt.map { _meta, file -> file })

    emit:
    cram = ch_cram
    crai = ch_crai
    multiqc_files
    versions
}
