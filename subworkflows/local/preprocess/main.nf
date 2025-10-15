/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BWA_MEM                      } from '../../../modules/nf-core/bwa/mem'
include { FASTQC                       } from '../../../modules/nf-core/fastqc'
include { GATK4_ADDORREPLACEREADGROUPS } from '../../../modules/nf-core/gatk4/addorreplacereadgroups'
include { PICARD_MERGESAMFILES         } from '../../../modules/nf-core/picard/mergesamfiles'
include { SAMTOOLS_SORT                } from '../../../modules/nf-core/samtools/sort'
include { TRIMMOMATIC                  } from '../../../modules/nf-core/trimmomatic'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PREPROCESS {
    take:
        samplesheet

    main:
    // Collect software versions
    ch_versions = Channel.empty()

    // Parse TSV and emit tuple(pair_id, path R1, path R2, read group info)
    parsed = samplesheet.splitCsv(header: true, sep: '\t')
        .map { row ->
            def rg = "--RGID ${row.RGID} --RGLB ${row.RGLB} --RGPL ${row.RGPL} --RGPU ${row.RGPU} --RGSM ${row.RGSM}"
            tuple(row.RGSM, file(row.file_name_one), file(row.file_name_two), rg)
        }

    // FastQC before trimming
    parsed.map { id, r1, r2, rg -> tuple(id+"_1", r1) }
          .concat(parsed.map { id, r1, r2, rg -> tuple(id+"_2", r2) })
          .set { raw_fastq }

    FASTQC(raw_fastq)
    FASTQC_RAW = FASTQC.out
    ch_versions = ch_versions.mix(FASTQC_RAW.versions.first())

    // Trimming
    parsed.map { id, r1, r2, rg -> tuple(id, [r1, r2]) } | TRIMMOMATIC
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions.first())

    // FastQC after trimming
    TRIMMOMATIC.out.trimmed_reads
        .map { id, reads -> tuple(id+"_1", reads[0]) }
        .concat(TRIMMOMATIC.out.trimmed_reads.map { id, reads -> tuple(id+"_2", reads[1]) })
        .set { trimmed_fastq }

    FASTQC(trimmed_fastq)
    FASTQC_TRIMMED = FASTQC.out

    // Combine both FASTQC results
    combined_qc = FASTQC_RAW.concat(FASTQC_TRIMMED)

    // Alignment
    TRIMMOMATIC.out.trimmed_reads
        .map { id, r1r2 -> tuple(id, r1r2[0], r1r2[1]) }
        .combine(parsed.map { id, r1, r2, rg -> tuple(id, rg) })
        .map { id, r1, r2, rg -> tuple(id, [r1, r2], rg) }
        | BWA_MEM
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    // Sort BAMs
    BWA_MEM.out.bam | SAMTOOLS_SORT
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    // Merge BAMs (if needed)
    SAMTOOLS_SORT.out | PICARD_MERGESAMFILES
    ch_versions = ch_versions.mix(PICARD_MERGESAMFILES.out.versions.first())


    PICARD_MERGESAMFILES.out
        .combine(parsed.map { id, r1, r2, rg -> tuple(id, rg) })
        .map { id, bam, rg -> tuple(id, bam, rg) }
        | GATK4_ADDORREPLACEREADGROUPS
    ch_versions = ch_versions.mix(GATK4_ADDORREPLACEREADGROUPS.out.versions.first())

    emit:
    sorted_bams = SAMTOOLS_SORT.out
    qc_reports = combined_qc
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
