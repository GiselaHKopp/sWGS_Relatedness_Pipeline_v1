/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//include { GATK4_ANALYZECOVARIATES                    } from '../../../modules/nf-core/gatk4/analyzecovariates'
include { GATK4_APPLYBQSR                            } from '../../../modules/nf-core/gatk4/applybqsr'
include { SAMTOOLS_INDEX                             } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_SCATTERED } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                             } from '../../../modules/nf-core/samtools/merge/main'

include { CRAM_BASERECALIBRATOR                                      } from '../cram_baserecalibrator'
include { CRAM_BASERECALIBRATOR as CRAM_BASERECALIBRATOR_SECOND_PASS } from '../cram_baserecalibrator'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow BASE_QUALITY_SCORE_RECALIBRATION {
    take:
    fasta
    fai
    dict
    intervals
    cram
    crai
    vcf
    tbi

    main:
    versions = channel.empty()
    multiqc_files = channel.empty()

    // Combine CRAM with intervals
    ch_cram = cram.join(crai)
    .combine(intervals)
    .map { cram_meta, cram_file, crai_file, interval_meta, interval_file, num_intervals ->
        // Construct new ID: sampleID_intervalName
        def new_id = "${cram_meta.id}_${interval_meta.interval_name}"

        // Merge metadata and overwrite id
        def meta = cram_meta + interval_meta + [ id: new_id ] + [ num_intervals: num_intervals ]

        tuple(meta, cram_file, crai_file, interval_file)
    }

    CRAM_BASERECALIBRATOR(fasta, fai, dict, ch_cram, vcf, tbi)
    versions = versions.mix(CRAM_BASERECALIBRATOR.out.versions)

    // Combine CRAM with BQSR table
    ch_cram_with_table = ch_cram
        .map { meta, cram_file, crai_file, interval_file ->
            def key = meta.RGSM
            tuple(key, meta, cram_file, crai_file, interval_file)
        }
        .join(
            CRAM_BASERECALIBRATOR.out.ch_table_bqsr.map { meta, table ->
                def key = meta.RGSM
                tuple(key, table)
            },
            by: 0
        )
        .map { _key, meta, cram_file, crai_file, interval_file, table ->
            tuple(meta, cram_file, crai_file, table, interval_file)
        }

    // Run ApplyBQSR
    GATK4_APPLYBQSR(
        ch_cram_with_table,
        fasta.map { _meta, fasta_file -> [fasta_file] },
        fai.map { _meta, fai_file -> [fai_file] },
        dict.map { _meta, dict_file -> [dict_file] },
    )
    versions = versions.mix(GATK4_APPLYBQSR.out.versions)

    // Index recalibrated CRAMs (per interval)
    SAMTOOLS_INDEX_SCATTERED(GATK4_APPLYBQSR.out.cram)
    versions = versions.mix(SAMTOOLS_INDEX_SCATTERED.out.versions)

    // Combine recalibrated CRAM with intervals for second pass
    ch_cram_second_pass = GATK4_APPLYBQSR.out.cram.join(SAMTOOLS_INDEX_SCATTERED.out.crai)
    .combine(intervals)
    .map { cram_meta, cram_file, crai_file, _interval_meta, interval_file, _num_intervals ->
        tuple(cram_meta, cram_file, crai_file, interval_file)
    }

    // Run BaseRecalibrator (second pass, for quality control)
    CRAM_BASERECALIBRATOR_SECOND_PASS(
        fasta,
        fai,
        dict,
        ch_cram_second_pass,
        vcf.map{ _meta, files -> [['id' : 'known_sites'], files]},
        tbi.map{ _meta, files -> [['id' : 'known_sites'], files]}
    )

    // Merge recalibrated CRAMs if needed
    ch_cram_branch = GATK4_APPLYBQSR.out.cram.map{ meta, table -> [ groupKey(meta, meta.num_intervals), table ] }.groupTuple()
        .branch { tuple ->
            single:   tuple[0].num_intervals == 1
            multiple: tuple[0].num_intervals > 1
        }

    // Merge CRAMs if multiple intervals
    SAMTOOLS_MERGE(
        ch_cram_branch.multiple,
        fasta,
        fai,
        [[id: 'no_gzi'],[]]
    )
    versions = versions.mix(SAMTOOLS_MERGE.out.versions)

    // Mix intervals and no_intervals channels together
    ch_recalibrated_cram = SAMTOOLS_MERGE.out.cram.mix(ch_cram_branch.single)

    // Index CRAM
    SAMTOOLS_INDEX(ch_recalibrated_cram)
    versions = versions.mix(SAMTOOLS_INDEX.out.versions)

    emit:
    recalibrated_cram = ch_recalibrated_cram
    recalibrated_crai = SAMTOOLS_INDEX.out.crai
    multiqc_files
    versions
}
