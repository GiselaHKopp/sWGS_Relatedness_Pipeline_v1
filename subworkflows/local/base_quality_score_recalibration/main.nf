/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GATK4_ANALYZECOVARIATES                    } from '../../../modules/local/gatk4/analyzecovariates'
include { GATK4_APPLYBQSR                            } from '../../../modules/nf-core/gatk4/applybqsr'
include { SAMTOOLS_INDEX                             } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_SCATTERED } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                             } from '../../../modules/nf-core/samtools/merge/main'

include { COMBINE_CRAM_CRAI_INTERVALS                                } from '../combine_cram_crai_intervals'
include { CRAM_BASERECALIBRATOR                                      } from '../cram_baserecalibrator'
include { CRAM_BASERECALIBRATOR as CRAM_BASERECALIBRATOR_SECOND_PASS } from '../cram_baserecalibrator'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow BASE_QUALITY_SCORE_RECALIBRATION {
    take:
    fasta       // channel: [ meta, fasta]
    fai         // channel: [ meta, fai]
    dict        // channel: [ meta, dict]
    intervals   // channel: [ meta, intervals, number_of_intervals]
    cram        // channel: [ meta, cram]
    crai        // channel: [ meta, crai]
    vcf         // channel: [ meta, vcf]
    tbi         // channel: [ meta, tbi]

    main:
    versions = channel.empty()
    multiqc_files = channel.empty()

    // Combine CRAM with intervals
    COMBINE_CRAM_CRAI_INTERVALS(intervals, cram, crai)
    COMBINE_CRAM_CRAI_INTERVALS.out.cram_crai_intervals
        .map { meta, cram_file, crai_file, interval_file ->
            def new_id = meta.id + (meta.bootstrapping_round ? "_${meta.bootstrapping_round}" : "")
            def new_meta = meta + [ id: new_id ]
            tuple(new_meta, cram_file, crai_file, interval_file)
        }.dump(tag: 'BQSR (combined_cram_crai_intervals)')
        .set { combined_cram_crai_intervals }

    // Run BaseRecalibrator
    CRAM_BASERECALIBRATOR(fasta, fai, dict, combined_cram_crai_intervals, vcf, tbi)
    versions = versions.mix(CRAM_BASERECALIBRATOR.out.versions)

    CRAM_BASERECALIBRATOR.out.table_bqsr.dump(tag: 'BQSR (CRAM_BASERECALIBRATOR.out.table_bqsr)')

    // Combine CRAM with BQSR table
    ch_cram_with_table = combined_cram_crai_intervals
        .combine(CRAM_BASERECALIBRATOR.out.table_bqsr).dump(tag: 'BQSR (combined_cram_crai_intervals.combine())')
        .filter { meta_cc, _cram_file, _crai_file, _interval_file, meta_tab, _table ->
            // only keep pairs where sample IDs match
            meta_cc.RGSM == meta_tab.RGSM
        }.dump(tag: 'BQSR (combined_cram_crai_intervals.combine().filter())')
        .map { meta_cram, cram_file, crai_file, interval_file, _meta_table, table ->
            tuple(meta_cram, cram_file, crai_file, table, interval_file)
        }.dump(tag: 'BQSR (ch_cram_with_table)')

    // Run ApplyBQSR
    GATK4_APPLYBQSR(
        ch_cram_with_table,
        fasta.map { _meta, fasta_file -> [fasta_file] },
        fai.map { _meta, fai_file -> [fai_file] },
        dict.map { _meta, dict_file -> [dict_file] },
    )
    versions = versions.mix(GATK4_APPLYBQSR.out.versions)

/*
    // TODO: Do we need quality control per-sample or per-sample-per-interval
    // Index recalibrated CRAMs (per interval)
    SAMTOOLS_INDEX_SCATTERED(GATK4_APPLYBQSR.out.cram)
    versions = versions.mix(SAMTOOLS_INDEX_SCATTERED.out.versions)

    SAMTOOLS_INDEX_SCATTERED.out.crai.dump(tag: 'BQSR (SAMTOOLS_INDEX_SCATTERED.out.crai)')

    // Combine recalibrated CRAM with intervals for second pass
    ch_cram_second_pass = GATK4_APPLYBQSR.out.cram.join(SAMTOOLS_INDEX_SCATTERED.out.crai).dump(tag: 'BQSR (GATK4_APPLYBQSR.out.cram.join())')
        .combine(intervals).dump(tag: 'BQSR (GATK4_APPLYBQSR.out.cram.join().combine())')
        .map { cram_meta, cram_file, crai_file, _interval_meta, interval_file, _num_intervals ->
            // Change id to get distinct filenames
            def new_meta = cram_meta + [ id: "${cram_meta.id}.after" ]
            tuple(new_meta, cram_file, crai_file, interval_file)
        }.dump(tag: 'BQSR (ch_cram_second_pass)')

    // Run BaseRecalibrator (second pass, for quality control)
    CRAM_BASERECALIBRATOR_SECOND_PASS(
        fasta,
        fai,
        dict,
        ch_cram_second_pass,
        vcf.map{ _meta, files -> [['id' : 'known_sites'], files]},
        tbi.map{ _meta, files -> [['id' : 'known_sites'], files]}
    )

    ch_bqsr_first = CRAM_BASERECALIBRATOR.out.table_bqsr.dump(tag: 'BQSR (CRAM_BASERECALIBRATOR.out.table_bqsr) ')
        .map { meta, table ->
            tuple(meta.RGSM ?: meta.id, [meta, table])
        }.dump(tag: 'BQSR (ch_bqsr_first)')
    ch_bqsr_second = CRAM_BASERECALIBRATOR_SECOND_PASS.out.table_bqsr.dump(tag: 'BQSR (CRAM_BASERECALIBRATOR_SECOND_PASS.out.table_bqsr) ')
        .map { meta, table ->
            tuple(meta.RGSM ?: meta.id, [meta, table])
        }.dump(tag: 'BQSR (ch_bqsr_second)')
    ch_bqsr_tables = ch_bqsr_first
        .join(ch_bqsr_second).dump(tag: 'BQSR (ch_bqsr_first.join(ch_bqsr_second))')
        .map { _key, first, second ->
            def meta = first[0]
            def table_before = first[1]
            def table_after  = second[1]

            tuple(meta, table_before, table_after)
        }.dump(tag: 'BQSR (ch_bqsr_tables)')

    // Run AnalyzeCovariates
    GATK4_ANALYZECOVARIATES(ch_bqsr_tables)
    versions = versions.mix(GATK4_ANALYZECOVARIATES.out.versions)
*/
    // Merge recalibrated CRAMs if needed
    ch_cram_branch = GATK4_APPLYBQSR.out.cram
        .map{ meta, cram_file ->
            def new_id = (meta.RGSM ?: meta.id.split('_')[0]) + (meta.bootstrapping_round ? "_${meta.bootstrapping_round}" : "")
            def new_meta = meta + [ id: new_id ] - meta.subMap('interval_name')
            tuple(new_meta, cram_file)
        }
        .groupTuple()
        .branch { tuple ->
            single:   tuple[0].num_intervals <= 1
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
        .map{ meta, cram_file ->
            // Use sample name as key, ensure num_intervals is available
            def key = meta.RGSM ?: meta.id.split('_')[0]

            // Remove interval_name from meta in order to group by sample only
            tuple(meta - meta.subMap('interval_name') - meta.subMap('num_intervals') + [ id: key ], cram_file)
        }

    // Index CRAM
    SAMTOOLS_INDEX(ch_recalibrated_cram)
    versions = versions.mix(SAMTOOLS_INDEX.out.versions)

    emit:
    recalibrated_cram = ch_recalibrated_cram
    recalibrated_crai = SAMTOOLS_INDEX.out.crai
    multiqc_files
    versions
}
