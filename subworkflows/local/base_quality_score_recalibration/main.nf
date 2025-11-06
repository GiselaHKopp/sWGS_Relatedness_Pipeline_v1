/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//include { GATK4_ANALYZECOVARIATES } from '../../../modules/nf-core/gatk4/analyzecovariates'
include { GATK4_APPLYBQSR         } from '../../../modules/nf-core/gatk4/applybqsr'
include { GATK4_BASERECALIBRATOR  } from '../../../modules/nf-core/gatk4/baserecalibrator'
include { GATK4_GATHERBQSRREPORTS } from '../../../modules/nf-core/gatk4/gatherbqsrreports/main'
include { SAMTOOLS_INDEX          } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE          } from '../../../modules/nf-core/samtools/merge/main'

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

    // Combine CRAM with intervals (scatter)
    ch_cram = cram.join(crai)
    .combine(intervals)
    .map { cram_meta, cram_file, crai_file, interval_meta, interval_file, num_intervals ->
        // Construct new ID: sampleID_intervalName
        def new_id = "${cram_meta.id}_${interval_meta.interval_name}"

        // Merge metadata and overwrite id
        def meta = cram_meta + interval_meta + [ id: new_id ] + [ num_intervals: num_intervals ]

        tuple(meta, cram_file, crai_file, interval_file)
    }

    // Run BaseRecalibrator
    GATK4_BASERECALIBRATOR(
        ch_cram,
        fasta,
        fai,
        dict,
        vcf.map{ _meta, files -> [['id' : 'known_sites'], files]},
        tbi.map{ _meta, files -> [['id' : 'known_sites'], files]}
    )
    versions = versions.mix(GATK4_BASERECALIBRATOR.out.versions)

    // Figuring out if there is one or more table(s) from the same sample
    ch_table_to_merge = GATK4_BASERECALIBRATOR.out.table
        .map{ meta, table ->
            // Use sample name as key, ensure num_intervals is available
            def key = meta.RGSM ?: meta.sample_id ?: meta.id.split('_')[0]

            // Remove interval_name from meta in order to group by sample only
            tuple(meta - meta.subMap('interval_name') + [ id: key ], table)
        }
        .groupTuple()
        .branch{ tuple ->
        // Use meta.num_intervals to asses number of intervals
        single:   tuple[0].num_intervals <= 1
        multiple: tuple[0].num_intervals > 1
    }

    // Only when using intervals
    GATK4_GATHERBQSRREPORTS(ch_table_to_merge.multiple)
    versions = versions.mix(GATK4_GATHERBQSRREPORTS.out.versions)

    // Mix intervals and no_intervals channels together
    ch_table_bqsr = GATK4_GATHERBQSRREPORTS.out.table.mix(ch_table_to_merge.single.map{ meta, table -> [ meta, table[0] ] })
        // Remove no longer necessary field: num_intervals
        .map{ meta, table -> [ meta - meta.subMap('num_intervals'), table ] }.dump(tag: 'Final BQSR Tables')

    // Combine CRAM with BQSR table
    ch_cram_with_table = ch_cram.combine(ch_table_bqsr).dump(tag: 'ch_cram.combine(ch_table_bqsr)')
        .map{ meta, cram_file, crai_file, interval_file, _meta_table, bqsr_table ->
            tuple( meta, cram_file, crai_file, bqsr_table, interval_file )
        }.dump(tag: 'BQSR Input CRAM with BQSR Table')

    // Run ApplyBQSR
    GATK4_APPLYBQSR(
        ch_cram_with_table,
        fasta.map { _meta, fasta_file -> [fasta_file] },
        fai.map { _meta, fai_file -> [fai_file] },
        dict.map { _meta, dict_file -> [dict_file] },
    )
    versions = versions.mix(GATK4_APPLYBQSR.out.versions)

    // Merge recalibrated CRAMs if needed
    ch_cram_branch = GATK4_APPLYBQSR.out.cram.map{ meta, table -> [ groupKey(meta, meta.num_intervals), table ] }.groupTuple()
        .dump(tag: 'Recalibrated CRAMs, joined with BAIs')
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
