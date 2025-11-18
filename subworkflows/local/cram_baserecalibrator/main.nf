/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GATK4_BASERECALIBRATOR  } from '../../../modules/nf-core/gatk4/baserecalibrator'
include { GATK4_GATHERBQSRREPORTS } from '../../../modules/nf-core/gatk4/gatherbqsrreports/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow CRAM_BASERECALIBRATOR {
    take:
    fasta   // channel: [ meta, fasta]
    fai     // channel: [ meta, fai]
    dict    // channel: [ meta, dict]
    ch_cram // channel: [ meta, cram, crai, intervals ]
    vcf     // channel: [ meta, vcf]
    tbi     // channel: [ meta, tbi]

    main:
    versions = channel.empty()

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
            def key = meta.RGSM ?: meta.id.split('_')[0]

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
        .map{ meta, table -> [ meta - meta.subMap('num_intervals'), table ] }

    emit:
    ch_table_bqsr
    versions
}
