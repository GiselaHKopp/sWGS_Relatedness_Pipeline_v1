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
    cram    // channel: [ meta, cram, crai, intervals ]
    vcf     // channel: [ meta, vcf]
    tbi     // channel: [ meta, tbi]

    main:
    versions = channel.empty()
    cram.dump(tag: 'CRAM_BASERECALIBRATOR (cram)')

    // Run BaseRecalibrator
    GATK4_BASERECALIBRATOR(
        cram,
        fasta,
        fai,
        dict,
        vcf.map{ _meta, files -> [['id' : 'known_sites'], files]},
        tbi.map{ _meta, files -> [['id' : 'known_sites'], files]}
    )
    versions = versions.mix(GATK4_BASERECALIBRATOR.out.versions)
    GATK4_BASERECALIBRATOR.out.table.dump(tag: 'CRAM_BASERECALIBRATOR (out.table)')
    // Figuring out if there is one or more table(s) from the same sample
    ch_table_to_merge = GATK4_BASERECALIBRATOR.out.table
        .map{ meta, table ->
            // Use sample name as key, ensure num_intervals is available
            def key = meta.RGSM ?: meta.id.split('_')[0]
            def new_meta = meta - meta.subMap('interval_name')
            new_meta.id = "${key}_${meta.bootstrapping_round}"
            // Remove interval_name from meta in order to group by sample only
            tuple(new_meta, table)
        }.dump(tag: 'CRAM_BASERECALIBRATOR (GATK4_BASERECALIBRATOR.out.table.map())')
        .groupTuple().dump(tag: 'CRAM_BASERECALIBRATOR (GATK4_BASERECALIBRATOR.out.table.map().groupTuple())')
        .branch{ tuple ->
            // Use meta.num_intervals to asses number of intervals
            single:   tuple[0].num_intervals <= 1
            multiple: tuple[0].num_intervals > 1
        }

    // Only when using intervals
    GATK4_GATHERBQSRREPORTS(ch_table_to_merge.multiple)
    versions = versions.mix(GATK4_GATHERBQSRREPORTS.out.versions)

    // Mix intervals and no_intervals channels together
    table_bqsr = GATK4_GATHERBQSRREPORTS.out.table.mix(ch_table_to_merge.single.map{ meta, table -> [ meta, table[0] ] })
        // Remove no longer necessary field: num_intervals
        .map{ meta, table -> [ meta - meta.subMap('num_intervals'), table ] }

    emit:
    table_bqsr
    versions
}
