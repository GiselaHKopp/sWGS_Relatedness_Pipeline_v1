/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GAWK as BUILD_INTERVALS } from '../../../modules/nf-core/gawk'
include { SPLIT_INTERVALS         } from '../../../modules/local/splitintervals'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PREPARE_INTERVALS {
    take:
    fai // channel: [ meta, fai]

    main:
    versions = channel.empty()

    // Build intervals from FASTA index
    BUILD_INTERVALS(fai, [], false)
    versions = versions.mix(BUILD_INTERVALS.out.versions)

    intervals_combined = BUILD_INTERVALS.out.output
    .map { meta, intervals ->
        def num_intervals = intervals.readLines().size()
        tuple(meta + [ reference_fasta: meta.id ], intervals, num_intervals)
    }

    // Split intervals into separate files
    SPLIT_INTERVALS(BUILD_INTERVALS.out.output)
    versions = versions.mix(SPLIT_INTERVALS.out.versions)

    intervals_split = SPLIT_INTERVALS.out.bed
        .flatMap { meta, beds ->
            def list = beds instanceof List ? beds : [beds]
            def count = list.size()
            list.collect { bed ->
                def contig = bed.baseName
                def fasta_file = meta.id
                def new_meta = meta + [ id: contig ] + [ interval_name: contig ] + [ reference_fasta: fasta_file ]
                tuple(new_meta, bed, count)
            }
        }

    emit:
    intervals_combined  // [[id:'reference_fasta', reference_fasta: 'reference_fasta']], interval.bed, number_of_intervals]
    intervals_split     // [[id:'scaffold', interval_name:'scaffold', reference_fasta: 'reference_fasta'], interval.bed, number_of_intervals]
    versions
}
