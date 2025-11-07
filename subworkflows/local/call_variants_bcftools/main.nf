/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_CALL    } from '../../../modules/nf-core/bcftools/call/main'
include { BCFTOOLS_CONCAT  } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_INDEX   } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_MPILEUP } from '../../../modules/nf-core/bcftools/mpileup/main'
include { SAMTOOLS_CONVERT } from '../../../modules/nf-core/samtools/convert/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow CALL_VARIANTS_BCFTOOLS {
    take:
    fasta       // tuple(meta2, path_to_fasta)                          e.g. [ id: 'ref' ], ref.fasta
    fai         // tuple(meta, path_to_fasta.fai)                       e.g. [ id: 'ref' ], ref.fasta.fai
    intervals   // tuple(meta, path_to_intervals, number_of_intervals)  e.g. [[ id: 'ref', interval_name:'scaffold'], intervals.bed, number_of_intervals]
    cram        // tuple(meta, path_to_cram)                            e.g. [ id: 'sample1' ], sample1.cram
    crai        // tuple(meta, path_to_crai)                            e.g. [ id: 'sample1' ], sample1.cram.crai

    main:
    versions = channel.empty()
    multiqc_files = channel.empty()

    // Convert CRAM/CAI to BAM/BAI
    ch_cram_crai_to_convert = cram.join(crai)
    SAMTOOLS_CONVERT(ch_cram_crai_to_convert, fasta, fai)
    versions = versions.mix(SAMTOOLS_CONVERT.out.versions)

    // Combine BAM with BAI and intervals
    ch_bam_bai_intervals = SAMTOOLS_CONVERT.out.bam
    .combine(intervals)
    .map { bam_meta, bam_file, interval_meta, interval_file, _num_intervals ->
        // Construct new ID: sampleID_intervalName
        def new_id = "${bam_meta.id}_${interval_meta.interval_name}"

        // Merge metadata and overwrite id
        def meta = bam_meta + interval_meta + [ id: new_id ]

        tuple(meta, bam_file, interval_file)
    }

    // Run Bcftools mpileup
    BCFTOOLS_MPILEUP(ch_bam_bai_intervals, fasta, false)
    versions = versions.mix(BCFTOOLS_MPILEUP.out.versions)
    multiqc_files = multiqc_files.mix(BCFTOOLS_MPILEUP.out.stats.map { tuple -> tuple[1] })

    // Run Bcftools call
    ch_vcf_tbi_call = BCFTOOLS_MPILEUP.out.vcf.join(BCFTOOLS_MPILEUP.out.tbi).dump(tag: 'ch_vcf_tbi_call')
    BCFTOOLS_CALL(ch_vcf_tbi_call, channel.empty(), channel.empty(), channel.empty())
    versions = versions.mix(BCFTOOLS_CALL.out.versions)

    // Run Bcftools concat
    ch_vcf_tbi_concat = BCFTOOLS_CALL.out.vcf.join(BCFTOOLS_CALL.out.tbi)
        .collect().dump(tag: 'ch_vcf_tbi_concat')
    BCFTOOLS_CONCAT(ch_vcf_tbi_concat)
    versions = versions.mix(BCFTOOLS_CONCAT.out.versions)

    // Run Bcftools index
    BCFTOOLS_INDEX(BCFTOOLS_CONCAT.out.vcf)
    versions = versions.mix(BCFTOOLS_INDEX.out.versions)

    emit:
    vcf = BCFTOOLS_CONCAT.out.vcf
    tbi = BCFTOOLS_INDEX.out.tbi
    multiqc_files
    versions
}
