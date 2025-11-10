/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_MPILEUP } from '../../../modules/nf-core/bcftools/mpileup/main'
include { GATK4_MERGEVCFS  } from '../../../modules/nf-core/gatk4/mergevcfs'
include { SAMTOOLS_CONVERT } from '../../../modules/nf-core/samtools/convert/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow CALL_VARIANTS_BCFTOOLS {
    take:
    fasta       // channel: [ meta, fasta]
    fai         // channel: [ meta, fai]
    dict        // channel: [ meta, dict]
    intervals   // channel: [ meta, intervals, number_of_intervals]
    cram        // channel: [ meta, cram]
    crai        // channel: [ meta, crai]

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
    keep_bcftools_mpileup = false
    BCFTOOLS_MPILEUP(ch_bam_bai_intervals, fasta, keep_bcftools_mpileup)
    versions = versions.mix(BCFTOOLS_MPILEUP.out.versions)
    multiqc_files = multiqc_files.mix(BCFTOOLS_MPILEUP.out.stats.map { tuple -> tuple[1] })

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_mpileup = BCFTOOLS_MPILEUP.out.vcf.branch { tuple ->
        single: tuple[0].num_intervals > 1
        multiple: tuple[0].num_intervals <= 1
    }

    // Merge VCF
    vcf_to_merge = vcf_mpileup.single.map { meta, vcf -> [groupKey(meta, meta.num_intervals), vcf] }.groupTuple()
    GATK4_MERGEVCFS(vcf_to_merge, dict)

    // Mix single and multiple channels together
    vcf = GATK4_MERGEVCFS.out.vcf
        .mix(vcf_mpileup.multiple)
        .map { meta, vcf -> [meta - meta.subMap('num_intervals') + [variantcaller: 'bcftools'], vcf] }

    // Merge TBI
    tbi = GATK4_MERGEVCFS.out.tbi
        .mix(BCFTOOLS_MPILEUP.out.tbi.filter { meta, _tbi -> meta.num_intervals <= 1 })
        .map { meta, tbi -> [meta - meta.subMap('num_intervals') + [variantcaller: 'bcftools'], tbi] }

    emit:
    vcf
    tbi
    multiqc_files
    versions
}
