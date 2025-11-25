/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_MPILEUP } from '../../../modules/local/bcftools/mpileup/main'
include { GATK4_MERGEVCFS  } from '../../../modules/nf-core/gatk4/mergevcfs'
include { SAMTOOLS_CONVERT } from '../../../modules/nf-core/samtools/convert/main'

include { COMBINE_CRAM_INTERVALS } from '../combine_cram_intervals'

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

    // Add variant callerrapping metadata to channels
    fasta.map { meta, fasta_file ->
        tuple( meta + [variantcaller: 'bcftools'], fasta_file ) }
        .set { fasta }
    fai.map { meta, fai_file ->
        tuple( meta + [variantcaller: 'bcftools'], fai_file ) }
        .set { fai }
    dict.map { meta, dict_file ->
        tuple( meta + [variantcaller: 'bcftools'], dict_file ) }
        .set { dict }
    intervals.map { meta, interval_file, num_intervals ->
        tuple( meta + [variantcaller: 'bcftools'], interval_file, num_intervals ) }
        .set { intervals }
    cram.map { meta, cram_file ->
        tuple( meta + [variantcaller: 'bcftools'], cram_file ) }
        .set { cram }
    crai.map { meta, crai_file ->
        tuple( meta + [variantcaller: 'bcftools'], crai_file ) }
        .set { crai }

    // Convert CRAM/CAI to BAM/BAI
    ch_cram_crai_to_convert = cram.join(crai)
    SAMTOOLS_CONVERT(ch_cram_crai_to_convert, fasta, fai)
    versions = versions.mix(SAMTOOLS_CONVERT.out.versions)

    // Collect a list of all BAM files
    ch_bams = SAMTOOLS_CONVERT.out.bam.dump(tag: 'CALL_VARIANTS_BCFTOOLS (SAMTOOLS_CONVERT.out.bam)')
        .map { _meta, bam -> file(bam) }.dump(tag: 'CALL_VARIANTS_BCFTOOLS (SAMTOOLS_CONVERT.out.map())')
        .collect().dump(tag: 'CALL_VARIANTS_BCFTOOLS (SAMTOOLS_CONVERT.out.map().collect())')

    ch_intervals = intervals
        .map { meta, interval_file, _num_intervals ->
            tuple(meta, interval_file)
        }.dump(tag: 'CALL_VARIANTS_BCFTOOLS (ch_intervals)')

    // Run Bcftools mpileup
    keep_bcftools_mpileup = false
    BCFTOOLS_MPILEUP(ch_intervals, ch_bams, fasta, keep_bcftools_mpileup)
    versions = versions.mix(BCFTOOLS_MPILEUP.out.versions)

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_mpileup = BCFTOOLS_MPILEUP.out.vcf.dump(tag: 'CALL_VARIANTS_BCFTOOLS (BCFTOOLS_MPILEUP.out.vcf)')
        .branch { tuple ->
            single: tuple[0].num_intervals <= 1
            multiple: tuple[0].num_intervals > 1
        }

    emit:
    //vcf = vcf_mpileup.multiple
    vcf = BCFTOOLS_MPILEUP.out.vcf
    multiqc_files
    versions
}
