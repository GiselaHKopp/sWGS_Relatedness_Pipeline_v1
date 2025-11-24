/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_MPILEUP } from '../../../modules/nf-core/bcftools/mpileup/main'
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

    // Combine BAM and intervals
    COMBINE_CRAM_INTERVALS(intervals, SAMTOOLS_CONVERT.out.bam)
    ch_bam_intervals = COMBINE_CRAM_INTERVALS.out.cram_intervals.dump(tag: 'CALL_VARIANTS_BCFTOOLS (COMBINE_CRAM_INTERVALS.out.cram_intervals)')
        .map { meta, bam_file, interval_file ->
            def sample_name = meta.RGSM ?: meta.id.split('_')[0]
            def new_meta = meta + [ id: sample_name ]
            tuple(new_meta, bam_file, interval_file)
        }.dump(tag: 'CALL_VARIANTS_BCFTOOLS (ch_bam_intervals)')

    // Run Bcftools mpileup
    keep_bcftools_mpileup = false
    BCFTOOLS_MPILEUP(ch_bam_intervals, fasta, keep_bcftools_mpileup)
    versions = versions.mix(BCFTOOLS_MPILEUP.out.versions)
    multiqc_files = multiqc_files.mix(BCFTOOLS_MPILEUP.out.stats.map { tuple -> tuple[1] })

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_mpileup = BCFTOOLS_MPILEUP.out.vcf.dump(tag: 'CALL_VARIANTS_BCFTOOLS (BCFTOOLS_MPILEUP.out.vcf)')
        .branch { tuple ->
            single: tuple[0].num_intervals <= 1
            multiple: tuple[0].num_intervals > 1
        }
    tbi_mpileup = BCFTOOLS_MPILEUP.out.tbi.dump(tag: 'CALL_VARIANTS_BCFTOOLS (BCFTOOLS_MPILEUP.out.tbi)')
        .branch { tuple ->
            single: tuple[0].num_intervals <= 1
            multiple: tuple[0].num_intervals > 1
        }

    // Merge VCF
    vcf_to_merge = vcf_mpileup.multiple.dump(tag: 'CALL_VARIANTS_BCFTOOLS (vcf_mpileup.multiple)')
        .map { meta, vcf ->
            def new_id = meta.id + "_merged" + ".${meta.variantcaller}"
            def new_meta = meta + [ id: new_id ] - meta.subMap('interval_name', 'num_intervals')
            tuple(new_meta, vcf)
        }.dump(tag: 'CALL_VARIANTS_BCFTOOLS (vcf_mpileup.multiple.map())')
        .groupTuple().dump(tag: 'CALL_VARIANTS_BCFTOOLS (vcf_to_merge)')
    GATK4_MERGEVCFS(vcf_to_merge, dict)

    // Mix single and multiple channels together
    vcf = GATK4_MERGEVCFS.out.vcf.dump(tag: 'CALL_VARIANTS_BCFTOOLS (GATK4_MERGEVCFS.out.vcf)')
        .mix(vcf_mpileup.single)
        .dump(tag: 'CALL_VARIANTS_BCFTOOLS (vcf)')
    tbi = GATK4_MERGEVCFS.out.tbi.dump(tag: 'CALL_VARIANTS_BCFTOOLS (GATK4_MERGEVCFS.out.tbi)')
        .mix(tbi_mpileup.single)
        .dump(tag: 'CALL_VARIANTS_BCFTOOLS (tbi)')

    emit:
    vcf
    tbi
    multiqc_files
    versions
}
