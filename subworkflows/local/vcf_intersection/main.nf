/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_ISEC } from '../../../modules/nf-core/bcftools/isec'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow VCF_INTERSECTION {
    take:
    vcf_tool1 // channel: [ meta1, vcf]
    tbi_tool1 // channel: [ meta1, tbi]
    vcf_tool2 // channel: [ meta2, vcf]
    tbi_tool2 // channel: [ meta2, tbi]

    main:
    versions = channel.empty()

    // Prepare input for BCFTOOLS_ISEC
    ch_isec_input = vcf_tool1.join(tbi_tool1)
        .join(vcf_tool2.join(tbi_tool2))
        .map { meta1, vcf1, tbi1, meta2, vcf2, tbi2 ->
            def meta = meta1 + [ id: "${meta1.id}_vs_${meta2.id}" ]
            tuple(meta, [vcf1, vcf2])
    }
    .dump(tag: 'BAM_INTERSECTION - isec input')
    // Run BCFTOOLS_ISEC
    BCFTOOLS_ISEC(ch_isec_input)
    versions = versions.mix(BCFTOOLS_ISEC.out.versions)

    // Collect intersection output
    intersection = BCFTOOLS_ISEC.out.results
        .map { meta, dir ->
            def file_common = file("${dir}/0002.vcf.gz")
            def newname = "${meta.id}.isec.common.vcf.gz"
            tuple(meta, file_common.renameTo(newname))
        }.dump(tag: 'BAM_INTERSECTION - intersection output')

    emit:
    intersection
    versions
}
