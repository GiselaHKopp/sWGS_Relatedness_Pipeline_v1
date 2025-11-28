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
    tbi_tool2 // channel: [ meta1, tbi]

    main:
    versions = channel.empty()

    vcf_tool1_prepared = vcf_tool1
        .map { meta, vcf ->
            tuple(meta + [id: 'called_variants'] - meta.subMap('variantcaller'), vcf)
        }
    vcf_tool2_prepared = vcf_tool2
        .map { meta, vcf ->
            tuple(meta + [id: 'called_variants'] - meta.subMap('variantcaller'), vcf)
        }
    tbi_tool1_prepared = tbi_tool1
        .map { meta, tbi ->
            tuple(meta + [id: 'called_variants'] - meta.subMap('variantcaller'), tbi)
        }
    tbi_tool2_prepared = tbi_tool2
        .map { meta, tbi ->
            tuple(meta + [id: 'called_variants'] - meta.subMap('variantcaller'), tbi)
        }

    ch_vcfs = vcf_tool1_prepared
        .mix(vcf_tool2_prepared)
        .groupTuple()
        .map { meta, vcfs -> tuple(meta, vcfs) }

    ch_tbis = tbi_tool1_prepared
        .mix(tbi_tool2_prepared)
        .groupTuple()
        .map { meta, tbis -> tuple(meta, tbis) }

    ch_isec_input = ch_vcfs
        .join(ch_tbis)
        .map { meta, vcfs, tbis -> tuple(meta, vcfs, tbis) }

    // Run BCFTOOLS_ISEC
    BCFTOOLS_ISEC(ch_isec_input)
    versions = versions.mix(BCFTOOLS_ISEC.out.versions)

    // Collect intersection output
    intersection = BCFTOOLS_ISEC.out.results
        .map { meta, dir ->
            def file_common = file("${dir}/0002.vcf")
            tuple(meta, file_common)
        }

    emit:
    intersection
    versions
}
