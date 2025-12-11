/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GATK4_SELECTVARIANTS    } from '../../../modules/nf-core/gatk4/selectvariants'
include { GATK4_VARIANTFILTRATION } from '../../../modules/nf-core/gatk4/variantfiltration'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow FILTER_VARIANTS {
    take:
    fasta   // channel: [ meta, fasta]
    fai     // channel: [ meta, fai
    dict    // channel: [ meta, dict]
    vcf     // channel: [ meta, vcf]
    tbi     // channel: [ meta, tbi]

    main:
    // Collect software versions and QC reports
    versions = channel.empty()
    multiqc_files = channel.empty()

    vcf_branches = vcf.branch {
        filter    : !params.skip_filter_variants
        no_filter :  params.skip_filter_variants
    }

    tbi_branches = tbi.branch {
        filter    : !params.skip_filter_variants
        no_filter :  params.skip_filter_variants
    }

    // Filter variants to exclude low-quality calls
    ch_filtered_input = vcf_branches.filter.join(tbi_branches.filter)
        .map { meta, vcf_unfiltered, tbi_unfiltered ->
            def new_meta = meta + [ id: "${meta.id}.filtered" ]
            tuple(new_meta, vcf_unfiltered, tbi_unfiltered)
        }

    GATK4_VARIANTFILTRATION(ch_filtered_input, fasta, fai, dict, [[id: 'no_gzi'], []])
    versions = versions.mix(GATK4_VARIANTFILTRATION.out.versions)

    // Select only passing variants
    ch_selected_input = GATK4_VARIANTFILTRATION.out.vcf.join(GATK4_VARIANTFILTRATION.out.tbi)
        .map { meta, vcf_filtered, tbi_filtered ->
            def new_meta = meta + [ id: "${meta.id}.selected" ]
            tuple(new_meta, vcf_filtered, tbi_filtered, [])
        }
    GATK4_SELECTVARIANTS(ch_selected_input)
    versions = versions.mix(GATK4_SELECTVARIANTS.out.versions)

    vcf_output = vcf_branches.no_filter.mix(GATK4_SELECTVARIANTS.out.vcf)
    tbi_output = vcf_branches.no_filter.mix(GATK4_SELECTVARIANTS.out.tbi)

    emit:
    vcf = vcf_output
    tbi = tbi_output
    multiqc_files
    versions
}
