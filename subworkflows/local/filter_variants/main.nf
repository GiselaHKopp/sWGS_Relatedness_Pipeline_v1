include { GATK4_MERGEVCFS         } from '../../../modules/nf-core/gatk4/mergevcfs'
include { GATK4_SELECTVARIANTS    } from '../../../modules/nf-core/gatk4/selectvariants'
include { GATK4_VARIANTFILTRATION } from '../../../modules/nf-core/gatk4/variantfiltration'

workflow FILTER_VARIANTS {
    take:
    vcf         // tuple(meta, path_to_vcf)         e.g. [ id: 'sample1' ], sample1.vcf
    tbi         // tuple(meta, path_to_tbi)         e.g. [ id: 'sample1' ], sample1.vcf.tbi
    fai         // tuple(meta, path_to_fasta.fai)   e.g. [ id: 'ref' ], ref.fasta.fai
    dict        // tuple(meta, path_to_dict)        e.g. [ id: 'ref' ], ref.dict
    ref_fasta   // tuple(meta2, path_to_fasta)      e.g. [ id: 'ref' ], ref.fasta

    main:
    // Collect software versions and QC reports
    versions = channel.empty()
    multiqc_files = channel.empty()

    // Filter variants to exclude low-quality calls
    vcf.join(tbi).map { meta, vcf_unfiltered, tbi_unfiltered ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}.filtered"
            tuple(new_meta, vcf_unfiltered, tbi_unfiltered)
            }.set { ch_filtered_input }

    GATK4_VARIANTFILTRATION(ch_filtered_input, ref_fasta, fai, dict, [[id: 'no_gzi'], []])
    versions = versions.mix(GATK4_VARIANTFILTRATION.out.versions)

    // Select only passing variants
    GATK4_VARIANTFILTRATION.out.vcf.join(GATK4_VARIANTFILTRATION.out.tbi)
    .map { meta, vcf_filtered, tbi_filtered ->
        def new_meta = meta.clone()
        new_meta.id = "${meta.id}.selected"
        tuple(new_meta, vcf_filtered, tbi_filtered, [])
        }
        .set { ch_selected_input }
    GATK4_SELECTVARIANTS(ch_selected_input)
    versions = versions.mix(GATK4_SELECTVARIANTS.out.versions)

    // Gather VCFs into single file
    GATK4_SELECTVARIANTS.out.vcf
        .map { _meta, vcf_selected -> vcf_selected } // discard meta
        .collect()
        .map { vcfs_selected ->
            def meta = [ id: "merged" ]
            tuple(meta, vcfs_selected)
        }
        .set { ch_merge_vcf_input }
    GATK4_MERGEVCFS(ch_merge_vcf_input, dict)
    versions = versions.mix(GATK4_MERGEVCFS.out.versions)

    emit:
    vcf = GATK4_MERGEVCFS.out.vcf
    tbi = GATK4_MERGEVCFS.out.tbi
    multiqc_files
    versions
}
