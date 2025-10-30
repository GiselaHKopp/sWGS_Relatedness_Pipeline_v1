/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_SORT          } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_STATS         } from '../../../modules/nf-core/bcftools/stats'
include { GATK4_GENOMICSDBIMPORT } from '../../../modules/nf-core/gatk4/genomicsdbimport'
include { GATK4_GENOTYPEGVCFS    } from '../../../modules/nf-core/gatk4/genotypegvcfs'
include { GATK4_HAPLOTYPECALLER  } from '../../../modules/nf-core/gatk4/haplotypecaller'
include { GATK4_MERGEVCFS        } from '../../../modules/nf-core/gatk4/mergevcfs'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow CALL_VARIANTS {
    take:
    fasta       // tuple(meta2, path_to_fasta)                          e.g. [ id: 'ref' ], ref.fasta
    fai         // tuple(meta, path_to_fasta.fai)                       e.g. [ id: 'ref' ], ref.fasta.fai
    dict        // tuple(meta, path_to_dict)                            e.g. [ id: 'ref' ], ref.dict
    intervals   // tuple(meta, path_to_intervals, number_of_intervals)  e.g. [[ id: 'ref', interval_name:'scaffold'], intervals.bed, number_of_intervals]
    bam         // tuple(meta, path_to_bam)                             e.g. [ id: 'sample1' ], sample1.bam
    bai         // tuple(meta, path_to_bai)                             e.g. [ id: 'sample1' ], sample1.bam.bai

    main:
    versions = channel.empty()
    multiqc_files = channel.empty()

    // Combine BAM with intervals (scatter)
    ch_inputs = bam.join(bai)
    .combine(intervals)
    .map { bam_meta, bamfile, baifile, interval_meta, interval_file, _num_intervals ->
        // Construct new ID: sampleID_intervalName
        def new_id = "${bam_meta.id}_${interval_meta.interval_name}"

        // Merge metadata and overwrite id
        def meta = bam_meta + interval_meta + [ id: new_id ]

        tuple(meta, bamfile, baifile, interval_file, [])
    }

    // Run GATK HaplotypeCaller
    GATK4_HAPLOTYPECALLER(ch_inputs, fasta, fai, dict, [[id: 'no_dbsnp'], []], [[id: 'no_dbsnp_tbi'], []])
    versions = versions.mix(GATK4_HAPLOTYPECALLER.out.versions)

    // Prepare for GenomicsDBImport (scatter)
    ch_gvcfs = GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi)
        .map { meta, vcf, tbi -> tuple(meta.interval_name, vcf, tbi) }

    // Key intervals by interval_name
    ch_intervals_keyed = intervals.map { meta, bed, num_intervals ->
        tuple(meta.interval_name, bed, num_intervals)
    }

    // Prepare GenomicsDBImport input by grouping GVCFs by interval_name
    ch_gdb_input = ch_gvcfs
        .groupTuple()
        .join(ch_intervals_keyed)
        .map { interval_name, vcfs, tbis, bed, _num_intervals ->
            def meta = [id: "joint_${interval_name}", interval_name: interval_name]
            tuple(
                meta,
                vcfs,
                tbis,
                bed,
                [],
                file('.')
            )
        }

    // Run GATK GenomicsDBImport
    GATK4_GENOMICSDBIMPORT(ch_gdb_input, false, false, false)
    versions = versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)

    // Run GATK GenotypeGVCFs
    ch_gtp_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map { meta, genomicsdb -> tuple(meta, genomicsdb, [], [], []) }
    GATK4_GENOTYPEGVCFS(ch_gtp_input, fasta, fai, dict, [[id: 'no_dbsnp'], []], [[id: 'no_dbsnp_tbi'], []])
    versions = versions.mix(GATK4_GENOTYPEGVCFS.out.versions)

    // Run BCFtools stats
    GATK4_GENOTYPEGVCFS.out.vcf.join(GATK4_GENOTYPEGVCFS.out.tbi)
    .map { meta, vcf, tbi -> tuple(meta, vcf, tbi) }
    .set { ch_vcf_tbi }
    BCFTOOLS_STATS(ch_vcf_tbi, [[id: 'no_regions'], []], [[id: 'no_targets'], []], [[id: 'no_samples'], []], [[id: 'no_exons'], []], fasta)
    multiqc_files = multiqc_files.mix(BCFTOOLS_STATS.out.stats.map { tuple -> tuple[1] })
    versions = versions.mix(BCFTOOLS_STATS.out.versions)

    // Sort each scaffold VCF before merging
    GATK4_GENOTYPEGVCFS.out.vcf
        .map { meta, vcf ->
            def new_meta = meta + [ id: "${meta.id}.sorted" ]
            tuple(new_meta, vcf) }
        .set { ch_vcfs }

    BCFTOOLS_SORT(ch_vcfs)
    versions = versions.mix(BCFTOOLS_SORT.out.versions)

    // Collect sorted VCFs into one tuple for merging
    BCFTOOLS_SORT.out.vcf
        .map { _meta, vcf -> vcf }
        .collect()
        .map { vcfs -> tuple([id: "joint_merged"], vcfs) }
        .set { ch_merge_vcfs }

    // Merge all intervals into one VCF
    GATK4_MERGEVCFS(ch_merge_vcfs, dict)
    versions = versions.mix(GATK4_MERGEVCFS.out.versions)

    emit:
    vcf = GATK4_MERGEVCFS.out.vcf
    tbi = GATK4_MERGEVCFS.out.tbi
    multiqc_files
    versions
}
