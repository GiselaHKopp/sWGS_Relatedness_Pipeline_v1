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

include { COMBINE_CRAM_CRAI_INTERVALS } from '../combine_cram_crai_intervals'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow CALL_VARIANTS_GATK {
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

    // Combine CRAM with CRAI and intervals
    COMBINE_CRAM_CRAI_INTERVALS(intervals, cram, crai)

    // Prepare HaplotypeCaller input
    ch_haplotypecaller_input = COMBINE_CRAM_CRAI_INTERVALS.out.cram_crai_intervals
        .map { meta, cram_file, crai_file, interval_file ->
            def new_meta = meta + [variantcaller: 'gatk']
            tuple(new_meta, cram_file, crai_file, interval_file, [])
        }

    // Run GATK HaplotypeCaller
    GATK4_HAPLOTYPECALLER(ch_haplotypecaller_input, fasta, fai, dict, [[id: 'no_dbsnp'], []], [[id: 'no_dbsnp_tbi'], []])
    versions = versions.mix(GATK4_HAPLOTYPECALLER.out.versions)

    // Prepare for GenomicsDBImport
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
    ch_vcf_tbi = GATK4_GENOTYPEGVCFS.out.vcf.join(GATK4_GENOTYPEGVCFS.out.tbi)
    .map { meta, vcf, tbi -> tuple(meta, vcf, tbi) }.dump(tag: 'CALL_VARIANTS_GATK (ch_vcf_tbi)')
    BCFTOOLS_STATS(ch_vcf_tbi, [[id: 'no_regions'], []], [[id: 'no_targets'], []], [[id: 'no_samples'], []], [[id: 'no_exons'], []], fasta)
    multiqc_files = multiqc_files.mix(BCFTOOLS_STATS.out.stats.map { tuple -> tuple[1] })
    versions = versions.mix(BCFTOOLS_STATS.out.versions)

    // Sort each interval VCF before merging
    ch_vcfs = GATK4_GENOTYPEGVCFS.out.vcf
        .map { meta, vcf ->
            def new_meta = meta + [ id: "${meta.id}.sorted" ] + [ variantcaller: 'gatk' ]
            tuple(new_meta, vcf)
        }

    BCFTOOLS_SORT(ch_vcfs)
    versions = versions.mix(BCFTOOLS_SORT.out.versions)

    // Collect sorted VCFs into one tuple for merging
    ch_merge_vcfs = BCFTOOLS_SORT.out.vcf
        .map { _meta, vcf -> vcf }.dump(tag: 'CALL_VARIANTS_GATK (ch_merge_vcfs)')
        .collect()
        .map { vcfs ->
            def new_meta = [id: "joint_merged", variantcaller: 'gatk']
            tuple(new_meta, vcfs)
        }

    // Merge all intervals into one VCF
    GATK4_MERGEVCFS(ch_merge_vcfs, dict)
    versions = versions.mix(GATK4_MERGEVCFS.out.versions)

    // Extract the bootstrapping round from any CRAM meta
    ch_bootstrap_round = cram.map { meta, _cram_file -> meta.bootstrapping_round ?: null }.first()

    ch_final_vcf = GATK4_MERGEVCFS.out.vcf
        .combine(ch_bootstrap_round)
        .map { meta, vcf, round ->
            if (round) {
                def new_meta = meta + [
                    id: meta.id + "_${round}",
                    bootstrapping_round: round
                ]
                return tuple(new_meta, vcf)
            }
            return tuple(meta, vcf)
        }

    ch_final_tbi = GATK4_MERGEVCFS.out.tbi
        .combine(ch_bootstrap_round)
        .map { meta, tbi, round ->
            if (round) {
                def new_meta = meta + [
                    id: meta.id + "_${round}",
                    bootstrapping_round: round
                ]
                return tuple(new_meta, tbi)
            }
            return tuple(meta, tbi)
        }

    emit:
    vcf = ch_final_vcf
    tbi = ch_final_tbi
    multiqc_files
    versions
}
