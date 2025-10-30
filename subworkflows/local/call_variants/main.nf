include { BCFTOOLS_STATS          } from '../../../modules/nf-core/bcftools/stats'
include { GATK4_GENOMICSDBIMPORT  } from '../../../modules/nf-core/gatk4/genomicsdbimport'
include { GATK4_GENOTYPEGVCFS     } from '../../../modules/nf-core/gatk4/genotypegvcfs'
include { GATK4_HAPLOTYPECALLER   } from '../../../modules/nf-core/gatk4/haplotypecaller'
include { GATK4_MERGEVCFS         } from '../../../modules/nf-core/gatk4/mergevcfs'

workflow CALL_VARIANTS {
    take:
    bam         // tuple(meta, path_to_bam)         e.g. [ id: 'sample1' ], sample1.bam
    bam_index   // tuple(meta, path_to_bai)         e.g. [ id: 'sample1' ], sample1.bam.bai
    fai         // tuple(meta, path_to_fasta.fai)   e.g. [ id: 'ref' ], ref.fasta.fai
    dict        // tuple(meta, path_to_dict)        e.g. [ id: 'ref' ], ref.dict
    ref_fasta   // tuple(meta2, path_to_fasta)      e.g. [ id: 'ref' ], ref.fasta

    main:
    // Collect software versions and QC reports
    versions = channel.empty()
    multiqc_files = channel.empty()

    // Join BAM and BAI channels
    ch_inputs = bam.join(bam_index)
        .map { meta, bamfile, baifile ->
            tuple(meta, bamfile, baifile, [], [])
        }

    // Run GATK HaplotypeCaller
    GATK4_HAPLOTYPECALLER(ch_inputs, ref_fasta, fai, dict, [[id: 'no_dbsnp'], []], [[id: 'no_dbsnp_tbi'], []])
    versions = versions.mix(GATK4_HAPLOTYPECALLER.out.versions)

    // GExtract scaffolds from reference FASTA
    ch_scaffolds = ref_fasta.map { tuple -> tuple[1] }
        .map { path -> file(path).readLines() }
        .flatten()
        .filter { line -> line.startsWith(">") }
        .map { line -> line.split()[0][1..-1] }

    // Join VCF and TBI channels
    ch_vcf_and_tbi = GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi)
        .map { meta, vcf, tbi -> tuple(meta, vcf, tbi) }
        .toList()
        .map { records ->
            def vcfs = records.collect { tuple -> tuple[1] } // vcf paths
            def tbis = records.collect { tuple -> tuple[2] } // tbi paths
            tuple(vcfs, tbis)
        }

    // Prepare input for GATK GenomicsDBImport
    ch_gdb_input = ch_scaffolds
        .map { scaffold -> [id: scaffold] }
        .combine(ch_vcf_and_tbi)
        .map { meta, vcfs, tbis ->
            tuple(meta, vcfs, tbis, [], meta.id, file('.'))
        }

    // Run GATK GenomicsDBImport
    GATK4_GENOMICSDBIMPORT(ch_gdb_input, false, false, false)
    versions = versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)

    // Run GATK GenotypeGVCFs
    ch_gtp_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map { meta, gdb_path -> tuple(meta, gdb_path, [], [], []) }
    GATK4_GENOTYPEGVCFS(ch_gtp_input, ref_fasta, fai, dict, [[id: 'no_dbsnp'], []], [[id: 'no_dbsnp_tbi'], []])
    versions = versions.mix(GATK4_GENOTYPEGVCFS.out.versions)

    // Run BCFtools stats
    GATK4_GENOTYPEGVCFS.out.vcf.join(GATK4_GENOTYPEGVCFS.out.tbi)
    .map { meta, vcf, tbi -> tuple(meta, vcf, tbi) }
    .set { ch_vcf_tbi }
    BCFTOOLS_STATS(ch_vcf_tbi, [[id: 'no_regions'], []], [[id: 'no_targets'], []], [[id: 'no_samples'], []], [[id: 'no_exons'], []], ref_fasta)
    multiqc_files = multiqc_files.mix(BCFTOOLS_STATS.out.stats.map { tuple -> tuple[1] })
    versions = versions.mix(BCFTOOLS_STATS.out.versions)

    emit:
    vcf = GATK4_GENOTYPEGVCFS.out.vcf
    tbi = GATK4_GENOTYPEGVCFS.out.tbi
    multiqc_files
    versions
}
