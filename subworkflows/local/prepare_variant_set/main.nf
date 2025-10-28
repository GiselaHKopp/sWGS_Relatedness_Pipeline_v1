include { BCFTOOLS_STATS          } from '../../../modules/nf-core/bcftools/stats'
include { GATK4_GENOMICSDBIMPORT  } from '../../../modules/nf-core/gatk4/genomicsdbimport'
include { GATK4_GENOTYPEGVCFS     } from '../../../modules/nf-core/gatk4/genotypegvcfs'
include { GATK4_HAPLOTYPECALLER   } from '../../../modules/nf-core/gatk4/haplotypecaller'
include { GATK4_SELECTVARIANTS    } from '../../../modules/nf-core/gatk4/selectvariants'
include { GATK4_VARIANTFILTRATION } from '../../../modules/nf-core/gatk4/variantfiltration'

workflow PREPARE_VARIANT_SET {
    take:
    bam          // tuple(meta, path_to_bam)         e.g. [ id: 'sample1' ], sample1.bam
    bam_index    // tuple(meta, path_to_bai)         e.g. [ id: 'sample1' ], sample1.bam.bai
    fai          // tuple(meta, path_to_fasta.fai)   e.g. [ id: 'ref' ], ref.fasta.fai
    dict         // tuple(meta, path_to_dict)        e.g. [ id: 'ref' ], ref.dict
    ref_fasta    // tuple(meta2, path_to_fasta)      e.g. [ id: 'ref' ], ref.fasta

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
    ch_vcf_tbi = GATK4_GENOTYPEGVCFS.out.vcf.join(GATK4_GENOTYPEGVCFS.out.tbi).map { meta, vcf, tbi -> tuple(meta, vcf, tbi) }
    BCFTOOLS_STATS(ch_vcf_tbi, [[id: 'no_regions'], []], [[id: 'no_targets'], []], [[id: 'no_samples'], []], [[id: 'no_exons'], []], ref_fasta)
    multiqc_files = multiqc_files.mix(BCFTOOLS_STATS.out.stats.map { tuple -> tuple[1] })
    versions = versions.mix(BCFTOOLS_STATS.out.versions)

    // Append "filtered" to avoid input/output name clashes
    ch_filtered_input = ch_vcf_tbi.map { meta, vcf, tbi ->
        def new_meta = meta.clone()
        new_meta.id = "${meta.id}.filtered"
        tuple(new_meta, vcf, tbi)
    }

    // Filter variants to exclude low-quality calls
    GATK4_VARIANTFILTRATION(ch_filtered_input, ref_fasta, fai, dict, [[id: 'no_gzi'], []])
    versions = versions.mix(GATK4_VARIANTFILTRATION.out.versions)

/*
    // Select only passing variants
    GATK4_SELECTVARIANTS(GATK4_VARIANTFILTRATION.out.vcf)


    // Gather VCFs into single file
    GATK4_GATHERVCFS(GATK4_SELECTVARIANTS.out.vcf)
*/
    emit:
    vcf = GATK4_HAPLOTYPECALLER.out.vcf
    tbi = GATK4_HAPLOTYPECALLER.out.tbi
    multiqc_files
    versions
    ch_scaffolds
}
