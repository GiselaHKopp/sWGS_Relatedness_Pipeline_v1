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
    versions = Channel.empty()
    multiqc_files = Channel.empty()

    // --- Combine BAM + BAI ---
    ch_inputs = bam.join(bam_index)
        .map { meta, bamfile, baifile ->
            tuple(meta, bamfile, baifile, [], [])
        }

    GATK4_HAPLOTYPECALLER(ch_inputs, ref_fasta, fai, dict, [[id: 'no_dbsnp'], []], [[id: 'no_dbsnp_tbi'], []])
    versions = versions.mix(GATK4_HAPLOTYPECALLER.output.versions)

    ch_scaffolds = ref_fasta.map { it[1] }
        .map { file(it).readLines() }
        .flatten()
        .filter { it.startsWith(">") }
        .map { it.split()[0][1..-1] }
        .dump()
    /*
    GATK4_GENOMICSDBIMPORT(haplotype_vcf_and_tbi, false, false, false)
    versions = versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)

    // --- 3. Joint genotyping ---
    GATK4_GENOTYPEGVCFS(GATK4_GENOMICSDBIMPORT.out.workspace, ref_fasta)
    versions = versions.mix(GATK4_GENOTYPEGVCFS.out.versions)

    // --- 4. Hard filtering ---
    GATK4_VARIANTFILTRATION(GATK4_GENOTYPEGVCFS.out.vcf)
    GATK4_SELECTVARIANTS(GATK4_VARIANTFILTRATION.out.vcf)

    // --- 5. Gather to single VCF ---
    GATK4_GATHERVCFS(GATK4_SELECTVARIANTS.out.vcf)
*/
    emit:
    vcf = GATK4_HAPLOTYPECALLER.out.vcf
    tbi = GATK4_HAPLOTYPECALLER.out.tbi
    multiqc_files
    versions
    ch_scaffolds
}
