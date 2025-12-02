/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PLINK  } from '../../../modules/local/plink/'
include { READv2 } from '../../../modules/local/readv2'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow RELATEDNESS_READ {
    take:
    vcf // channel: [ meta, vcf ]

    main:
    versions = channel.empty()
    multiqc_files = channel.empty()


    vcf.dump(tag: 'RELATEDNESS_READ (vcf)')
    PLINK(vcf)
    versions = versions.mix(PLINK.out.versions)

    input_read = PLINK.out.bed.dump(tag: 'RELATEDNESS_READ (PLINK.out.bed))')
      .join(PLINK.out.bim).dump(tag: 'RELATEDNESS_READ (VCFTOOLS.out.bed.join(VCFTOOLS.out.bim))')
      .join(PLINK.out.fam).dump(tag: 'RELATEDNESS_READ (input_read)')

    READv2(input_read)
    versions = versions.mix(READv2.out.versions)

    READv2.out.tsv.dump(tag: 'RELATEDNESS_READ (READv2.out.tsv)')

    emit:
    read_tsv = READv2.out.tsv
    read_pdf = READv2.out.pdf
    multiqc_files
    versions
}
