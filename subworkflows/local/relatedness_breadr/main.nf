/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//include { BREADR    } from '../../../modules/local/breadr/'
include { PGDSPIDER } from '../../../modules/local/pgdspider/'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow RELATEDNESS_BREADR {
    take:
    vcf // channel: [ meta, vcf ]

    main:
    versions = channel.empty()
    multiqc_files = channel.empty()


    vcf.dump(tag: 'RELATEDNESS_BREADR (vcf)')
    PGDSPIDER(vcf)
    versions = versions.mix(PGDSPIDER.out.versions)

    PGDSPIDER.out.ind.dump(tag: 'RELATEDNESS_BREADR (PGDSPIDER.out.ind))')

    //BREADR(input_read)
    //versions = versions.mix(BREADR.out.versions)

    emit:
    test = PGDSPIDER.out.ind
    multiqc_files
    versions
}
