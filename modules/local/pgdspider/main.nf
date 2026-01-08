process PGDSPIDER {
    tag "$meta.id"
    label 'process_small'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pgdspider%3A2.1.1.5--hdfd78af_1' :
        'community.wave.seqera.io/library/pgdspider:2.1.1.5--5f780e1a319ff195' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.ind") , emit: ind
    tuple val(meta), path("*.snp") , emit: snp
    tuple val(meta), path("*.geno"), emit: geno
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def file_ending = args.contains("-outputformat EIGENSOFT") ? "eigenstrat" : ""

    """
    PGDSpider2-cli \\
    -inputfile ${vcf} \\
    -inputformat VCF \\
    -outputfile ${meta.id}.${file_ending} \\
    ${args} \\
    -spid ${workflow.projectDir}/assets/pgdspider_vcf_to_eigenstrat.spid

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PGDSpider: "2.1.1.5"
    END_VERSIONS
    """
}
