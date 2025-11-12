process ANGSD_NGSRELATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ngsrelate:2.0--hea85c65_0':
          'biocontainers/ngsrelate:2.0--hea85c65_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.ngsrelate.results"), emit: plots
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id
    def suffix = "ngsrelate.results"
    """
    ngsRelate \\
      -O ${prefix}.${suffix}  \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngsrelate: "2.0"
    END_VERSIONS
    """
}
