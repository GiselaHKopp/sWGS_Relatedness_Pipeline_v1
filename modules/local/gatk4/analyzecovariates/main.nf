process GATK4_ANALYZECOVARIATES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(table1), path(table2)

    output:
    tuple val(meta), path("*.pdf"),  emit: plots
    tuple val(meta), path("*.csv"),  emit: data
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gatk AnalyzeCovariates \\
      -before ${table1} \\
      -after ${table2} \\
      -csv ${meta.id}.csv \\
      -plots ${meta.id}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk AnalyzeCovariates --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
