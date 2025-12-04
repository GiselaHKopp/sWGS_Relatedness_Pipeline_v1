process GATK4_ANALYZECOVARIATES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(before_table), path(after_table), path(table3)

    output:
    tuple val(meta), path("${meta.id}.pdf"), emit: plots
    tuple val(meta), path("${meta.id}.csv"), emit: data
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def third_table = table3 ? "-bqsr ${table3}" : ""

    """
    gatk AnalyzeCovariates \\
      -before ${before_table} \\
      -after ${after_table} \\
      ${third_table} \\
      -csv ${meta.id}.csv \\
      -plots ${meta.id}.pdf \\
      ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk AnalyzeCovariates --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
