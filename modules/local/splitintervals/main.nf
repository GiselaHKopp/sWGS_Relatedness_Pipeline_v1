process SPLIT_INTERVALS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(intervals)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml",            emit: versions

    script:
    def pad_width = params.target_number_of_intervals.toString().size()

    """
    # Number of scaffolds
    S=\$(wc -l < ${intervals})

    # Target number of interval files
    N=${params.target_number_of_intervals}

    # Ceil(S / N)
    K=\$(( (S + N - 1) / N ))

    awk -v K=\$K '
        {
            file = int((NR - 1) / K) + 1
            printf "%s\\t%s\\t%s\\n", \$1, \$2, \$3 >> sprintf("interval_%0${pad_width}d.bed", file)
        }
    ' ${intervals}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
