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
    # Sort scaffolds by length (column 3) descending
    sort -k3,3nr ${intervals} | \
    awk -v N=${params.target_number_of_intervals} '
        {
            bin = (NR - 1) % N + 1
            bins[bin] = bins[bin] sprintf("%s\\t%s\\t%s\\n", \$1, \$2, \$3)
        }
        END {
            for (i = 1; i <= N; i++) {
                if (bins[i] != "") {
                    fname = sprintf("interval_%0${pad_width}d.bed", i)
                    printf "%s", bins[i] > fname
                }
            }
        }
    '

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version | sed '1!d; s/.* //')
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
