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
    # Number of target bins
    N=${params.target_number_of_intervals}

    # Sort contigs by length (column 3) descending
    sort -k3,3nr ${intervals} | \
    awk -v N=\$N '
        BEGIN {
            # initialize bins
            for (i = 1; i <= N; i++) {
                sum[i] = 0
                data[i] = ""
            }
        }

        {
            # find currently lightest bin
            best = 1
            for (i = 2; i <= N; i++) {
                if (sum[i] < sum[best])
                    best = i
            }

            # assign contig to that bin
            data[best] = data[best] sprintf("%s\\t%s\\t%s\\n", \$1, \$2, \$3)
            sum[best] += \$3
        }

        END {
            for (i = 1; i <= N; i++) {
                if (data[i] != "") {
                    fname = sprintf("interval_%0${pad_width}d.bed", i)
                    printf "%s", data[i] > fname
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
