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
    N=${params.target_number_of_intervals}
    MIN_WINDOW=1000000   # 1 Mb safeguard

    TOTAL=\$(awk '{ sum += \$3 - \$2 } END { print sum }' ${intervals})
    TARGET=\$(( TOTAL / N ))
    WINDOW=\$(( TARGET / 5 ))
    if [ \$WINDOW -lt \$MIN_WINDOW ]; then WINDOW=\$MIN_WINDOW; fi

    awk -v TARGET=\$TARGET -v WINDOW=\$WINDOW -v N=\$N '
        BEGIN {
            idx = 1
            chunk = 0
            fname = sprintf("interval_%0${pad_width}d.bed", idx)
        }

        {
            contig = \$1
            start  = \$2
            end    = \$3

            for (pos = start; pos < end; pos += WINDOW) {
                win_start = pos
                win_end   = (pos + WINDOW < end) ? pos + WINDOW : end
                len = win_end - win_start

                if (chunk > 0 && (chunk + len) > TARGET && idx < N) {
                    idx++
                    chunk = 0
                    fname = sprintf("interval_%0${pad_width}d.bed", idx)
                }

                printf "%s\\t%d\\t%d\\n", contig, win_start, win_end >> fname
                chunk += len
            }
        }
    ' ${intervals}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
