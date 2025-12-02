process READv2 {
    tag "$meta.id"
    label 'process_small'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kinship-read:2.1.1--pyh7cba7a3_0' :
        'community.wave.seqera.io/library/kinship-read:2.1.1--69f4dcbacf410d9c' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("Read_Results.tsv")     , emit: tsv
    tuple val(meta), path("READ_results_plot.pdf"), emit: pdf
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}"
    def prefix = bed.getName() - ".bed"
    println "prefix   = ${prefix}"
    """
    kinship-read \\
        -i ${prefix} \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        READv2: \$(READ2 --version 2>&1 | head -n1 | sed 's/^.*READv2 //; s/ .*\$//')
    END_VERSIONS
    """
}
