process BCFTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data':
        'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f' }"

    input:
    tuple val(meta), path(intervals)
    val(bams)
    tuple val(meta3), path(fasta)
    val save_mpileup   // boolean

    output:
    tuple val(meta), path("${meta.id}.mpileup.bcf")     , emit: bcf, optional: true
    tuple val(meta), path("${meta.id}.mpileup.vcf.gz")  , emit: vcf, optional: true
    tuple val(meta), path("${meta.id}.mpileup.gz")      , emit: raw_mpileup, optional: true
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix       = task.ext.prefix ?: meta.num_intervals <= 1 ? meta.id : meta.id + meta.interval_name
    def args         = task.ext.args ?: ''
    def intervalsOpt = intervals ? "-T ${intervals}" : ""
    bam_args = bams.collect{ tuple -> tuple.toString() }.join(" ")

    // Allowed output types: b, u, z, v, v0 … etc.
    def outputType = task.ext.output_type ?: "z"

    // Optional saving of the raw textual mpileup
    def raw_mpileup_cmd = save_mpileup ? "| tee ${prefix}.mpileup" : ""
    def compress_raw_mpileup = save_mpileup ? "bgzip -f ${prefix}.mpileup" : ""

    // Output filename depends on format
    def outputFile =
        (outputType == 'b') ? "${prefix}.mpileup.bcf" :
        (outputType.startsWith('v') || outputType.startsWith('z')) ? "${prefix}.mpileup.vcf.gz" :
        "${prefix}.mpileup.bcf"

    """
    bcftools mpileup \
        --fasta-re ${fasta} \
        ${intervalsOpt} \
        ${args} \
        -O ${outputType} \
        -o ${outputFile} \
        ${bam_args} \
        ${raw_mpileup_cmd}

    ${compress_raw_mpileup}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
