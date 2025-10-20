/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: SAMPLESHEET_TO_CHANNEL
    Reads a user-provided samplesheet (CSV) and emits a channel with tuples of:
    [ meta, [file1, file2] ]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SAMPLESHEET_TO_CHANNEL {
    take:
    samplesheet

    main:
    ch_samplesheet = Channel
        .of(file(samplesheet).text)
        .splitCsv(header:true)
        .map { row ->
            def meta = [
                id: row.sample,
                RGID: row.RGID,
                RGLB: row.RGLB,
                RGPL: row.RGPL,
                RGPU: row.RGPU,
                RGSM: row.RGSM
            ]
            tuple(meta, [ file(row.fastq_1), file(row.fastq_2) ])
        }

    emit:
    ch_samplesheet
}
