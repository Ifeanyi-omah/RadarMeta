process KRAKEN2_BRACKEN {
    tag "$meta.id"
    label 'large'
    publishDir "${params.outdir}/04_classification/kraken2/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.k2report"),         emit: report
    tuple val(meta), path("${meta.id}.bracken"),          emit: bracken
    tuple val(meta), path("${meta.id}.kraken2"),          emit: classified

    script:
    def r1 = reads[0]; def r2 = reads[1]
    """
    kraken2 --db ${params.kraken2_db} \\
        --threads ${task.cpus} \\
        --gzip-compressed --paired \\
        --confidence 0.1 \\
        --report-minimizer-data \\
        --report ${meta.id}.k2report.rich \\
        --output ${meta.id}.kraken2 \\
        ${r1} ${r2}

    # strip minimizer cols so Bracken can parse
    cut -f1-3,6- ${meta.id}.k2report.rich > ${meta.id}.k2report

    bracken -d ${params.kraken2_db} \\
        -i ${meta.id}.k2report \\
        -o ${meta.id}.bracken \\
        -w ${meta.id}.bracken.k2report \\
        -r ${params.bracken_readlen} -l S
    """
}
