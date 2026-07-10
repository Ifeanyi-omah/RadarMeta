process DIAMOND_BLASTX {
    tag "$meta.id"
    label 'large'
    publishDir "${params.outdir}/04_classification/diamond_nr/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}_diamond_nr.m9"), emit: hits

    script:
    """
    diamond blastx \\
        -d ${params.diamond_nr} \\
        -q ${contigs} \\
        -o ${meta.id}_diamond_nr.m9 \\
        --evalue ${params.diamond_evalue} \\
        --threads ${task.cpus} \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen \\
                   qstart qend sstart send evalue bitscore \\
                   staxids sscinames stitle
    """
}
