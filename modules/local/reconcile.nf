process RECONCILE {
    tag "$meta.id"
    label 'small'
    publishDir "${params.outdir}/05_reconciled/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(contigs), path(rdrpcatch), path(diamond), path(blastnt), path(ml_calls)
    path pathogens_list

    output:
    tuple val(meta), path("${meta.id}_calls.tsv"),         emit: calls
    tuple val(meta), path("${meta.id}_viral_genomes.fa"),  optional: true, emit: viral_fa
    tuple val(meta), path("${meta.id}_alerts.json"),       optional: true, emit: alerts
    tuple val(meta), path("${meta.id}_novel.tsv"),         optional: true, emit: novel

    script:
    """
    python3 ${projectDir}/bin/reconcile.py \\
        --sample ${meta.id} \\
        --contigs ${contigs} \\
        --rdrpcatch ${rdrpcatch} \\
        --diamond ${diamond} \\
        --blastnt ${blastnt} \\
        --ml_calls ${ml_calls} \\
        --pathogens ${pathogens_list} \\
        --knn_min_cos_sim ${params.knn_min_cos_sim} \\
        --calls_out ${meta.id}_calls.tsv \\
        --viral_fa_out ${meta.id}_viral_genomes.fa \\
        --alerts_out ${meta.id}_alerts.json \\
        --novel_out ${meta.id}_novel.tsv
    """
}
