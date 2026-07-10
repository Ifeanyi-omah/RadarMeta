process RDRPCATCH {
    tag "$meta.id"
    label 'large'
    publishDir "${params.outdir}/04_classification/rdrpcatch/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}_rdrpcatch_output_annotated.tsv"), emit: hits
    tuple val(meta), path("${meta.id}_rdrpcatch_orfs.faa"), optional: true, emit: orfs
    path "${meta.id}_rdrpcatch.log"

    script:
    """
    rdrpcatch scan \\
        -i ${contigs} \\
        -o ${meta.id}_rdrpcatch \\
        -db-dir ${params.rdrpcatch_db} \\
        -seq-type nuc \\
        -length-thr ${params.min_contig_len} \\
        -cpus ${task.cpus} \\
        --overwrite 2>&1 | tee ${meta.id}_rdrpcatch.log

    # surface the annotated table at the top level
    cp ${meta.id}_rdrpcatch/*_rdrpcatch_output_annotated.tsv \\
       ${meta.id}_rdrpcatch_output_annotated.tsv

    # collect predicted ORFs if RdRpCATCH emitted them
    if ls ${meta.id}_rdrpcatch/*orfs*.faa 1>/dev/null 2>&1; then
        cat ${meta.id}_rdrpcatch/*orfs*.faa > ${meta.id}_rdrpcatch_orfs.faa
    fi
    """
}
