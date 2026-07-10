/*
 * ML novelty classifier.
 *
 * Inputs:  contigs with an RdRpCATCH hit
 * Outputs: per-ORF embedding + kNN call against a precomputed
 *          embedding database of known viral RdRps (family-labelled).
 *
 * Tier-3 logic in reconcile.py uses the cosine similarity to the
 * nearest neighbour to decide whether to accept a novel-virus call.
 */
process ML_NOVELTY {
    tag "$meta.id"
    label 'gpu'
    publishDir "${params.outdir}/04_classification/ml_novelty/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(rdrpcatch_hits), path(contigs)

    output:
    tuple val(meta), path("${meta.id}_ml_calls.tsv"), emit: calls
    tuple val(meta), path("${meta.id}_embeddings.npz"), optional: true, emit: embeddings

    script:
    """
    python3 ${projectDir}/bin/embed_and_knn.py \\
        --contigs ${contigs} \\
        --hmm_hits ${rdrpcatch_hits} \\
        --reference_db ${params.esm2_embeddings} \\
        --k ${params.knn_neighbours} \\
        --output ${meta.id}_ml_calls.tsv \\
        --embeddings_out ${meta.id}_embeddings.npz \\
        --sample ${meta.id}
    """
}
