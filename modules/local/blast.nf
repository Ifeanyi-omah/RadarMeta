process BLAST_NT {
    tag "$meta.id"
    label 'xlarge_mem'
    publishDir "${params.outdir}/04_classification/blast_nt/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}_blast_nt.m9"), emit: hits

    script:
    """
    export BLASTDB=\$(dirname ${params.blast_nt})

    # chunk to keep BLAST memory bounded
    seqkit split2 -p 10 ${contigs} -O chunks

    for chunk in chunks/*.fa; do
        blastn -task dc-megablast \\
            -db \$(basename ${params.blast_nt}) \\
            -num_threads ${task.cpus} \\
            -query \$chunk \\
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname stitle" \\
            -evalue 1e-5 \\
            -out \${chunk}.m9
    done

    cat chunks/*.m9 > ${meta.id}_blast_nt.m9
    rm -rf chunks
    """
}
