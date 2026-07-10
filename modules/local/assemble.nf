process ASSEMBLE_MEGAHIT {
    tag "$meta.id"
    label 'large'
    publishDir "${params.outdir}/03_assembly/${meta.id}", mode: 'copy', pattern: 'final.contigs.fa'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.contigs.fa"), emit: contigs
    path "${meta.id}_megahit.log",                  emit: log

    script:
    def r1 = reads[0]; def r2 = reads[1]
    """
    megahit \\
        -1 ${r1} -2 ${r2} \\
        -o ${meta.id}_megahit \\
        --num-cpu-threads ${task.cpus} \\
        --min-contig-len ${params.min_contig_len} \\
        2>&1 | tee ${meta.id}_megahit.log

    # rename with sample prefix and put at top level
    awk -v s=${meta.id} '/^>/{print ">"s"_"substr(\$0,2); next}{print}' \\
        ${meta.id}_megahit/final.contigs.fa > ${meta.id}.contigs.fa

    # tar intermediates to save space
    tar -czf ${meta.id}_megahit_intermediates.tar.gz \\
        -C ${meta.id}_megahit intermediate_contigs 2>/dev/null || true
    rm -rf ${meta.id}_megahit/intermediate_contigs
    """
}

process MAP_BACK {
    tag "$meta.id"
    label 'medium'
    publishDir "${params.outdir}/03_assembly/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(contigs), path(reads)

    output:
    tuple val(meta), path("${meta.id}_coverage.tsv"), emit: coverage

    script:
    def r1 = reads[0]; def r2 = reads[1]
    """
    bowtie2-build ${contigs} ${meta.id}_idx
    bowtie2 -x ${meta.id}_idx -1 ${r1} -2 ${r2} \\
        --very-sensitive --threads ${task.cpus} 2> ${meta.id}_mapback.log \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.bam -
    samtools index ${meta.id}.bam
    samtools idxstats ${meta.id}.bam | \\
        awk 'BEGIN{OFS="\\t"; print "contig_id","length","mapped","unmapped"} {print}' \\
        > ${meta.id}_coverage.tsv
    """
}
