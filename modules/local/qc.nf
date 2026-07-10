process FASTP {
    tag "$meta.id"
    label 'medium'
    publishDir "${params.outdir}/01_qc/fastp", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed_R{1,2}.fastq.gz"), emit: reads
    path "*.json",                                      emit: json
    path "*.html",                                      emit: html
    path "*.log",                                       emit: log

    script:
    def r1 = reads[0]; def r2 = reads[1]
    """
    fastp \\
        -i ${r1} -I ${r2} \\
        -o ${meta.id}_trimmed_R1.fastq.gz \\
        -O ${meta.id}_trimmed_R2.fastq.gz \\
        --detect_adapter_for_pe \\
        --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \\
        --length_required 50 \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html 2> ${meta.id}.fastp.log
    """
}

process MULTIQC {
    label 'small'
    publishDir "${params.outdir}/01_qc/multiqc", mode: 'copy'

    input:
    path '*'

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc . --force
    """
}
