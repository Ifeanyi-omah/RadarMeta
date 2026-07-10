/*
 * Host-removal policy:
 *   human       → aggressive (bowtie2 GRCh38 + Kraken2 human) + verification
 *   vertebrate  → tag host, retain reads for downstream species ID
 *   invertebrate→ no-op (handled in main.nf by direct pass-through)
 *   wastewater  → no-op
 *
 * The verification step is critical for clinical samples: we will not
 * release reads from a human sample unless residual host fraction is
 * below params.dehost_verify_thresh.
 */

process DEHOST_HUMAN {
    tag "$meta.id"
    label 'large'
    publishDir "${params.outdir}/02_dehost/${meta.id}", mode: 'copy', pattern: '*.{log,tsv}'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_dehosted_R{1,2}.fastq.gz"), emit: reads
    path "${meta.id}_dehost_stats.tsv",                            emit: stats
    path "${meta.id}_dehost.log",                                  emit: log

    script:
    def r1 = reads[0]; def r2 = reads[1]
    """
    set -euo pipefail

    # 1. bowtie2 vs GRCh38, keep unmapped pairs
    bowtie2 \\
        -x ${params.host_grch38} \\
        -1 ${r1} -2 ${r2} \\
        --very-sensitive --threads ${task.cpus} \\
        --un-conc-gz ${meta.id}_step1_R%.fastq.gz \\
        -S /dev/null 2>> ${meta.id}_dehost.log

    # 2. second-pass Kraken2 against human DB to catch divergent host reads
    kraken2 --db ${params.kraken2_db}/human \\
        --threads ${task.cpus} --gzip-compressed --paired \\
        --unclassified-out ${meta.id}_unclass_R#.fastq \\
        --classified-out   ${meta.id}_humanhit_R#.fastq \\
        --output /dev/null --report ${meta.id}_kraken_host.report \\
        ${meta.id}_step1_R1.fastq.gz ${meta.id}_step1_R2.fastq.gz \\
        2>> ${meta.id}_dehost.log

    # 3. shuffle read names to break re-identification (Scylla pattern)
    python3 ${projectDir}/bin/shuffle_read_names.py \\
        --in1 ${meta.id}_unclass_R_1.fastq \\
        --in2 ${meta.id}_unclass_R_2.fastq \\
        --out1 ${meta.id}_dehosted_R1.fastq.gz \\
        --out2 ${meta.id}_dehosted_R2.fastq.gz \\
        --sample ${meta.id}

    # 4. verify residual host fraction
    python3 ${projectDir}/bin/verify_dehost.py \\
        --kraken ${meta.id}_kraken_host.report \\
        --threshold ${params.dehost_verify_thresh} \\
        --output ${meta.id}_dehost_stats.tsv \\
        --sample ${meta.id}

    # clean up intermediate uncompressed fastq
    rm -f ${meta.id}_step1_R*.fastq.gz ${meta.id}_unclass_R_*.fastq ${meta.id}_humanhit_R_*.fastq
    """
}

process HOST_TAG_VERTEBRATE {
    tag "$meta.id"
    label 'medium'
    publishDir "${params.outdir}/02_dehost/${meta.id}", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads),                           emit: reads
    path "${meta.id}_host_id.tsv", optional: true,          emit: host_id

    script:
    """
    # vertebrate samples: retain all reads, tag candidate host via mito BLAST
    # (downstream after assembly; here we just pass through and write a placeholder)
    echo -e "sample_id\\tnote" > ${meta.id}_host_id.tsv
    echo -e "${meta.id}\\thost retained for species ID after assembly" >> ${meta.id}_host_id.tsv
    """
}
