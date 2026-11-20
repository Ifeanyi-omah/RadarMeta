process NAMING_AUDIT {
    tag "$meta.id"
    label 'small'
    publishDir "${params.outdir}/05_reconciled/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(calls)
    path vmr
    path history      // assets/NO_FILE when not supplied
    path literature   // assets/NO_FILE when not supplied

    output:
    tuple val(meta), path("${meta.id}_name_audit.tsv"),          emit: audit
    tuple val(meta), path("${meta.id}_name_audit_summary.json"), emit: summary
    tuple val(meta), path("${meta.id}_literature_rejected.tsv"), emit: rejected

    script:
    def hist_arg = history.name != 'NO_FILE' ? "--history ${history}" : ''
    def lit_arg  = literature.name != 'NO_FILE' ? "--literature ${literature}" : ''
    def rel_arg  = params.ictv_release ? "--ictv_release ${params.ictv_release}" : ''
    """
    python3 ${projectDir}/bin/audit_names.py \\
        --calls ${calls} \\
        --vmr ${vmr} \\
        ${rel_arg} ${hist_arg} ${lit_arg} \\
        --audit_out ${meta.id}_name_audit.tsv \\
        --summary_out ${meta.id}_name_audit_summary.json \\
        --rejected_out ${meta.id}_literature_rejected.tsv
    """
}
