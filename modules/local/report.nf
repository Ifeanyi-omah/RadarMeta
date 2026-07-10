process REPORT {
    label 'small'
    publishDir "${params.outdir}/06_report", mode: 'copy'

    input:
    path call_files

    output:
    path "cohort_summary.html"
    path "per_sample/*.html"

    script:
    """
    mkdir -p per_sample
    python3 ${projectDir}/bin/render_report.py \\
        --calls ${call_files} \\
        --cohort_out cohort_summary.html \\
        --per_sample_dir per_sample
    """
}
