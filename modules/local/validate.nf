process VALIDATE_SAMPLESHEET {
    label 'small'

    input:
    path samplesheet

    output:
    path 'samplesheet_validated.csv'

    script:
    """
    python3 ${projectDir}/bin/validate_samplesheet.py \\
        --input ${samplesheet} \\
        --output samplesheet_validated.csv \\
        --from_contigs ${params.from_contigs}
    """
}
