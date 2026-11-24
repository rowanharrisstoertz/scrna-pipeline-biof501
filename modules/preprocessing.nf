
process PREPROCESS {

    tag "$sample_id"

    conda "${projectDir}/envs/scrna.yml"

    input:
    tuple val(sample_id), path(matrix), path(features), path(barcodes)

    output:
    tuple val(sample_id),
          path("${sample_id}_filtered.rds"),
          emit: processed

    script:
    """
    Rscript ${projectDir}/scripts/preprocessing.R \
        . . ${sample_id} \
        ${params.min_features} \
        ${params.max_features} \
        ${params.max_percent_mt}
    """
}
