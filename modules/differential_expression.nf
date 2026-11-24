
process DIFFERENTIAL_EXPRESSION {

    tag "$sample_id"

    conda "${projectDir}/envs/scrna.yml"

    input:
    tuple val(sample_id), path(seurat_rds)

    output:
    tuple val(sample_id),
          path("${sample_id}_cluster_markers.csv"),
          emit: markers

    script:
    """
    Rscript ${projectDir}/scripts/differential_expression.R ${seurat_rds} . ${sample_id}
    """
}
