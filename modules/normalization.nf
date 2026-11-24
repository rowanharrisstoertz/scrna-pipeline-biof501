
process NORMALIZE {

    tag "$sample_id"

    conda "${projectDir}/envs/scrna.yml"

    input:
    tuple val(sample_id), path(seurat_rds)

    output:
    tuple val(sample_id),
          path("${sample_id}_normalized.rds"),
          emit: normalized

    script:
    """
    Rscript ${projectDir}/scripts/normalization.R \
        ${seurat_rds} \
        . \
        ${sample_id} \
        ${params.n_pcs}
    """
}
