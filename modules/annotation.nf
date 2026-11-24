
process ANNOTATION {

    tag "$sample_id"

    conda "${projectDir}/envs/scrna.yml"

    publishDir { "${params.outdir}/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id),
          path(seurat_rds),
          path(marker_file),
          path(annotation_script)

output:
tuple val(sample_id),
      path("${sample_id}_annotation_scores.csv"),
      path("${sample_id}_annotated_umap.png"),
      path("${sample_id}_annotated.rds"),
      emit: diagnostics

    script:
    """
    Rscript ${annotation_script} \
        ${seurat_rds} \
        . \
        ${sample_id} \
        ${marker_file} \
        ${params.annotation_min_score} \
        ${params.annotation_min_difference}
    """
}
