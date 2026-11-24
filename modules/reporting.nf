
process REPORT {

    tag "$sample_id"

    conda "${projectDir}/envs/scrna.yml"

    publishDir { "${params.outdir}/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(seurat_rds), path(markers_csv)

    output:
        tuple val(sample_id),
        path("${sample_id}_umap.png"),
        path("${sample_id}_qc_violin.png"),
          path("${sample_id}_marker_heatmap.png"),
        path("${sample_id}_summary.txt"),
        emit: report

    script:
    """
    Rscript ${projectDir}/scripts/reporting.R \
        ${seurat_rds} \
        ${markers_csv} \
        . \
        ${sample_id}
    """
}
