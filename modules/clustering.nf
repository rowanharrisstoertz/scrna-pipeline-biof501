
process CLUSTER {

    tag "$sample_id"

    conda "${projectDir}/envs/scrna.yml"

    input:
    tuple val(sample_id), path(seurat_rds)

    output:
    tuple val(sample_id),
          path("${sample_id}_clustered.rds"),
          emit: clustered

    script:
        """
        Rscript ${projectDir}/scripts/clustering.R \
            ${seurat_rds} \
            . \
            ${sample_id} \
            ${params.findneighbors_dims} \
            ${params.resolution} \
            ${params.umap_dims} \
            ${params.umap_min_dist}
        """
}
