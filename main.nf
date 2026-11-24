
include { PREPROCESS } from './modules/preprocessing'
include { NORMALIZE } from './modules/normalization'
include { CLUSTER } from './modules/clustering'
include { DIFFERENTIAL_EXPRESSION } from './modules/differential_expression'
include { ANNOTATION } from './modules/annotation'
include { REPORT } from './modules/reporting'

workflow {
    // check parameters
    if (params.findneighbors_dims.toInteger() > params.n_pcs.toInteger()) {
        error "findneighbors_dims cannot be greater than n_pcs"
    }
    if (params.umap_dims.toInteger() > params.n_pcs.toInteger()) {
        error "umap_dims cannot be greater than n_pcs"
    }
    if (params.resolution.toDouble() <= 0) {
        error "Clustering resolution must be greater than 0"
    }
    if (params.umap_min_dist.toDouble() < 0 || params.umap_min_dist.toDouble() > 1) {
        error "umap_min_dist must be between 0 and 1"
    }
    if (params.annotation_min_score.toDouble() < 0) {
        error "annotation_min_score must be non-negative"
    }

    if (params.annotation_min_difference.toDouble() < 0) {
        error "annotation_min_difference must be non-negative"
    }

        input_files = [
        "${params.input_dir}/${params.sample_id}/matrix.mtx.gz",
        "${params.input_dir}/${params.sample_id}/features.tsv.gz",
        "${params.input_dir}/${params.sample_id}/barcodes.tsv.gz"
    ]

    missing_files = input_files.findAll { !file(it).exists() }

    if (missing_files) {
        error "Missing input files for sample ${params.sample_id}: ${missing_files.join(', ')}"
    }
    raw_input = Channel.of(
        tuple(
            params.sample_id,
            file("${params.input_dir}/${params.sample_id}/matrix.mtx.gz"),
            file("${params.input_dir}/${params.sample_id}/features.tsv.gz"),
            file("${params.input_dir}/${params.sample_id}/barcodes.tsv.gz")
        )
    )

    PREPROCESS(raw_input)

    normalized_input = PREPROCESS.out.processed.map {
        sample_id, rds ->
        tuple(sample_id, rds)
    }

    NORMALIZE(normalized_input)

    cluster_input = NORMALIZE.out.normalized

    CLUSTER(cluster_input)

    marker_reference = file(params.marker_file)
    annotation_script = file("${projectDir}/scripts/annotation.R")

    annotation_input = CLUSTER.out.clustered.map {
        sample_id, rds ->
        tuple(sample_id, rds, marker_reference, annotation_script)
    }

    ANNOTATION(annotation_input)

    marker_input = CLUSTER.out.clustered.map {
        sample_id, rds ->
        tuple(sample_id, rds)
    }

    DIFFERENTIAL_EXPRESSION(marker_input)

report_input = CLUSTER.out.clustered.map {
    sample_id, rds ->
    tuple(sample_id, rds)
}.join(
    DIFFERENTIAL_EXPRESSION.out.markers
)

REPORT(report_input)
}