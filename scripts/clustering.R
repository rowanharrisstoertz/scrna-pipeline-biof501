library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 7) {
    stop("Enter: input_file, output_dir, sample_name, findneighbors_dims, resolution, umap_dims, umap_min_dist")
}

input_file <- args[1]
output_dir <- args[2]
sample_name <- args[3]
findneighbors_dims <- as.integer(args[4])
resolution <- as.numeric(args[5])
umap_dims <- as.integer(args[6])
umap_min_dist <- as.numeric(args[7])

# ensure parameters are valid
if (anyNA(c(findneighbors_dims, resolution, umap_dims, umap_min_dist))) {
    stop("Clustering and UMAP parameters must be valid numeric values")
}
if (findneighbors_dims < 1 || umap_dims < 1) {
    stop("Dimension counts must be at least 1")
}
if (resolution <= 0) {
    stop("Clustering resolution must be greater than 0")
}
if (umap_min_dist < 0 || umap_min_dist > 1) {
    stop("UMAP min.dist must be between 0 and 1")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load the filtered Seurat object
seurat_obj <- readRDS(input_file)

# find neighbors using the specified number of principal components
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:findneighbors_dims)

cat("\nFindNeighbors completed for", sample_name, "using", findneighbors_dims, "principal components\n")

# perform clustering using the Louvain algorithm
seurat_obj <- FindClusters(seurat_obj, resolution = resolution)

cat("\nFindClusters completed for", sample_name, "with resolution", resolution, "\n")

# run UMAP
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:umap_dims, min.dist = umap_min_dist)

cat("\nUMAP completed for", sample_name, "with", umap_dims, "dimensions and min.dist =", umap_min_dist, "\n")

# Save clustered object
output_file <- file.path(
    output_dir,
    paste0(sample_name, "_clustered.rds")
)
saveRDS(seurat_obj, output_file)

cat("\nClustering and UMAP completed for", sample_name, "\n")