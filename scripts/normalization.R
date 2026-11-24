library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
    stop("Enter: input_file, output_dir, sample_name, n_pcs")
}

input_file <- args[1]
output_dir <- args[2]
sample_name <- args[3]
n_pcs <- as.integer(args[4])

if (is.na(n_pcs) || n_pcs < 1) {
    stop("n_pcs must be a positive integer")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load the filtered Seurat object
seurat_obj <- readRDS(input_file)

# Normalize and scale the data and find variable features using SCTransform
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)

cat("\nNormalization and scaling completed for", sample_name, "\n")

# Run PCA on the Seurat object
seurat_obj <- RunPCA(seurat_obj, npcs = n_pcs, assay = "SCT", seed.use = 42)

cat("\nPCA completed for", sample_name, "with", n_pcs, "principal components\n")

# Save the normalized Seurat object to an RDS file
output_file <- file.path(output_dir, paste0(sample_name, "_normalized.rds"))
saveRDS(seurat_obj, output_file)

cat("\nNormalization and PCA completed for", sample_name, "\n")