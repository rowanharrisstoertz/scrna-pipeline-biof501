library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop("Enter: input_dir, output_dir, sample_name")
}

input_file <- args[1]
output_dir <- args[2]
sample_name <- args[3]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load the filtered Seurat object
seurat_obj <- readRDS(input_file)

# set the default assay to SCT
DefaultAssay(seurat_obj) <- "SCT"
seurat_obj <- PrepSCTFindMarkers(seurat_obj)

# find markers
markers <- FindAllMarkers(seurat_obj, assay = "SCT", only.pos = TRUE)

# Save markers to CSV
output_file <- file.path(output_dir, paste0(sample_name, "_cluster_markers.csv"))
write.csv(markers, output_file, row.names = FALSE)

cat("\nDifferential expression analysis completed for", sample_name, "\n")
