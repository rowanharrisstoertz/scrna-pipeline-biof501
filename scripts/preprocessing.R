
#Processing and QC script for a single scRNA-seq sample

library(Seurat)

# Get directories and sample name
args <- commandArgs(trailingOnly = TRUE)
# Check if the correct number of arguments is provided
if (length(args) != 6) {
  stop("Enter: input_dir, output_dir, sample_name, min_features, max_features, max_percent_mt")
}
input_dir <- args[1]
output_dir <- args[2]
sample_name <- args[3]
min_features <- as.numeric(args[4])
max_features <- as.numeric(args[5])
max_percent_mt <- as.numeric(args[6])

if (!dir.exists(input_dir)) {
    stop("Input directory does not exist: ", input_dir)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load the raw count matrix from the input directory
cat("Loading raw count matrix\n")
counts <- Read10X(
  data.dir = input_dir,
  gene.column = 2,
  unique.features = TRUE
)

cat("\nOriginal matrix dimensions:", dim(counts), "\n")

# Create Seurat object from the raw counts
cat("Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = sample_name,
  min.cells = 0,
  min.features = 0
)



# Calculate Mitochondrial percentage
cat("Calculating mitochondrial percentage...\n")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Filter the Seurat object based on the filtering thresholds
filtered_obj <- subset(seurat_obj, subset = nFeature_RNA > min_features & nFeature_RNA < max_features & percent.mt < max_percent_mt)
cat("\nCells before filtering:", ncol(seurat_obj), "\n")
cat("Cells after filtering:", ncol(filtered_obj), "\n")

# Save the filtered Seurat object
saveRDS(
  filtered_obj,
  file = file.path(output_dir, paste0(sample_name, "_filtered.rds"))
)

cat("\nPreprocessing completed for", sample_name, "\n")
