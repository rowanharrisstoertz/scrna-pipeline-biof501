# cell type annotation script
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
    stop("Enter: input_file, output_dir, sample_name, marker_file, annotation_min_score, annotation_min_difference")
}

input_file <- args[1]
output_dir <- args[2]
sample_name <- args[3]
marker_file <- args[4]
annotation_min_score <- as.numeric(args[5])
annotation_min_difference <- as.numeric(args[6])

# ensure annotation thresholds are valid
if (anyNA(c(annotation_min_score, annotation_min_difference))) {
    stop("Annotation thresholds must be numeric")
}
if (annotation_min_score < 0) {
    stop("annotation_min_score must be non-negative")
}
if (annotation_min_difference < 0) {
    stop("annotation_min_difference must be non-negative")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load the filtered Seurat object
seurat_obj <- readRDS(input_file)

DefaultAssay(seurat_obj) <- "SCT"
Idents(seurat_obj) <- "seurat_clusters"

# load the sample markers CSV file
markers <- read.csv(
    marker_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

cat("\nLoaded", nrow(markers), "markers from", marker_file, "\n")

# check for correct columns in the markers file
if (!all(c("cell_type", "gene") %in% colnames(markers))) {
    stop("Marker file must contain 'cell_type' and 'gene' columns")
}

# remove any other columns from the markers file, keeping only 'cell_type' and 'gene'
markers <- markers %>% 
  select(cell_type, gene)

# remove whitespace from cell_type and gene columns
markers <- markers %>%
  mutate(cell_type = trimws(cell_type),
         gene = trimws(gene))

# remove rows with NA or empty values in cell_type or gene columns
markers <- markers %>%
  filter(!is.na(cell_type) & cell_type != "" & !is.na(gene) & gene != "")

markers <- unique(markers)

cat("\nLoaded", nrow(markers), "unique markers\n")

# get list of available genes in the Seurat object
available_genes <- rownames(seurat_obj[["SCT"]])

# add a new column to the markers data frame indicating whether each gene is present in the Seurat object
present_markers <- markers %>%
    filter(gene %in% available_genes)
    
if (nrow(present_markers) == 0) {
    stop("No marker genes present in dataset")
}

#group genes into marker panels
marker_panels <- split(present_markers$gene, present_markers$cell_type)

cat("\nCalculating average expression of marker genes in each cluster...\n")

# calculate average expression of each marker gene in each cluster
average_expression <- AverageExpression(seurat_obj, features = unique(present_markers$gene), assays = "SCT", group.by = "seurat_clusters")$SCT

cat("\nScoring marker genes for each cluster...\n")

# Score each cell type in each cluster
marker_scores <- lapply(names(marker_panels), function(cell_type) {
    genes <- intersect(marker_panels[[cell_type]], rownames(average_expression))

    if (length(genes) == 0) {
        return(NULL)
        }

    data.frame(cluster = colnames(average_expression), cell_type = cell_type, score = colMeans(average_expression[genes, , drop = FALSE]))
})
marker_scores <- do.call(rbind, marker_scores)


# Create a cell type score table
marker_score_table <- marker_scores %>%
    pivot_wider(
        names_from = cell_type,
        values_from = score
    ) %>%
    arrange(match(cluster, colnames(average_expression)))

# Save annotation scores
write.csv(
    marker_score_table,
    file = file.path(output_dir, paste0(sample_name, "_annotation_scores.csv")),
    row.names = FALSE
)

cat("\nAnnotation score table saved\n")


cat("\nAssigning cell types to clusters based on scores...\n")

# Assign cell types to clusters based on scores
annotation_summary <- marker_scores %>%
  group_by(cluster) %>%
  arrange(desc(score), .by_group = TRUE) %>%
  summarise(
    best_score = first(score),
    second_score = nth(score, 2, default = 0),
    best_type = first(cell_type),
    cell_type = case_when(
      best_score < annotation_min_score ~ "Unassigned",
      best_score - second_score < annotation_min_difference ~ paste("Ambiguous — possible", best_type),
      TRUE ~ best_type
    ),
    .groups = "drop"
  ) %>%
  select(cluster, cell_type, best_score)


cat("\nCell type assignment completed\n")

# add cell types to the Seurat object metadata
annotation_map <- setNames(
    annotation_summary$cell_type,
    sub("^g", "", annotation_summary$cluster)
)

seurat_obj$cell_type <- unname(
    annotation_map[as.character(seurat_obj$seurat_clusters)]
)

cat("\nCell type counts:\n")
print(table(seurat_obj$cell_type, useNA = "always"))

#create annotated UMAP plot
umap_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) + ggtitle(paste(sample_name, "— cell annotations")) + theme_minimal()

# Save annotated UMAP plot
ggsave(file.path(output_dir, paste0(sample_name, "_annotated_umap.png")), plot = umap_plot, width = 10, height = 6, dpi = 300)

# Save the annotated Seurat object
saveRDS(
  seurat_obj,
  file = file.path(output_dir, paste0(sample_name, "_annotated.rds"))
)


cat("\nAnnotation completed for", sample_name, "\n")

