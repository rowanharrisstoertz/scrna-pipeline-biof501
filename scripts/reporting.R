library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
    stop("Enter: input_dir, markers.csv, output_dir, sample_name")
}

input_file <- args[1]
markers_file <- args[2]
output_dir <- args[3]
sample_name <- args[4]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# Load the filtered Seurat object and markers
seurat_obj <- readRDS(input_file)
markers <- read.csv(markers_file)

DefaultAssay(seurat_obj) <- "SCT"

cluster_counts <- as.data.frame(table(cluster = Idents(seurat_obj)))

# UMAP plot
umap_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle(paste(sample_name, "- Cell Clusters")) + theme_minimal()

# Save UMAP plot
ggsave(file.path(output_dir, paste0(sample_name, "_umap.png")), plot = umap_plot, width = 8, height = 6, dpi = 300)

# QC violin plots for nFeature_RNA, nCount_RNA, and percent.mt
qc_plot1 <- VlnPlot(seurat_obj, features = "nFeature_RNA") + ggtitle("Number of Unique Genes")+ theme_minimal()
qc_plot2 <- VlnPlot(seurat_obj, features = "nCount_RNA") + ggtitle("Total Gene Counts") + NoLegend()+ theme_minimal()
qc_plot3 <- VlnPlot(seurat_obj, features = "percent.mt") + ggtitle("Mitochondrial Percentage") + NoLegend()+ theme_minimal()

qc_plot <- qc_plot1 + qc_plot2 + qc_plot3 + plot_annotation(title = paste(sample_name, "- QC Metrics"))

# Save QC violin plots
ggsave(file.path(output_dir, paste0(sample_name, "_qc_violin.png")), plot = qc_plot, width = 10, height = 8, dpi = 300)

# marker gene heatmap
top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
    ungroup()

heatmap_plot <- DoHeatmap(seurat_obj, features = top_markers$gene) + ggtitle(paste(sample_name, "- Top Marker Genes Heatmap")) + theme_minimal() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())

# Save marker gene heatmap
ggsave(file.path(output_dir, paste0(sample_name, "_marker_heatmap.png")), plot = heatmap_plot, width = 10, height = 8, dpi = 300)


# summary report 
summary_file <- file.path(
    output_dir,
    paste0(sample_name, "_summary.txt")
)

sink(summary_file)

cat("\nSingle Cell Analysis Report for sample", sample_name, "\n")
cat("Number of cells after filtering:", ncol(seurat_obj), "\n")
cat("Total genes detected:", nrow(seurat_obj), "\n")
cat("Number of clusters:", length(unique(Idents(seurat_obj))), "\n")
cat("Number of marker genes identified:", nrow(markers), "\n")
cat("Cell counts per cluster:\n")
print(cluster_counts)
cat("Top marker genes per cluster:\n")
for (cluster_id in unique(top_markers$cluster)) {
    cat("\nCluster", cluster_id, ":\n")
    cluster_markers <- top_markers %>%
        filter(cluster == cluster_id) %>%
        select(gene, avg_log2FC, p_val_adj)
    print(as.data.frame(cluster_markers), row.names = FALSE)
}

sink()