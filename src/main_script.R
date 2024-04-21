library(Seurat)
library(harmony)
library(DoubletFinder)
library(dplyr)
library(DropletUtils)
library(purrr)
library(tidyverse)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(cowplot)
library(topGO)

### Add info about packages versions used

# Source
source("functions.R")

# List directories
local.path <- "/ASC_SC/MEX/"

directories <- list.dirs(local.path, full.names = FALSE, recursive = FALSE)

# Preprocess raw data
seurat.list <- list()
for (dir in directories) {
  sample.name <- toupper(dir)
  print(dir)  
  seurat.obj <- read_and_mt(file.path(local.path, dir, fsep = ""))
  seurat.obj <- filter_data(seurat.obj)
  seurat.obj <- AddMetaData(seurat.obj, metadata = dir, col.name = "Sample")
  seurat.list[[sample.name]] <- seurat.obj
}

# Remove doublets
result_list <- list()
for (i in seq_along(seurat.list)) {
  result_list[[names(seurat.list)[i]]] <- find_doublet(seurat.list[[i]])
}

# Save results of doublet removal
for (i in seq_along(result.list)) {
  seurat_name <- names(result.list)[i]
  print(seurat_name)
  count_table <- result.list[[i]]@assays$RNA$counts
  directory_path <- file.path("/Volumes/Extreme SSD/scRNA_MAB/TEST", seurat_name)
  write10xCounts(directory_path, count_table)
}

## Create QC plots
seurat.list.qc <- list()
for (dir in directories) {
  sample.name <- toupper(dir)
  print(dir)  
  seurat.obj <- read_and_mt(file.path(local.path, dir, fsep = ""))
  seurat.obj <- AddMetaData(seurat.obj, metadata = dir, col.name = "Sample")
  seurat.list.qc[[sample.name]] <- seurat.obj
}

# Fibroblasts
fib_seu <- seurat.list.qc[grep("FIB", names(seurat.list.qc))]
fib_names <- fib_seu %>% names()

fib_qc_metrics <- map2(fib_seu, fib_names, ~ {
  as_tibble(.x@meta.data[c("nCount_RNA", "nFeature_RNA", "mitoPercent")], rownames = "Cell.Barcode") %>%
    mutate(patient = .y)
})

fib_qc <- bind_rows(fib_qc_metrics)
fib_qc_plots <- map2(fib_qc_metrics, fib_names, ~ create_qc_plot(.x, .y))
fib_qc_panel <- plot_grid(fib_qc_plots[[1]], fib_qc_plots[[2]], fib_qc_plots[[3]], ncol = 2)

fib.counts <- violin_plot(fib_qc, 100000, y = "nCount_RNA")
fib.mito <- violin_plot(fib_qc, 40, y = "mitoPercent")
fib.genes <- violin_plot(fib_qc, 10000, y = "nFeature_RNA")

fib_vln_panel <- plot_grid(fib.counts, fib.genes, fib.mito, ncol = 3, nrow = 1, rel_widths = c(1, 1, 1))
fib_panel <- plot_grid(fib_qc_panel, NULL, fib_vln_panel, ncol = 1, rel_heights = c(2, 0.1, 1), 
                         labels = c("A", "", "B"), label_fontfamily = "sans", label_size = 16)

# SASC
sasc_seu <- seurat.list.qc[grep("SASC", names(seurat.list.qc))]
sasc_names <- sasc_seu %>% names()

sasc_qc_metrics <- map2(sasc_seu, sasc_names, ~ {
  as_tibble(.x@meta.data[c("nCount_RNA", "nFeature_RNA", "mitoPercent")], rownames = "Cell.Barcode") %>%
    mutate(patient = .y)
})

sasc_qc <- bind_rows(sasc_qc_metrics)
sasc_qc_plots <- map2(sasc_qc_metrics, sasc_names, ~ create_qc_plot(.x, .y))
sasc_qc_panel <- plot_grid(sasc_qc_plots[[1]], sasc_qc_plots[[2]], sasc_qc_plots[[3]], 
                           sasc_qc_plots[[4]], ncol = 2)

sasc.counts <- violin_plot(sasc_qc, 100000, y = "nCount_RNA")
sasc.mito <- violin_plot(sasc_qc, 40, y = "mitoPercent")
sasc.genes <- violin_plot(sasc_qc, 10000, y = "nFeature_RNA")

sasc_vln_panel <- plot_grid(sasc.counts, sasc.genes, sasc.mito, ncol = 3, nrow = 1, rel_widths = c(1, 1, 1))
sasc_panel <- plot_grid(sasc_qc_panel, NULL, sasc_vln_panel, ncol = 1, rel_heights = c(2, 0.1, 1), 
                       labels = c("A", "", "B"), label_fontfamily = "sans", label_size = 16)

# VASC
vasc_seu <- seurat.list.qc[grep("VASC", names(seurat.list.qc))]
vasc_names <- vasc_seu %>% names()

vasc_qc_metrics <- map2(vasc_seu, vasc_names, ~ {
  as_tibble(.x@meta.data[c("nCount_RNA", "nFeature_RNA", "mitoPercent")], rownames = "Cell.Barcode") %>%
    mutate(patient = .y)
})

vasc_qc <- bind_rows(vasc_qc_metrics)
vasc_qc_plots <- map2(vasc_qc_metrics, vasc_names, ~ create_qc_plot(.x, .y))
vasc_qc_panel <- plot_grid(vasc_qc_plots[[1]], vasc_qc_plots[[2]], vasc_qc_plots[[3]], 
                           vasc_qc_plots[[4]], ncol = 2)

vasc.counts <- violin_plot(vasc_qc, 100000, y = "nCount_RNA")
vasc.mito <- violin_plot(vasc_qc, 40, y = "mitoPercent")
vasc.genes <- violin_plot(vasc_qc, 10000, y = "nFeature_RNA")

vasc_vln_panel <- plot_grid(vasc.counts, vasc.genes, vasc.mito, ncol = 3, nrow = 1, rel_widths = c(1, 1, 1))
vasc_panel <- plot_grid(vasc_qc_panel, NULL, vasc_vln_panel, ncol = 1, rel_heights = c(2, 0.1, 1), 
                        labels = c("A", "", "B"), label_fontfamily = "sans", label_size = 16)

# Integrate Fibroblast samples
fib.seurat.list.nodoub <- list()
fib.directories <- directories[grep("FIB", directories)]

for (dir in fib.directories) {
  sample.name <- toupper(dir)
  print(dir)  
  seurat.obj <- read_and_mt(file.path(local.path, dir, fsep = ""))
  seurat.obj <- AddMetaData(seurat.obj, metadata = dir, col.name = "Sample")
  fib.seurat.list.nodoub[[sample.name]] <- seurat.obj
}

fib.integr <- merge_seurat_list(fib.seurat.list.nodoub) %>% run_integration()
cluster_plot(DimPlot(fib.integr, reduction = "umap", group.by = "RNA_snn_res.0.4", label = T))

# Integrate SASC samples
sasc.seurat.list.nodoub <- list()
sasc.directories <- directories[grep("SASC", directories)]

for (dir in sasc.directories) {
  sample.name <- toupper(dir)
  print(dir)  
  seurat.obj <- read_and_mt(file.path(local.path, dir, fsep = ""))
  seurat.obj <- AddMetaData(seurat.obj, metadata = dir, col.name = "Sample")
  sasc.seurat.list.nodoub[[sample.name]] <- seurat.obj
}

sasc.integr <- merge_seurat_list(sasc.seurat.list.nodoub) %>% run_integration()
cluster_plot(DimPlot(sasc.integr, reduction = "umap", group.by = "RNA_snn_res.0.4", label = T))

# Integrate VASC samples
vasc.seurat.list.nodoub <- list()
vasc.directories <- directories[grep("VASC", directories)]

for (dir in vasc.directories) {
  sample.name <- toupper(dir)
  print(dir)  
  seurat.obj <- read_and_mt(file.path(local.path, dir, fsep = ""))
  seurat.obj <- AddMetaData(seurat.obj, metadata = dir, col.name = "Sample")
  vasc.seurat.list.nodoub[[sample.name]] <- seurat.obj
}

vasc.integr <- merge_seurat_list(vasc.seurat.list.nodoub) %>% run_integration()
cluster_plot(DimPlot(vasc.integr, reduction = "umap", group.by = "RNA_snn_res.0.4", label = T))

# Integrate all samples
all.seurat.list.nodoub <- list()

for (dir in directories) {
  sample.name <- toupper(dir)
  print(dir)  
  seurat.obj <- read_and_mt(file.path(local.path, dir, fsep = ""))
  seurat.obj <- AddMetaData(seurat.obj, metadata = dir, col.name = "Sample")
  all.seurat.list.nodoub[[sample.name]] <- seurat.obj
}

all.integr <- merge_seurat_list(all.seurat.list.nodoub) %>% run_integration()
cluster_plot(DimPlot(all.integr, reduction = "umap", group.by = "RNA_snn_res.0.4", label = T))

# Differential expression analysis
all.integr@meta.data$Type <- ifelse(grepl("^FIB", all.integr@meta.data$Sample), "FIBR",
                                           ifelse(grepl("^SASC", all.integr@meta.data$Sample), "SASC",
                                                  ifelse(grepl("^VASC", all.integr@meta.data$Sample), "VASC", NA)))

markers.fib.sasc <- all.integr %>% FindMarkers(ident.1 = "FIBR", ident.2 = "SASC", verbose = FALSE, group.by = "Type")
markers.fib.vasc <- all.integr %>% FindMarkers(ident.1 = "FIBR", ident.2 = "VASC", verbose = FALSE, group.by = "Type")
markers.vasc.sasc <- all.integr %>% FindMarkers(ident.1 = "VASC", ident.2 = "SASC", verbose = FALSE, group.by = "Type")

# Gene enrichment
fib.sasc.go <- bind_rows(
  generate_go_data(markers.fib.sasc, rownames(all.integr), "up"),
  generate_go_data(markers.fib.sasc, rownames(all.integr), "down")
)

fib.vasc.go <- bind_rows(
  generate_go_data(markers.fib.vasc, rownames(all.integr), "up"),
  generate_go_data(markers.fib.vasc, rownames(all.integr), "down")
)

vasc.sasc.go <- bind_rows(
  generate_go_data(markers.vasc.sasc, rownames(all.integr), "up"),
  generate_go_data(markers.vasc.sasc, rownames(all.integr), "down")
)
