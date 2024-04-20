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

# Source
source("/Users/mouse/GitRepo/MAB_script/src/functions.R")

# List directories
local.path <- "/Volumes/Extreme SSD/scRNA_MAB/Raws/"

directories <- list.dirs(local.path, full.names = FALSE, recursive = FALSE)
directories <- directories[grep("SASC", directories)] # Remove this

seurat.list <- list()

# Read data. filter out doublets, save

# Preprocess raw data
for (dir in directories) {
  sample.name <- toupper(dir)
  print(dir)  
  seurat.obj <- process_data_simple(file.path(local.path, dir, fsep = ""))
  seurat.obj <- AddMetaData(seurat.obj, metadata = dir, col.name = "Sample")
  seurat.list[[sample.name]] <- seurat.obj
}

result_list <- list()

# Remove doublets
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
  seurat.obj <- process_data_qc(file.path(local.path, dir, fsep = ""))
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

