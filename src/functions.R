# Process sample
process_data_simple <- function(dir_path) {
  data <- Read10X(data.dir = file.path(dir_path, "filtered_feature_bc_matrix/"))
  seurat_obj <- CreateSeuratObject(counts = data, min.cells = 2)
  seurat_obj$mitoPercent <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- subset(seurat_obj, subset = nCount_RNA > 2000 & nFeature_RNA > 1000 & mitoPercent < 10)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 5000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  #seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
  #seurat_obj <- FindClusters(seurat_obj, resolution = c(0.2, 0.4, 0.6, 0.8, 0.9, 1.2)) 
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
  return(seurat_obj)
}

process_data_qc <- function(dir_path) {
  data <- Read10X(data.dir = file.path(dir_path, "filtered_feature_bc_matrix/"))
  seurat_obj <- CreateSeuratObject(counts = data, min.cells = 2)
  seurat_obj$mitoPercent <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  return(seurat_obj)
}

# Calculate multiplet rate
get_mparam <- function(cells) {
  cell_values <- c(1000, 5000, 10000)
  mparam_values <- c(0.008, 0.04, 0.08)
  
  # Linear interpolation
  if (cells <= min(cell_values)) {
    return(min(mparam_values))
  } else if (cells >= max(cell_values)) {
    return(max(mparam_values))
  } else {
    approx(cell_values, mparam_values, xout = cells)$y
  }
}

# Find doublets
find_doublet <- function(sce){
  sweep.sce <- paramSweep(sce, PCs = 1:40, sct = F)
  sweep.sce.stat <- summarizeSweep(sweep.sce, GT = F)
  bcmvn <- find.pK(sweep.sce.stat)
  pK <- bcmvn %>% top_n(., 1, BCmetric) %>% 
    pull(pK) %>% as.character() %>% as.numeric()
  m.rate <- get_mparam(sce %>% ncol())
  annot <- sce@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annot)
  nExp_poi <- round(m.rate*nrow(sce@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  sce.db <- doubletFinder(sce, PCs = 1:40, pN=0.25, pK=pK, 
                          nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
  cells.use <- colnames(sce.db)[which(sce.db[[sce.db@meta.data %>% ncol()]] == "Singlet")]
  sce.db.sub <- subset(sce.db, cells = cells.use)
  return(sce.db.sub)
}

# Merge Seurat objects, normalize, transform
merge_seurat_list <- function(seurat_list) {
  seu.merged <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)], add.cell.ids = seurat_list %>% names())
  seu.merged <- subset(seu.merged, subset = nCount_RNA > 2000 & nFeature_RNA > 1000 & mitoPercent < 10) %>% 
    NormalizeData() %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>% 
    CellCycleScoring(g2m.features=cc.genes.updated.2019$g2m.genes, s.features=cc.genes.updated.2019$s.genes,) %>% 
    ScaleData(vars.to.regress = c("S.Score", "G2M.Score"))
  return(seu.merged)
}

# Run integration
run_integration <- function(seurat_merged) {
  mer_seu <- RunPCA(seurat_merged, npcs = 40)
  harmonized_seurat <- RunHarmony(mer_seu, assay = "SCT",
                                  group.by.vars = c("Sample"),
                                  reduction = "pca", reduction.save = "harmony",
                                  ncores = 20, max_iter = 60)
  harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony")
  harmonized_seurat <- FindClusters(harmonized_seurat, graph.name = "SCT_snn", resolution = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.6, 0.8, 1.0, 1.2))
  harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", dims = 1:12, n.epochs = 5000)
  return(harmonized_seurat)
}

# Create QC scatter plot
create_qc_plot <- function(qc_metrics, title) {
  col.pal <- colorRampPalette(brewer.pal(11, "RdYlGn"))(21)
  
  ggplot(qc_metrics, aes(nCount_RNA, nFeature_RNA, fill = mitoPercent)) + 
    geom_point(colour = "black", pch = 21, size = 1.5) + 
    scale_fill_gradientn(colours = rev(col.pal), limits = c(0, 100)) +
    ggtitle("") + xlab("Counts per cell") + ylab("Unique genes per cell") + labs(fill = '% mt') +
    scale_x_continuous(labels = comma, limits = c(0, 180000), breaks = seq(0, 320000, 80000), minor_breaks = seq(0, 320000, 20000)) +
    scale_y_continuous(labels = comma, limits = c(0, 10000)) + 
    theme_bw() + 
    ggtitle(title) +
    theme(legend.key.width = unit(0.2, "cm"), plot.margin = margin(0.2, 0.1, 0, 0.2, "cm"), 
          plot.title = element_text(size = 10, margin = margin(t = 20, b = -20), hjust = 0.1))
}

# Create QC violin plots
violin_plot <- function(df, ylim, x = "patient", y) {
  ylab <- switch(
    y,
    "nCount_RNA" = "Counts per cell",
    "mitoPercent" = "% mt",
    "nFeature_RNA" = "Unique genes per cell",
    stop("Invalid column name for y argument")
  )
  
  ggplot(df, aes_string(x = x, y = y, fill = x)) +
    geom_jitter(shape = 16, position = position_jitter(0.4), size = 0.2) + 
    geom_violin() + geom_boxplot(width = 0.08, color = "black", alpha = 0.6, outlier.shape = NA, fill = "white") + 
    scale_y_continuous(labels = scales::comma, limits = c(0, ylim)) +
    scale_fill_brewer(palette = "Set2") + theme_bw() + theme(legend.position = "none") +
    xlab("") + ylab(ylab)
}

# Proceed with DE genes

