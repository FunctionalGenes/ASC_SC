# Process sample
process_data_simple <- function(dir_path) {
  data <- Read10X(data.dir = file.path(dir_path, "filtered_feature_bc_matrix/"))
  seurat_obj <- CreateSeuratObject(counts = data, min.cells = 2)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 5000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
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
  cells.use <- colnames(sce.db)[which(sce.db[[5]] == "Singlet")]
  sce.db.sub <- subset(sce.db, cells = cells.use)
  return(sce.db.sub)
}
