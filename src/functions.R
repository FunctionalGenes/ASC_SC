# Process sample
read_and_mt <- function(dir_path) {
  data <- Read10X(data.dir = file.path(dir_path, "filtered_feature_bc_matrix/"))
  seurat_obj <- CreateSeuratObject(counts = data, min.cells = 2)
  seurat_obj$mitoPercent <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  return(seurat_obj)
}

filter_data <- function(seurat_obj) {
  seurat_obj <- subset(seurat_obj, subset = nCount_RNA > 2000 & nFeature_RNA > 1000 & mitoPercent < 10)
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
  harmonized_seurat <- RunHarmony(mer_seu,
                                  group.by.vars = c("Sample"),
                                  reduction = "pca", reduction.save = "harmony",
                                  ncores = 20, max_iter = 60)
  harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony")
  harmonized_seurat <- FindClusters(harmonized_seurat, graph.name = "RNA_snn", resolution = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.6, 0.8, 1.0, 1.2))
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

# Gene enrichment
rearrange_degs <- function(degs_table, all_names, log2fc_filter) {
  if (log2fc_filter == "up") {
    degs <- degs_table[degs_table$avg_log2FC > 0, ]
  } else if (log2fc_filter == "down") {
    degs <- degs_table[degs_table$avg_log2FC < 0, ]
  } else {
    stop()
  }
  degs <- tibble::rownames_to_column(degs, "Gene.Name") %>% 
    dplyr::select(Gene.Name, p_val_adj)
  missing_genes <- setdiff(all_names, degs$Gene.Name)
  sim_p <- runif(length(missing_genes), 0.06, 1)
  missing_genes_df <- data.frame(Gene.Name = missing_genes,
                                 p_val_adj = sim_p)
  rearranged_degs <- rbind(degs, missing_genes_df)
  colnames(rearranged_degs)[2] <- "padj"
  return(rearranged_degs)
}

topgo_ready <- function(df, onto = "BP", db = "org.Hs.eg.db"){ 
  df.vector <- df[,2]
  names(df.vector) <- df[,1]
  df.vector <- df.vector[!is.na(df.vector)]
  selection <- function(allScore){ return(allScore < 0.05)}
  allGO2genes <- annFUN.org(whichOnto=onto, feasibleGenes=NULL, mapping=db, ID="symbol")
  GOdata <- new("topGOdata",
                ontology=onto,
                allGenes=df.vector,
                annot=annFUN.GO2genes,
                GO2genes=allGO2genes,
                geneSel=selection,
                nodeSize=8) 
  
  results.fisher <- runTest(GOdata, algorithm="weight01", statistic="fisher") 
  goEnrichment <- GenTable(GOdata, fisher=results.fisher, orderBy="fisher", topNodes=length(score(results.fisher)), numChar=1000)
  goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
  goEnrichment$fisher[is.na(goEnrichment$fisher)] <- as.numeric(1.1e-30)
  goEnrichment <- goEnrichment[goEnrichment$fisher<0.05,]
  annotGenes = lapply(goEnrichment$GO.ID, function(x) as.character(unlist(genesInTerm(object = GOdata, whichGO = x))))
  z <- as.data.frame(df) %>% filter(padj < 0.05)
  sigGenes = lapply(annotGenes, function(x) intersect(x, z[,1])) 
  goEnrichment <- goEnrichment[,c("GO.ID","Term","Annotated", "Significant", "fisher", "Expected")]
  goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
  temp.df <- goEnrichment
  temp.df$Identified_Genes <- vapply(sigGenes, paste, collapse = ", ", character(1L))
  temp.df <- temp.df %>% dplyr::mutate(FoldEnrichment = Significant / Expected)
  
  return(temp.df)
}

all_go <- function(df, db = "org.Hs.eg.db"){
  enrich.df <- data.frame()
  for(j in c("BP", "MF", "CC")){
    temp <- topgo_ready(df, onto = j, db)
    enrich.df <- rbind(enrich.df, data.frame(temp, Ontology = rep(j, dim(temp)[1])))
  }
  return(enrich.df)
}

generate_go_data <- function(markers, integr_names, log2fc_filter) {
  rearrange_degs(markers, integr_names, log2fc_filter = log2fc_filter) %>%
    all_go() %>%
    mutate(Activation = ifelse(log2fc_filter == "up", "positive", "negative"))
}

# Draw modified DimPlot
cluster_plot <- function(dim.plot){
  custom_cols <- c("firebrick", "gold", "steelblue", "forestgreen", "sienna",
              "orchid4", "orange", "palegreen", "steelblue1", "darkgrey", "pink", "blue")
  colnames(dim.plot$data)[3] <- "clust"
  ggplot(dim.plot$data, aes(x=UMAP_1, y = UMAP_2, group = clust)) + 
    geom_point(aes(fill=clust), shape=21, color="black", size=1.8, stroke=.13) + 
    scale_fill_manual(values=custom_cols) + theme_classic() + theme(legend.position = "none")
}
