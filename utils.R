#!/usr/bin/env Rscript

# Set a CRAN mirror
options(repos = "https://cloud.r-project.org")

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("yaml"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("EnsDb.Hsapiens.v86"))
suppressPackageStartupMessages(library("SingleR"))
suppressPackageStartupMessages(library("celldex"))
suppressPackageStartupMessages(library("Signac"))
suppressPackageStartupMessages(library("dplyr"))

# scRNA and scATAC data processing
preprocess <- function(name, config) {
  # Load 10X data
  path <- paste0(
    config$cellranger_output_path, "/outs/filtered_feature_bc_matrix.h5"
  )
  rawdata <- Read10X_h5(path)

  # extract RNA and ATAC data
  rna_counts <- rawdata$`Gene Expression`
  atac_counts <- rawdata$Peaks

  # Create Seurat object
  data <- CreateSeuratObject(counts = rna_counts)
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

  # Add ATAC-seq data
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- config$atac_annot_source
  genome(annotations) <- config$atac_annot_genome

  frag.file <- paste0(
    config$cellranger_output_path, "/outs/atac_fragments.tsv.gz"
  ) 
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = config$atac_annot_genome,
    fragments = frag.file,
    min.cells = 10,
    annotation = annotations
  )
  data[["ATAC"]] <- chrom_assay

  # Filterting
  data <- subset(
    data,
    subset = nCount_ATAC >= config$min_genes_atac & 
      nCount_ATAC <= config$max_genes_atac &
      nCount_RNA >= config$min_genes_rna & 
      nCount_RNA <= config$max_genes_rna &
      percent.mt < config$max_mt_percentage
  )

  return(data)
}


# Cluster and normalization
cluster <- function(data, config, output) {

  # RNA Cluster
  DefaultAssay(data) <- "RNA"
  gene_names <- rownames(data)
  data <- NormalizeData(data) %>%
    FindVariableFeatures(nfeatures=config$variable_features) %>%
    ScaleData(features = gene_names) %>%
    RunPCA(npcs=config$npcs) %>%
    RunUMAP(
      dims=1:config$npcs,
      reduction.name='umap_rna',
      reduction.key='rnaUMAP_'
      )
  
  # Annotation using singleR and celldexr
  ref <- HumanPrimaryCellAtlasData()
  test <- GetAssayData(data, assay='RNA', slot='counts')
  pred <- SingleR(
    test=test,
    ref=ref,
    labels=ref$label.main
  )
  idx = match(rownames(data@meta.data), rownames(pred))
  data$singleR.labels <- pred$labels[idx]

  # ATAC cluster
  # We exclude the first dimension as this is typically correlated with sequencing depth
  DefaultAssay(data) <- "ATAC"
  data <- RunTFIDF(data) %>%
    FindTopFeatures(min.cutoff = 'q0') %>%
    RunSVD() %>%
    RunUMAP(
      reduction = 'lsi', 
      dims = 2:config$npcs, 
      reduction.name = "umap_atac", 
      reduction.key = "atacUMAP_"
    )

  # Merged RNA and ATAC using WNN
  data <- FindMultiModalNeighbors(
      data, 
      reduction.list = list("pca", "lsi"), 
      dims.list = list(1:config$npcs, 2:config$npcs)
    ) %>%
    RunUMAP(
      nn.name = "weighted.nn", 
      reduction.name = "wnn_umap", 
      reduction.key = "wnnUMAP_"
    ) %>%
    FindClusters(
      graph.name = "wsnn", 
      algorithm = 3, 
    verbose = FALSE)
    
  # Draw figures
  # Draw RNA clusters
  png(paste0(output, "/rna_clusters_umap.png"), width=800, height=800)
  g <- DimPlot(
    data, 
    reduction="umap_rna",
    group.by='singleR.labels',
    label=T,
    repel=T
  )
  print(g)
  dev.off()

  # Draw score heatmap
  png(paste0(output, "/rna_score_heatmap_data.png"), width=600, height=600)
  g <- plotScoreHeatmap(pred)
  print(g)
  dev.off()

  # Draw delta distribution
  png(paste0(output, "/rna_delta_distribution.png"), width=600, height=600)
  g <- plotDeltaDistribution(pred)
  print(g)
  dev.off()

  # Draw ATAC clusters
  png(paste0(output, "/atac_clusters_umap.png"), width=800, height=800)
  g <- DimPlot(
    data, 
    reduction="umap_atac",
    group.by='singleR.labels',
    label=T,
    repel=T
  )
  print(g)
  dev.off()

  # Draw WNN clusters
  png(paste0(output, "/merged_clusters_umap.png"), width=800, height=800)
  g <- DimPlot(
    data, 
    reduction="wnn_umap",
    group.by='singleR.labels',
    label=T,
    repel=T
  )
  print(g)
  dev.off()

  return(data)
}
