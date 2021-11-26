# scRNA-seq analysis of larval CNS across time
# Creating figures components
# Author: Noah R Dillon 
# Date 11252021 

# Set directory 
setwd("Workspace/") 

# Load libraries 
{
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(glmGamPoi)
  library(ggplot2)
  library(sctransform)
  library(cowplot)
  library(plyr)
}

# Creating the whole atlas from Cell Ranger output
{
  # Load in data sets filtered matrix files 
  {
    atlas_1h_alh <- Read10X("/Users/agg_1h_alh/filtered_feature_bc_matrix")
    atlas_24h_alh <- Read10X("/Users/agg_24h_alh/filtered_feature_bc_matrix")
    atlas_48h_alh <- Read10X("/Users/agg_48h_alh/filtered_feature_bc_matrix")
    atlas_96h_alh <- Read10X("/Users/agg_96h_alh/filtered_feature_bc_matrix")
  }
  # Create Seurat objects for each time point
  {
    atlas_1h_alh <- CreateSeuratObject(counts = atlas_1h_alh, project = "atlas_1h_alh", min.cells = 3, min.features = 100)
    atlas_24h_alh <- CreateSeuratObject(counts = atlas_24h_alh, project = "atlas_24h_alh", min.cells = 3, min.features = 100)
    atlas_48h_alh <- CreateSeuratObject(counts = atlas_48h_alh, project = "atlas_48h_alh", min.cells = 3, min.features = 100)
    atlas_96h_alh <- CreateSeuratObject(counts = atlas_96h_alh, project = "atlas_96h_alh", min.cells = 3, min.features = 100)
  }
  # Filtering out high mitochondrial and low count cells 
  {
    atlas_1h_alh[["percent.mt"]] <- PercentageFeatureSet(atlas_1h_alh, pattern = "^mt:")
    atlas_24h_alh[["percent.mt"]] <- PercentageFeatureSet(atlas_24h_alh, pattern = "^mt:")
    atlas_48h_alh[["percent.mt"]] <- PercentageFeatureSet(atlas_48h_alh, pattern = "^mt:")
    atlas_1h_alh <- subset(atlas_1h_alh, subset = nFeature_RNA > 200 & percent.mt < 20)
    atlas_24h_alh <- subset(atlas_24h_alh, subset = nFeature_RNA > 200 & percent.mt < 20)
    atlas_48h_alh <- subset(atlas_48h_alh, subset = nFeature_RNA > 200 & percent.mt < 20)
    atlas_96h_alh <- subset(atlas_96h_alh, subset = nFeature_RNA > 200 & percent.mt < 20)
  }
  # Preparing for integration 
  {
    # Create list for integration
    data.list <- c(atlas_1h_alh, atlas_24h_alh, atlas_48h_alh, atlas_96h_alh)
    # normalize and identify variable features for each data set independently
    data.list <- lapply(X = data.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
    # select features that are repeatedly variable across data sets for integration
    features <- SelectIntegrationFeatures(object.list = data.list)
    # Pre run PCA to reduce computational time 
    data.list <- lapply(X = data.list, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = TRUE)
      x <- RunPCA(x, features = features, verbose = TRUE)
    })
  }
  # Integrated analysis 
  {
    atlas.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features, reduction = "rpca", dims = 1:50)
    atlas.combined <- IntegrateData(anchorset = atlas.anchors)
    
    DefaultAssay(atlas.combined) <- "integrated"
    # Run the standard work flow for visualization and clustering
    atlas.combined <- ScaleData(atlas.combined, verbose = TRUE)
    atlas.combined <- RunPCA(atlas.combined, npcs = 50, verbose = TRUE)
    atlas.combined <- RunUMAP(atlas.combined, reduction = "pca", dims = 1:50)
    atlas.combined <- FindNeighbors(atlas.combined, reduction = "pca", dims = 1:50)
    atlas.combined <- FindClusters(atlas.combined, resolution = 2.2) 
    DefaultAssay(atlas.combined) <- "RNA"
    
    # save progress
    saveRDS(atlas.combined, file = "atlas_combined.rds")
    # read RDS file
    atlas.combined <- readRDS("atlas_combined.rds")
  }
}

# Figure 1 analysis 
{
  # Whole atlas 
  {
    # Create UMAP (used in figure)
    atlas_UMAP <- DimPlot(atlas.combined, reduction = "umap", label = FALSE) + NoLegend() + NoAxes()
    atlas_UMAP
    ggsave(filename = "Figure_1/atlas_UMAP.tiff", plot = atlas_UMAP, device = "tiff",
           scale = 1, width = 35, height = 35, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    
    # Find cluster markers to identify cell identities 
    atlas.all.markers <- FindAllMarkers(object = atlas.combined, assay = "RNA", 
                                        logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                        min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
    write.csv(atlas.all.markers, file = "atlas_all_markers.csv")
    
    # Run feature plots 
    major.clusters.list <- c("dpn","Hey","nSyb","repo")
    atlas_clustering_features <- FeaturePlot(atlas.combined, features = major.clusters.list, min.cutoff = "q5", order = TRUE, 
                cols = c("light grey", "dark red"), ncol = 2,label = FALSE, combine = TRUE)
    atlas_clustering_features
    ggsave(filename = "Figure_1/atlas_clustering_features.tiff", plot = atlas_clustering_features, device = "tiff",
           scale = 1, width = 40, height = 40, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # Progenitors 
  {
    # subset re-cluster and label progenitors 
    {
      # subset 
      progenitor.list <- c(7,11,32,33,42,6,14,1,40,41,20,24,16,17,47,15)
      progenitors_unlabeled <- atlas.combined
      progenitors_unlabeled <- subset(atlas.combined, idents = progenitor.list)

      # re-cluster 
      DefaultAssay(progenitors_unlabeled) <- "integrated"
      progenitors_unlabeled <- RunUMAP(progenitors_unlabeled, reduction = "pca", dims = 1:50)
      progenitors_unlabeled <- FindNeighbors(progenitors_unlabeled, reduction = "pca", dims = 1:50)
      progenitors_unlabeled <- FindClusters(progenitors_unlabeled, resolution = 0.49) 
      DefaultAssay(progenitors_unlabeled) <- "RNA"
      
      # visualize re-clustering (plot used in figrue)
      Progenitors_UMAP <- DimPlot(progenitors_unlabeled, reduction = "umap", label = TRUE) + NoLegend() + NoAxes()
      ggsave(filename = "Figure_1/Progenitors_UMAP.tiff", plot = Progenitors_UMAP, device = "tiff",
             scale = 1, width = 25, height = 25, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
      
      # save progress
      saveRDS(progenitors_unlabeled, file = "progenitors_unlabeled.rds")
      # read RDS file
      progenitors_unlabeled <- readRDS("progenitors_unlabeled.rds")
      
      # label 
      progenitors_labeled <- progenitors_unlabeled
      new.cluster.ids <- c("Immature neurons - 0","GMCs - 1","Type II neuroblasts - 2", "INPs - 3", 
                           "Immature neurons - 4", "Immature neurons - 5", "Immature neurons - 6",
                           "Low quality - 7", "New-born neurons - 8", "Type I neuroblasts - 9",
                           "Immature neurons - 10", "Immature neurons - 11", "Quiescent neuroblasts - 12")
      names(new.cluster.ids) <- levels(progenitors_labeled)
      progenitors_labeled <- RenameIdents(progenitors_labeled, new.cluster.ids)
      # Visualize
      DimPlot(progenitors_labeled, reduction = "umap", label = TRUE, label.size = 6, label.color = "black", label.box = TRUE) + 
        NoLegend() + NoAxes()
      
      # order 
      order.list <- c("Quiescent neuroblasts - 12","Type I neuroblasts - 9","Type II neuroblasts - 2","INPs - 3", "GMCs - 1", 
                      "New-born neurons - 8","Immature neurons - 4","Immature neurons - 0", "Immature neurons - 10", "Immature neurons - 6",
                      "Immature neurons - 11", "Immature neurons - 5","Low quality - 7")
      progenitors_labeled@active.ident <- factor(x = progenitors_labeled@active.ident, levels = order.list) # change the order of the factor levels
      # Sub setting out low quality cluster 
      progenitors_labeled_noLQ <- subset(progenitors_labeled, idents = c("Quiescent neuroblasts - 12","Type I neuroblasts - 9","Type II neuroblasts - 2","INPs - 3", "GMCs - 1", 
                                                                         "New-born neurons - 8","Immature neurons - 4","Immature neurons - 0", "Immature neurons - 10", "Immature neurons - 6",
                                                                         "Immature neurons - 11", "Immature neurons - 5"))
      # visualize 
      DimPlot(progenitors_labeled_noLQ, reduction = "umap", label = TRUE, label.size = 6, label.color = "black", label.box = TRUE) + 
        NoLegend() + NoAxes()
      
      DefaultAssay(progenitors_labeled) <- "RNA"
      DefaultAssay(progenitors_labeled_noLQ) <- "RNA"
      
      # save
      saveRDS(progenitors_labeled, file = "progenitors_labeled.rds")
      saveRDS(progenitors_labeled_noLQ, file = "progenitors_labeled_noLQ.rds")
      # load
      progenitors_labeled <- readRDS("progenitors_labeled.rds")
      progenitors_labeled_noLQ <- progenitors_unlabeled <- readRDS("progenitors_labeled_noLQ.rds")
    }
    
    # Markers (used in subsequent analyses)
    {
      # Finding all markers for clusters of interest 
      progenitors.all.markers <- FindAllMarkers(object = progenitors_unlabeled, assay = "RNA", logfc.threshold = 0.1, test.use = "wilcox",
                                             min.pct = 0.1, min.diff.pct = 0.1, verbose = TRUE, return.thresh = 0.0501)
      write.csv(progenitors.all.markers, file = "progenitors.all.markers.csv")
    }
    
    # Dot plot 
    {
      Progenitor.identities.list <- c("trbl","dpn","mira","insc","CycE","wor","ase","pnt","tll","erm","ham","Phs",
                                      "tap","dap","Hey","E(spl)m6-BFM","fne","CadN","elav","nSyb","brp","Gad1","ChAT","VGlut")
      # plot used in figure 
      Progenitors_idents_plot <- DotPlot(progenitors_labeled, features = Progenitor.identities.list, cols = c("grey", "red"), 
              col.min = 0, cluster.idents = FALSE, dot.min = 0, scale.by = "size", dot.scale = 6) + RotatedAxis() 
      Progenitors_idents_plot
      ggsave(filename = "Figure_1/progenitor_idents_dotplot.tiff", plot = Progenitors_idents_plot, device = "tiff",
             scale = 1, width = 40, height = 11, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    }
    
    # Progenitor UMAP by time 
    {
      # temporal progenitor UMAP
      NPCs_temporal <- progenitors_unlabeled
      # sub-setting by time 
      NPCs_temporal_1h <-subset(x = NPCs_temporal, subset = orig.ident == "atlas_1h_alh" ) 
      NPCs_temporal_24h <-subset(x = NPCs_temporal, subset = orig.ident == "atlas_24h_alh") 
      NPCs_temporal_48h <-subset(x = NPCs_temporal, subset = orig.ident == "atlas_48h_alh") 
      # plotting UMAPs by time point
      plot_NPCs_temporal_1h <- DimPlot(NPCs_temporal_1h, reduction = "umap", label = FALSE) + NoLegend() + NoAxes()
      plot_NPCs_temporal_24h <- DimPlot(NPCs_temporal_24h, reduction = "umap", label = FALSE) + NoLegend() + NoAxes()
      plot_NPCs_temporal_48h <- DimPlot(NPCs_temporal_48h, reduction = "umap", label = FALSE) +  NoLegend() + NoAxes()
      # Plots used in figure 
      plot_NPCs_temporal_1h
      ggsave(filename = "Figure_1/plot_NPCs_temporal_1h.tiff", plot = plot_NPCs_temporal_1h, device = "tiff",
             scale = 1, width = 40, height = 40, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
      plot_NPCs_temporal_24h 
      ggsave(filename = "Figure_1/plot_NPCs_temporal_24h.tiff", plot = plot_NPCs_temporal_24h, device = "tiff",
             scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
      plot_NPCs_temporal_48h 
      ggsave(filename = "Figure_1/plot_NPCs_temporal_48h.tiff", plot = plot_NPCs_temporal_48h, device = "tiff",
             scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    }
  }
}

# Figure 2 analysis - qNBs
{
  # UMAPs 
  {
    # used in figure 
    cluster_12_UMAP <- DimPlot(progenitors_unlabeled, reduction = "umap", label = TRUE, cols = c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey",
                                                                              "grey","grey","#FF68A1"), label.size = 6, label.color = "black") + NoLegend() + NoAxes()
    ggsave(filename = "Figure_2/cluster_12_UMAP.tiff", plot = cluster_12_UMAP, device = "tiff",
           scale = 1, width = 25, height = 25, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    # sub-setting by time 
    cluster_12_1h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_1h_alh" ) 
    cluster_12_24h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_24h_alh") 
    cluster_12_48h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_48h_alh") 
    # plotting UMAPs by time point
    plot_cluster_12_1h <- DimPlot(cluster_12_1h, reduction = "umap", label = FALSE,  cols = c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey",
                                                                                              "grey","grey","#FF68A1")) + NoLegend() + NoAxes()
    plot_cluster_12_24h <- DimPlot(cluster_12_24h, reduction = "umap", label = FALSE, cols = c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey",
                                                                                               "grey","grey","#FF68A1")) + NoLegend() + NoAxes()
    plot_cluster_12_48h <- DimPlot(cluster_12_48h, reduction = "umap", label = FALSE,  cols = c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey",
                                                                                                "grey","grey","#FF68A1")) + NoLegend() + NoAxes()
    # used in figure 
    plot_cluster_12_1h
    ggsave(filename = "Figure_2/plot_cluster_12_1h.tiff", plot = plot_cluster_12_1h, device = "tiff",
           scale = 1, width = 25, height = 25, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    plot_cluster_12_24h 
    ggsave(filename = "Figure_2/plot_cluster_12_24h.tiff", plot = plot_cluster_12_24h, device = "tiff",
           scale = 1, width = 25, height = 25, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    plot_cluster_12_48h 
    ggsave(filename = "Figure_2/plot_cluster_12_48h.tiff", plot = plot_cluster_12_48h, device = "tiff",
           scale = 1, width = 25, height = 25, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # dot plots 
  {
    # dot plot top 10 and validated genes (used in figure)
    cluster12.genes.list <- c("Thor", "E23","CG1572", "CG44325","lncRNA:CR43334","lncRNA:CR31044","CG17646","Ccdc85","Chd64","lin-28",
                              "trbl","dpn","mira","insc","CycE","wor","ase","pnt","tll")
    cluster12_dotplot <- DotPlot(progenitors_labeled_noLQ, features = cluster12.genes.list, cols = c("grey", "red"), 
            cluster.idents = FALSE, col.min = 0, dot.min = 0.0, scale.by = "size", assay = "RNA") + RotatedAxis() 
    ggsave(filename = "Figure_2/cluster12_dotplot.tiff", plot = cluster12_dotplot, device = "tiff",
           scale = 1, width = 25, height = 11, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    
    # dot plot of model (used in figrue)
    siganling.markers.list <- c("InR","DENR","poly","lin-28",
                                "Akt1","Cul1","slmb","Roc1a","SkpA","eff","trbl",
                                "Tor","foxo","Thor", "REPTOR-BP", "REPTOR",
                                "stg", "CycE","S6k","Alk")
    signaling_plot <- DotPlot(progenitors_labeled_noLQ, features = siganling.markers.list, cols = c("grey", "red"), 
            cluster.idents = FALSE, col.min = 0, dot.min = 0.0, scale.by = "size", assay = "RNA") + RotatedAxis()
    ggsave(filename = "Figure_2/signaling_plot.tiff", plot = signaling_plot, device = "tiff",
           scale = 1, width = 23, height = 11, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # Model was informed from cluster 12 (qNB) expression found in the progenitors 
  # all markers .csv file 
  
  # Glial signaling to qNBs 
  {
    # repo positive (pan glial)
    DotPlot(atlas.combined, features = "repo", cols = c("grey", "red"), col.min = 0, cluster.idents = FALSE, dot.min = 0, scale.by = "size", dot.scale = 6) + RotatedAxis() 
    glial_unlabeled <- subset(atlas.combined, idents = c('9','19','25','26','29','36','37','44'))
    DimPlot(glial_unlabeled, reduction = "umap", label = TRUE, label.size = 6, label.color = "black") + NoLegend() + NoAxes()
    
    # UMAPing 
    DefaultAssay(glial_unlabeled) <- "integrated"
    glial_unlabeled <- RunUMAP(glial_unlabeled, reduction = "pca", dims = 1:50)
    glial_unlabeled <- FindNeighbors(glial_unlabeled, reduction = "pca", dims = 1:50)
    glial_unlabeled <- FindClusters(glial_unlabeled, resolution = 0.045)
    # used in figure 
    glial_UMAP <- DimPlot(glial_unlabeled, reduction = "umap", label = TRUE, label.size = 6, label.color = "black") + NoLegend() + NoAxes()
    ggsave(filename = "Figure_2/glial_UMAP.tiff", plot = glial_UMAP, device = "tiff",
           scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    
    DefaultAssay(glial_unlabeled) <- "RNA"
    
    # save
    saveRDS(glial_unlabeled, file = "glial_unlabeled.rds")
    # load 
    glial_unlabeled <- readRDS("glial_unlabeled.rds")
    
    glial.markers.list <- c("wrapper","hoe1",
                            "moody","AdamTS-A",
                            "CG6126","Indy","CG4797", 
                            "Gat","alrm")
    DotPlot(glial_unlabeled, features = glial.markers.list, cols = c("grey", "red"), col.min = 0, cluster.idents = FALSE, dot.min = 0, scale.by = "size", dot.scale = 6) + RotatedAxis() 
    
    # Markers for glial 
    # Finding all markers for clusters of interest (used in figure)
    glial.all.markers <- FindAllMarkers(object = glial_unlabeled, assay = "RNA", logfc.threshold = 0.1, test.use = "wilcox",
                                        min.pct = 0.1, min.diff.pct = 0.1, verbose = TRUE, return.thresh = 0.0501)
    write.csv(glial.all.markers, file = "glial_all_markers.csv")
    
    # label and order 
    glial_labeled <- glial_unlabeled
    DefaultAssay(glial_labeled) <- "RNA"
    new.cluster.ids <- c("Cortex/chiasm glia - 0", "Perineurial glia - 1", "Astrocytes/neuropil glia - 2", "Subperineurial glia - 3")
    names(new.cluster.ids) <- levels(glial_labeled)
    glial_labeled <- RenameIdents(glial_labeled, new.cluster.ids)
    order.list <- c("Cortex/chiasm glia - 0","Subperineurial glia - 3", "Perineurial glia - 1", "Astrocytes/neuropil glia - 2")
    glial_labeled@active.ident <- factor(x = glial_labeled@active.ident, levels = order.list) 
    # used in figure 
    glial_idents_plot <- DotPlot(glial_labeled, features = glial.markers.list, cols = c("grey", "red"), col.min = 0, cluster.idents = FALSE, dot.min = 0, scale.by = "size", dot.scale = 6) + RotatedAxis()
    glial_idents_plot
    ggsave(filename = "Figure_2/glial_idents_plot.tiff", plot = glial_idents_plot, device = "tiff",
           scale = 1, width = 17, height = 11, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    
    # save 
    saveRDS(glial_labeled, file = "glial_labeled.rds")
    # load 
    glial_labeled <- readRDS("glial_labeled.rds")
    
    # exclude 96h in temporal analysis (too few cells from samples)
    glial_labeled_temp_subset <- subset(glial_labeled, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))
    
    # save 
    saveRDS(glial_labeled_temp_subset, file = "glial_labeled_temp_subset.rds")
    # load 
    glial_labeled_temp_subset <- readRDS("glial_labeled_temp_subset.rds")
    
    signaling.glial.markers.list <- c("Ilp2","Ilp3","Ilp5","Ilp6","Ilp7","Ilp8","ana","trol")
    # used in figure 
    glial_signaling <- DotPlot(glial_labeled_temp_subset, features = signaling.glial.markers.list, cols = c("red","red","red","red"), col.min = 0, 
                               cluster.idents = FALSE, dot.min = 0, assay = "RNA", split.by = "orig.ident", scale.by = "size") + RotatedAxis()
    glial_signaling
    ggsave(filename = "Figure_2/glial_signaling.tiff", plot = glial_signaling, device = "tiff",
           scale = 1, width = 18, height = 10, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # Glial temporal component
  {
    Cortex_chiasm_glial <- subset(glial_labeled, idents = "Cortex/chiasm glial - 0")
    Subneurial_glial <- subset(glial_labeled, idents = "Subperineurial glial - 3")
    Perineurial_glial <- subset(glial_labeled, idents = "Perineurial glial - 1")
      
    # subset of time 
    {
        Cortex_chiasm_glial_temp <- subset(Cortex_chiasm_glial, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))

        Subperineurial_glial_temp <- subset(Subneurial_glial, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))

        Perineurial_glial_temp <- subset(Perineurial_glial, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))
    }
      
    # tables (used in figure)
    {
        Idents(Cortex_chiasm_glial_temp) <- "orig.ident"
        Cortex_chiasm_glial.temporal.markers <- FindAllMarkers(object = Cortex_chiasm_glial_temp, assay = "RNA", logfc.threshold = 0.1, test.use = "wilcox",
                                                               min.pct = 0.1, min.diff.pct = 0.1, verbose = TRUE, return.thresh = 0.0501)
        write.csv(Cortex_chiasm_glial.temporal.markers, file = "Cortex_chiasm_glial_temporal_markers.csv")
        
        Idents(Subperineurial_glial_temp) <- "orig.ident"
        Subperineurial_glial.temporal.markers <- FindAllMarkers(object = Subperineurial_glial_temp, assay = "RNA", logfc.threshold = 0.1, test.use = "wilcox",
                                                            min.pct = 0.1, min.diff.pct = 0.1, verbose = TRUE, return.thresh = 0.0501)
        write.csv(Subperineurial_glial.temporal.markers, file = "Subperineurial_glial_temporal_markers.csv")
        
        Idents(Perineurial_glial_temp) <- "orig.ident"
        Perineurial_glial.temporal.markers <- FindAllMarkers(object = Perineurial_glial_temp, assay = "RNA", logfc.threshold = 0.1, test.use = "wilcox",
                                                             min.pct = 0.1, min.diff.pct = 0.1, verbose = TRUE, return.thresh = 0.0501)
        write.csv(Perineurial_glial.temporal.markers, file = "Perineurial_glial_temporal_markers.csv")
    }
  }
}

# Figure 3 analysis - T1NBs
{
  # UMAPs
  {
    # used in figure 
    cluster_9_UMAP <- DimPlot(progenitors_unlabeled, reduction = "umap", label = TRUE, cols = c("grey","grey","grey","grey","grey","grey","grey","grey","grey","#8494FF","grey",
                                                                                       "grey","grey"), label.size = 6, label.color = "black") + NoLegend() + NoAxes()
    ggsave(filename = "Figure_3/cluster_9_UMAP.tiff", plot = cluster_9_UMAP, device = "tiff",
           scale = 1, width = 25, height = 25, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    # sub-setting by time 
    cluster_9_1h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_1h_alh" ) 
    cluster_9_24h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_24h_alh") 
    cluster_9_48h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_48h_alh") 
    # plotting UMAPs by time point
    plot_cluster_9_1h <- DimPlot(cluster_9_1h, reduction = "umap", label = FALSE, cols = c("grey","grey","grey","grey","grey","grey","grey","grey","grey","#8494FF","grey",
                                                                                           "grey","grey")) + NoLegend() + NoAxes()
    plot_cluster_9_24h <- DimPlot(cluster_9_24h, reduction = "umap", label = FALSE, cols = c("grey","grey","grey","grey","grey","grey","grey","grey","grey","#8494FF","grey",
                                                                                             "grey","grey")) + NoLegend() + NoAxes()
    plot_cluster_9_48h <- DimPlot(cluster_9_48h, reduction = "umap", label = FALSE, cols = c("grey","grey","grey","grey","grey","grey","grey","grey","grey","#8494FF","grey",
                                                                                             "grey","grey")) + NoLegend() + NoAxes()
    # used in figure 
    plot_cluster_9_1h
    ggsave(filename = "Figure_3/plot_cluster_9_1h.tiff", plot = plot_cluster_9_1h, device = "tiff",
           scale = 1, width = 30, height = 30, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    plot_cluster_9_24h
    ggsave(filename = "Figure_3/plot_cluster_9_24h.tiff", plot = plot_cluster_9_24h, device = "tiff",
           scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    plot_cluster_9_48h 
    ggsave(filename = "Figure_3/plot_cluster_9_48h.tiff", plot = plot_cluster_9_48h, device = "tiff",
           scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # dot plots
  {
    # dot plot top 10 and validated genes 
    cluster9.genes.list <- c("lncRNA:CR31386","grh","Syp","lncRNA:CR33938","Pen","lncRNA:CR46003","lncRNA:CR30009","CG13305",
                             "CycE","stg","wor","ase","dpn","mira","insc")
    # used in figure 
    cluster_9_dotplot <- DotPlot(progenitors_labeled_noLQ, features = cluster9.genes.list, cols = c("grey", "red"), 
            cluster.idents = FALSE, col.min = 0, dot.min = 0.0, scale.by = "size", assay = "RNA") + RotatedAxis() 
    ggsave(filename = "Figure_3/cluster_9_dotplot.tiff", plot = cluster_9_dotplot, device = "tiff",
           scale = 1, width = 25, height = 11, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    
    # temporal markers 
    cluster_9_subset <- subset(cluster_9, subset = (orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))
    Idents(cluster_9_subset) <- "orig.ident"
    cluster9_temporal.markers <- FindAllMarkers(object = cluster_9_subset, assay = "RNA", 
                                                logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                                min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
    write.csv(cluster9_temporal.markers, file = "cluster9_temporal_markers.csv")
    
    # TFs identified as Panther's GO analysis gene-specific transcriptional regulator protein class
    cluster9_sub_temp.TF.list <- c("SoxN","fru","Dp","wge","BtbVII","klu","bi","tara","ken","Crtc",
                                   "NC2alpha","luna","Pdp1","Mnt","mamo","Hmx","disco-r","ab","ttk",
                                   "cas","lin-28","svp","chinmo","Imp","EcR","br","Eip93F","Syp")
    #used in figrue 
    cluster_9_temp <- DotPlot(cluster_9_subset, features = cluster9_sub_temp.TF.list, cols = c("grey","red"), col.min = 0, 
            cluster.idents = FALSE, dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis() 
    cluster_9_temp
    ggsave(filename = "Figure_3/cluster_9_temp.tiff", plot = cluster_9_temp, device = "tiff",
           scale = 1, width = 25, height = 11, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
}

# Figure 4 analysis - T2NBs
{
  # Cluster 2 UMAPs
  {
    # UMAP (used in figure)
    cluster_2_UMAP <- DimPlot(progenitors_unlabeled, reduction = "umap", label = TRUE, cols = c("grey","grey","#CD9600","grey","grey","grey","grey","grey","grey","grey",
                                                                              "grey","grey","grey"), label.size = 6, label.color = "black") + NoLegend() + NoAxes()
    ggsave(filename = "Figure_4/cluster_2_UMAP.tiff", plot = cluster_4_UMAP, device = "tiff",
           scale = 1, width = 25, height = 25, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    # sub-setting by time 
    cluster_2_1h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_1h_alh" ) 
    cluster_2_24h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_24h_alh") 
    cluster_2_48h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_48h_alh") 
    # plotting UMAPs by time point
    plot_cluster_2_1h <- DimPlot(cluster_2_1h, reduction = "umap", label = FALSE, cols = c("grey","grey","#CD9600","grey","grey","grey","grey","grey","grey","grey",
                                                                                           "grey","grey","grey")) + NoLegend() + NoAxes()
    plot_cluster_2_24h <- DimPlot(cluster_2_24h, reduction = "umap", label = FALSE, cols = c("grey","grey","#CD9600","grey","grey","grey","grey","grey","grey","grey",
                                                                                             "grey","grey","grey")) + NoLegend() + NoAxes()
    plot_cluster_2_48h <- DimPlot(cluster_2_48h, reduction = "umap", label = FALSE, cols = c("grey","grey","#CD9600","grey","grey","grey","grey","grey","grey","grey",
                                                                                             "grey","grey","grey")) + NoLegend() + NoAxes()
    # used in figure 
    plot_cluster_2_1h
    ggsave(filename = "Figure_4/plot_cluster_2_1h.tiff", plot = plot_cluster_2_1h, device = "tiff",
           scale = 1, width = 25, height = 25, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    plot_cluster_2_24h 
    ggsave(filename = "Figure_4/plot_cluster_2_24h.tiff", plot = plot_cluster_2_24h, device = "tiff",
           scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    plot_cluster_2_48h 
    ggsave(filename = "Figure_4/plot_cluster_2_48h.tiff", plot = plot_cluster_2_48h, device = "tiff",
           scale = 1, width = 22, height = 22, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # Sub cluster analysis of cluster 2 
  {
    cluster2_sub_clustered <- subset(progenitors_unlabeled, idents = 2)
    
    DefaultAssay(cluster2_sub_clustered) <- "integrated"
    cluster2_sub_clustered <- RunUMAP(cluster2_sub_clustered, reduction = "pca", dims = 1:50)
    cluster2_sub_clustered <- FindNeighbors(cluster2_sub_clustered, reduction = "pca", dims = 1:50)
    cluster2_sub_clustered <- FindClusters(cluster2_sub_clustered, resolution = 0.1)
    DefaultAssay(cluster2_sub_clustered) <- "RNA"
    # visualize 
    DimPlot(cluster2_sub_clustered, reduction = "umap", label = TRUE) + NoLegend() + NoAxes()
    
    # save
    saveRDS(cluster2_sub_clustered, file = "cluster2_sub_clustered.rds")
    # load 
    cluster2_sub_clustered <- readRDS("cluster2_sub_clustered.rds")
    
    # T2NB specific markers  
    cluster2_sub_dotplot <- DotPlot(cluster2_sub_clustered, features = c("dpn","pnt","tll","ase", "CycE"), cols = c("grey", "red"), col.min = 0, cluster.idents = FALSE, dot.min = 0, scale.by = "size", dot.scale = 6) + RotatedAxis() 
    cluster2_sub_dotplot
    # In supplemental figure 
    ggsave(filename = "Figure_4/cluster2_sub_dotplot.tiff", plot = cluster2_sub_dotplot, device = "tiff",
           scale = 1, width = 15, height = 15, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    
    # markers 
    T2P.all.markers <- FindAllMarkers(object = cluster2_sub_clustered, assay = "RNA", logfc.threshold = 0.1, test.use = "wilcox",
                                        min.pct = 0.1, min.diff.pct = 0.1, verbose = TRUE, return.thresh = 0.0501)
    write.csv(T2P.all.markers, file = "T2P_all_markers.csv")
    
    # label 
    new.cluster.ids <- c("Nonannotated progenitors", "Type II neuroblasts")
    names(new.cluster.ids) <- levels(cluster2_sub_clustered)
    cluster2_sub_clustered <- RenameIdents(cluster2_sub_clustered, new.cluster.ids)
    
    # UMAP (used in figrue)
    cluster_2_sub_UMAP <- DimPlot(cluster2_sub_clustered, reduction = "umap", label = FALSE, label.size = 6, label.color = "Black") + 
      NoLegend() + NoAxes()
    ggsave(filename = "Figure_4/cluster_2_sub_UMAP.tiff", plot = cluster_2_sub_UMAP, device = "tiff",
           scale = 1, width = 15, height = 15, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    
    # sub-setting by time 
    cluster_2_subsetb_1h <-subset(x = cluster2_sub_clustered, subset = orig.ident == "atlas_1h_alh" ) 
    cluster_2_subsetb_24h <-subset(x = cluster2_sub_clustered, subset = orig.ident == "atlas_24h_alh") 
    cluster_2_subsetb_48 <-subset(x = cluster2_sub_clustered, subset = orig.ident == "atlas_48h_alh") 
    # plotting UMAPs by time point
    plot_cluster_2_1h_b <- DimPlot(cluster_2_subsetb_1h, reduction = "umap", label = FALSE, cols = c("grey","#00BFC4")) + NoLegend() + NoAxes()
    plot_cluster_2_24h_b <- DimPlot(cluster_2_subsetb_24h, reduction = "umap", label = FALSE, cols = c("grey","#00BFC4")) + NoLegend() + NoAxes()
    plot_cluster_2_48h_b <- DimPlot(cluster_2_subsetb_48, reduction = "umap", label = FALSE, cols = c("grey","#00BFC4")) + NoLegend() + NoAxes()
    # used in figure 
    plot_cluster_2_1h_b
    ggsave(filename = "Figure_4/plot_cluster_2_1h_b.tiff", plot = plot_cluster_2_1h_b, device = "tiff",
           scale = 1, width = 16, height = 16, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    plot_cluster_2_24h_b
    ggsave(filename = "Figure_4/plot_cluster_2_24h_b.tiff", plot = plot_cluster_2_24h_b, device = "tiff",
           scale = 1, width = 15, height = 15, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    plot_cluster_2_48h_b
    ggsave(filename = "Figure_4/plot_cluster_2_48h_b.tiff", plot = plot_cluster_2_48h_b, device = "tiff",
           scale = 1, width = 15, height = 15, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    
    # subset for T2NBs
    cluster2_sub_cluster_1b <- subset(cluster2_sub_clustered, idents = "Type II neuroblasts")
    # subset for 24h and 48h 
    cluster2_sub_cluster_1b <- subset(cluster2_sub_cluster_1b, subset = (orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))
    # markers 
    T2NBs.all.markers <- FindAllMarkers(object = cluster2_sub_cluster_1b, assay = "RNA", logfc.threshold = 0.1, test.use = "wilcox",
                                        min.pct = 0.1, min.diff.pct = 0.1, verbose = TRUE, return.thresh = 0.0501)
    write.csv(T2NBs.all.markers, file = "T2NBs_temp_markers.csv")
  }
  
  # Dot plots 
  {
    # dot plot top 10 and validated genes 
    cluster2.genes.list <- c("hid","ft","E(spl)mgamma-HLH","bowl","fru","CG4250", "E(spl)m4-BFM","jbug","Tis11","Brd",
                             "pnt","tll","ase")
    # used in figure 
    cluster_2_dotplot <- DotPlot(progenitors_labeled_noLQ, features = cluster2.genes.list, cols = c("grey", "red"), 
            cluster.idents = FALSE, col.min = 0, dot.min = 0.0, scale.by = "size", assay = "RNA") + RotatedAxis() 
    ggsave(filename = "Figure_4/cluster_2_dotplot.tiff", plot = cluster_2_dotplot, device = "tiff",
           scale = 1, width = 23, height = 11, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    
    # TFs identified as Panther's GO analysis gene-specific transcriptional regulator protein class
    cluster2_sub_temp.TF.list <- c("MTA1-like","BtbVII","fru","SoxN","kra",
                                   "cas","lin-28","svp","chinmo","Imp","Syp","br","EcR","Eip93F")
    # used in figure                
    cluster2_temp <- DotPlot(cluster2_sub_cluster_1b, features = cluster2_sub_temp.TF.list, cols = c("grey","red"), col.min = 0, cluster.idents = FALSE, 
            dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis()
    cluster2_temp
    ggsave(filename = "Figure_4/cluster2_temp.tiff", plot = cluster2_temp, device = "tiff",
           scale = 1, width = 20, height = 11, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
}

# Figure 5 analysis - INPs 
{
  # UMAPs
  {
    # used in figure 
    cluster_3_UMAP <- DimPlot(progenitors_unlabeled, reduction = "umap", label = TRUE, cols = c("grey","grey","grey","#ABA300","grey","grey","grey","grey","grey","grey",
                                                                                       "grey","grey","grey"), label.size = 6, label.color = "black") + NoLegend() + NoAxes()
    ggsave(filename = "Figure_5/cluster_3_UMAP.tiff", plot = cluster_3_UMAP, device = "tiff",
           scale = 1, width = 25, height = 25, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    # sub set out cluster 3
    cluster_3 <-subset(x = progenitors_unlabeled, idents = '3')
    # sub-setting by time 
    cluster_3_1h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_1h_alh" ) 
    cluster_3_24h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_24h_alh") 
    cluster_3_48h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_48h_alh") 
    # plotting UMAPs by time point
    plot_cluster_3_1h <- DimPlot(cluster_3_1h, reduction = "umap", label = FALSE, cols = c("grey","grey","grey","#ABA300","grey","grey","grey","grey","grey","grey",
                                                                                           "grey","grey","grey")) + NoLegend() + NoAxes()
    plot_cluster_3_24h <- DimPlot(cluster_3_24h, reduction = "umap", label = FALSE, cols = c("grey","grey","grey","#ABA300","grey","grey","grey","grey","grey","grey",
                                                                                             "grey","grey","grey")) + NoLegend() + NoAxes()
    plot_cluster_3_48h <- DimPlot(cluster_3_48h, reduction = "umap", label = FALSE, cols = c("grey","grey","grey","#ABA300","grey","grey","grey","grey","grey","grey",
                                                                                             "grey","grey","grey")) + NoLegend() + NoAxes()
    # used in figure 
    plot_cluster_3_1h
    ggsave(filename = "Figure_5/plot_cluster_3_1h.tiff", plot = plot_cluster_3_1h, device = "tiff",
           scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    plot_cluster_3_24h 
    ggsave(filename = "Figure_5/plot_cluster_3_24h.tiff", plot = plot_cluster_3_24h, device = "tiff",
           scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    plot_cluster_3_48h 
    ggsave(filename = "Figure_5/plot_cluster_3_48h.tiff", plot = plot_cluster_3_48h, device = "tiff",
           scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # Dot plots 
  {
    # dot plot top 10 and validated genes 
    cluster3.genes.list <- c("CycE","Mcm7","klu","trn","CG12708","hdly","CG32354","Pen","stg","lncRNA:CR30009",
                             "erm","ham")
    # used in figure 
    cluster_3_dotplot <- DotPlot(progenitors_labeled_noLQ, features = cluster3.genes.list, cols = c("grey", "red"), 
            cluster.idents = FALSE, col.min = 0, dot.min = 0.0, scale.by = "size", assay = "RNA") + RotatedAxis() 
    ggsave(filename = "Figure_5/cluster_3_dotplot.tiff", plot = cluster_3_dotplot, device = "tiff",
           scale = 1, width = 23, height = 11, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    
    # Find deferentially expressed genes by cell type 
    {
      # INPs - TypeII neuroblasts 
      cluster3_cluster2.markers <- FindMarkers(progenitors_unlabeled, ident.1 = 3, ident.2 = 2, min.pct = 0.1, logfc.threshold = 0.1)
      write.csv(cluster3_cluster2.markers, file = "cluster3_vs_cluster2_markers.csv")
      # INPs - GMCs
      cluster3_cluster1.markers <- FindMarkers(progenitors_unlabeled, ident.1 = 3, ident.2 = 1, min.pct = 0.1, logfc.threshold = 0.1)
      write.csv(cluster3_cluster1.markers, file = "cluster3_vs_cluster1_markers.csv")
      # INPs - T1NBs
      cluster3_cluster9.markers <- FindMarkers(progenitors_unlabeled, ident.1 = 3, ident.2 = 9, min.pct = 0.1, logfc.threshold = 0.1)
      write.csv(cluster3_cluster9.markers, file = "cluster3_vs_cluster9_markers.csv")
      
      # subseting clusters of interest 
      T2P_INP_GMC <- subset(progenitors_labeled, idents = c('Type II neuroblasts - 2','INPs - 3','GMCs - 1'))
      T1NB_INP <- subset(progenitors_labeled, idents = c("Type I neuroblasts - 9",'INPs - 3'))
      
      # INPs vs T2P
      INP.vs.T2P.up.list <- c("lncRNA:roX1","His3.3B","lncRNA:CR30009","Hsp27","lncRNA:CR45388","lncRNA:CR34335","CycE","klu","CG5359","SIFa")
      DotPlot(T2P_INP_GMC, features = INP.vs.T2P.up.list, cols = c("grey","red"), col.min = 0, cluster.idents = FALSE, dot.min = 0, assay = "RNA", scale.by = "size") + RotatedAxis() 
      INP.vs.T2P.down.list <- c("bowl","fru","ft","aop","dally","Tis11","RpL41","fog","Brd","heph")
      DotPlot(T2P_INP_GMC, features = INP.vs.T2P.down.list, cols = c("grey","red"), col.min = 0, cluster.idents = FALSE, dot.min = 0, assay = "RNA", scale.by = "size") + RotatedAxis() 

      # GMCs vs INP
      INP.vs.GMC.up.list <- c("Syp","trn","mamo","CycE","Hsp27","mnd","His3.3B","br","hdly","path")
      DotPlot(T2P_INP_GMC, features = INP.vs.GMC.up.list, cols = c("grey","red"), col.min = 0, cluster.idents = FALSE, dot.min = 0, assay = "RNA", scale.by = "size") + RotatedAxis() 
      INP.vs.GMC.down.list <- c("lncRNA:cherub","lncRNA:CR34335","E(spl)m8-HLH","grh","lncRNA:CR45388","lncRNA:CR45593","Lmpt","Nrg","lncRNA:CR44272","Bacc")
      DotPlot(T2P_INP_GMC, features = INP.vs.GMC.down.list, cols = c("grey","red"), col.min = 0, cluster.idents = FALSE, dot.min = 0, assay = "RNA", scale.by = "size") + RotatedAxis() 

      # T1NB vs INP
      INP.vs.GMC.up.list <- c("mamo","trn","His3.3B","fax","CG1124","RpL39","CG1552","CG1648","CG43658","aay")
      DotPlot(T1NB_INP, features = INP.vs.GMC.up.list, cols = c("grey","red"), col.min = 0, cluster.idents = FALSE, dot.min = 0, assay = "RNA", scale.by = "size") + RotatedAxis() 
      INP.vs.GMC.down.list <- c("lncRNA:cherub","lncRNA:CR31386","grh","Pex7","fru","wdp","lncRNA:CR33938","svp","Syp","vvl")
      DotPlot(T1NB_INP, features = INP.vs.GMC.down.list, cols = c("grey","red"), col.min = 0, cluster.idents = FALSE, dot.min = 0, assay = "RNA", scale.by = "size") + RotatedAxis()
      
      # INP differentiating dot plots 
      T2NBs.INPs.GMCs.markers.list <- c("fog","Brd","bowl","ft","aop","heph","path","fru","dally","Tis11","Syp","RpL41",
                                        "hdly","br","mnd","mamo","His3.3B","Hsp27","trn","CycE","klu","CG5359","SIFa",
                                        "lncRNA:CR30009","lncRNA:roX1",
                                        "lncRNA:CR34335","lncRNA:cherub","lncRNA:CR45388","Bacc","Nrg","grh","lncRNA:CR44272","lncRNA:CR45593","E(spl)m8-HLH","Lmpt")
      # used in figrue INP vs GMC and T2NBs 
      T2P_INP_GMC_dotplot <- DotPlot(T2P_INP_GMC, features = T2NBs.INPs.GMCs.markers.list, cols = c("grey","red"), col.min = 0, cluster.idents = FALSE, 
                                     dot.min = 0, assay = "RNA", scale.by = "size") + RotatedAxis()
      T2P_INP_GMC_dotplot
      ggsave(filename = "Figure_5/T2P_INP_GMC_dotplot.tiff", plot = T2P_INP_GMC_dotplot, device = "tiff",
             scale = 1, width = 29, height = 11, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
      
      T1NBs.INPs.markers.list <- c("lncRNA:cherub","lncRNA:CR31386","grh","Pex7","fru","wdp","lncRNA:CR33938","svp","Syp","vvl",
                                   "mamo","trn","His3.3B","fax","RpL39","CG43658","CG1124","CG1552","CG1648","aay")
      # used in figure INP vs T1NBs
      T1NB_INP_dotplot <- DotPlot(T1NB_INP, features = T1NBs.INPs.markers.list, cols = c("grey","red"), col.min = 0, cluster.idents = FALSE, 
                                  dot.min = 0, assay = "RNA", scale.by = "size") + RotatedAxis()
      T1NB_INP_dotplot
      ggsave(filename = "Figure_5/T1NB_INP_dotplot.tiff", plot = T1NB_INP_dotplot, device = "tiff",
             scale = 1, width = 23, height = 11, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    }
  }
}

# Figure 6 analysis - GMCs + New-born neurons + Immature neurons 
{
  # UMAPs - GMCs 
  {
    # used in figure
    cluster_1_UMAP <- DimPlot(progenitors_unlabeled, reduction = "umap", label = TRUE, cols = c("grey","#E68613","grey","grey","grey","grey","grey","grey","grey","grey",
                                                                                       "grey","grey","grey","grey"), label.size = 6, label.color = "black") + NoLegend() + NoAxes()
    ggsave(filename = "Figure_6/cluster_1_UMAP.tiff", plot = cluster_1_UMAP, device = "tiff",
           scale = 1, width = 25, height = 25, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    # sub-setting by time 
    cluster_1_1h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_1h_alh" ) 
    cluster_1_24h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_24h_alh") 
    cluster_1_48h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_48h_alh") 
    # plotting UMAPs by time point
    plot_cluster_1_1h <- DimPlot(cluster_1_1h, reduction = "umap", label = FALSE, cols = c("grey","#E68613","grey","grey","grey","grey","grey","grey","grey","grey",
                                                                                           "grey","grey","grey","grey")) + NoLegend() + NoAxes()
    plot_cluster_1_24h <- DimPlot(cluster_1_24h, reduction = "umap", label = FALSE, cols = c("grey","#E68613","grey","grey","grey","grey","grey","grey","grey","grey",
                                                                                             "grey","grey","grey","grey")) + NoLegend() + NoAxes()
    plot_cluster_1_48h <- DimPlot(cluster_1_48h, reduction = "umap", label = FALSE, cols = c("grey","#E68613","grey","grey","grey","grey","grey","grey","grey","grey",
                                                                                             "grey","grey","grey","grey")) + NoLegend() + NoAxes()
    # used in figure
    plot_cluster_1_1h
    ggsave(filename = "Figure_6/plot_cluster_1_1h.tiff", plot = plot_cluster_1_1h, device = "tiff",
           scale = 1, width = 22, height = 22, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    plot_cluster_1_24h
    ggsave(filename = "Figure_6/plot_cluster_1_24h.tiff", plot = plot_cluster_1_24h, device = "tiff",
           scale = 1, width = 15, height = 15, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    plot_cluster_1_48h 
    ggsave(filename = "Figure_6/plot_cluster_1_48h.tiff", plot = plot_cluster_1_48h, device = "tiff",
           scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # UMAPs - New-borns
  {
    # used in figure 
    cluster8_UMAP <- DimPlot(progenitors_unlabeled, reduction = "umap", label = TRUE, cols = c("grey","grey","grey","grey","grey","grey","grey","grey","#00A9FF","grey","grey",
                                                                                               "grey","grey"), label.size = 6, label.color = "black") + NoLegend() + NoAxes()
    ggsave(filename = "Figure_7/cluster8_UMAP.tiff", plot = cluster8_UMAP, device = "tiff",
           scale = 1, width = 25, height = 25, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    # sub-setting by time 
    cluster_8_1h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_1h_alh" ) 
    cluster_8_24h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_24h_alh") 
    cluster_8_48h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_48h_alh") 
    # plotting UMAPs by time point
    plot_cluster_8_1h <- DimPlot(cluster_8_1h, reduction = "umap", label = FALSE, cols = c("grey","grey","grey","grey","grey","grey","grey","grey","#00A9FF","grey","grey",
                                                                                           "grey","grey")) + NoLegend() + NoAxes()
    plot_cluster_8_24h <- DimPlot(cluster_8_24h, reduction = "umap", label = FALSE, cols = c("grey","grey","grey","grey","grey","grey","grey","grey","#00A9FF","grey","grey",
                                                                                             "grey","grey")) + NoLegend() + NoAxes()
    plot_cluster_8_48h <- DimPlot(cluster_8_48h, reduction = "umap", label = FALSE, cols = c("grey","grey","grey","grey","grey","grey","grey","grey","#00A9FF","grey","grey",
                                                                                             "grey","grey")) + NoLegend() + NoAxes()
    # used in figure 
    plot_cluster_8_1h
    ggsave(filename = "Figure_6/plot_cluster_8_1h.tiff", plot = plot_cluster_8_1h, device = "tiff",
           scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    plot_cluster_8_24h 
    ggsave(filename = "Figure_6/plot_cluster_8_24h.tiff", plot = plot_cluster_8_24h, device = "tiff",
           scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    plot_cluster_8_48h 
    ggsave(filename = "Figure_6/plot_cluster_8_48h.tiff", plot = plot_cluster_8_48h, device = "tiff",
           scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # UMAPs - Immature neurons 
  {
    # used in figure 
    cluster_imn_UMAP <- DimPlot(progenitors_unlabeled, reduction = "umap", label = TRUE, cols = c("#F8766D","grey","grey","grey","#0CB702","#00C19A","#00BFC4","grey","grey","grey",
                                                                                                  "#ED68ED","#FF61CC","grey"), label.size = 6, label.color = "black") + NoLegend() + NoAxes()
    ggsave(filename = "Figure_8/cluster_imn_UMAP.tiff", plot = cluster_imn_UMAP, device = "tiff",
           scale = 1, width = 25, height = 25, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    # sub-setting by time 
    cluster_imn_1h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_1h_alh") 
    cluster_imn_24h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_24h_alh") 
    cluster_imn_48h <-subset(x = progenitors_unlabeled, subset = orig.ident == "atlas_48h_alh") 
    # plotting UMAPs by time point
    plot_cluster_imn_1h <- DimPlot(cluster_imn_1h, reduction = "umap", label = FALSE, cols = c("#F8766D","grey","grey","grey","#0CB702","#00C19A","#00BFC4","grey","grey","grey",
                                                                                               "#ED68ED","#FF61CC","grey")) + NoLegend() + NoAxes()
    plot_cluster_imn_24h <- DimPlot(cluster_imn_24h, reduction = "umap", label = FALSE, cols = c("#F8766D","grey","grey","grey","#0CB702","#00C19A","#00BFC4","grey","grey","grey",
                                                                                                 "#ED68ED","#FF61CC","grey")) + NoLegend() + NoAxes()
    plot_cluster_imn_48h <- DimPlot(cluster_imn_48h, reduction = "umap", label = FALSE, cols = c("#F8766D","grey","grey","grey","#0CB702","#00C19A","#00BFC4","grey","grey","grey",
                                                                                                 "#ED68ED","#FF61CC","grey")) + NoLegend() + NoAxes()
    # used in figure 
    plot_cluster_imn_1h
    ggsave(filename = "Figure_6/plot_cluster_imn_1h.tiff", plot = plot_cluster_imn_1h, device = "tiff",
           scale = 1, width = 30, height = 30, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    plot_cluster_imn_24h 
    ggsave(filename = "Figure_6/plot_cluster_imn_24h.tiff", plot = plot_cluster_imn_24h, device = "tiff",
           scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    plot_cluster_imn_48h
    ggsave(filename = "Figure_6/plot_cluster_imn_48h.tiff", plot = plot_cluster_imn_48h, device = "tiff",
           scale = 1, width = 20, height = 20, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # Dot plots - GMCs 
  {
    # dot plot top 10 and validated genes 
    cluster1.genes.list <- c("lncRNA:CR45388","lncRNA:cherub","grh","sprt","Galphai","E(spl)m8-HLH","cas","edl","stg",
                             "tap","dap","Phs")
    # used in figure 
    cluster1_dotplot <- DotPlot(progenitors_labeled_noLQ, features = cluster1.genes.list, cols = c("grey", "red"), 
            cluster.idents = FALSE, col.min = 0, dot.min = 0.0, scale.by = "size", assay = "RNA") + RotatedAxis()
    cluster1_dotplot
    ggsave(filename = "Figure_6/cluster1_dotplot.tiff", plot = cluster1_dotplot, device = "tiff",
           scale = 1, width = 23, height = 11, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # Dot plots - New-borns
  {
    # dot plot top 10 and validated genes 
    cluster8.genes.list <- c("lncRNA:CR44272","robo2","lncRNA:cherub","l(3)neo38","Liprin-gamma","CG12605", "Antp", "Nrt",
                             "Hey","E(spl)m6-BFM")
    # used in figure 
    cluster8_dotplot <- DotPlot(progenitors_labeled_noLQ, features = cluster8.genes.list, cols = c("grey", "red"), 
                                cluster.idents = FALSE, col.min = 0, dot.min = 0.0, scale.by = "size", assay = "RNA") + RotatedAxis() 
    ggsave(filename = "Figure_6/cluster8_dotplot.tiff", plot = cluster8_dotplot, device = "tiff",
           scale = 1, width = 23, height = 11, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # Dot plots - Immature neurons 
  {
    # dot plot top 5 and validated genes from each cluster
    clusterimn.genes.list <- c("Abd-B","Antp","stai","beat-Ia","Imp",
                               "Trim9","trv","mspo","Alk",
                               "fz2",
                               "ara","caup","Fas2","mirr",
                               "Ggamma30A","sn","Hsp23","miple1","Myo81F",
                               "CG12071","Fas3","CG4467","hwt","Reph",
                               "pros","jim","lncRNA:noe","pdm3","lncRNA:roX2","unk","hdc",
                               "fne","CadN","elav","nSyb","brp",
                               "Gad1","ChAT","VGlut")
    # fz2 shared between cluster 0 and 10 
    # used in figure 
    cluster_imn_dotplot <- DotPlot(progenitors_labeled_noLQ, features = clusterimn.genes.list, cols = c("grey", "red"), 
                                   cluster.idents = FALSE, col.min = 0, dot.min = 0.0, scale.by = "size", assay = "RNA") + RotatedAxis()
    cluster_imn_dotplot
    ggsave(filename = "Figure_6/cluster_imn_dotplot.tiff", plot = cluster_imn_dotplot, device = "tiff",
           scale = 1, width = 45, height = 11, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # Comparison between GMCs to T1NBs New-born neurons and immature neurons
  {
    # Find deferentially expressed genes 
    cluster1_cluster8.markers <- FindMarkers(progenitors_unlabeled, ident.1 = 1, ident.2 = 8, min.pct = 0.1, logfc.threshold = 0.1)
    write.csv(cluster1_cluster8.markers, file = "cluster1_cluster8_markers.csv")
    cluster1_cluster9.markers <- FindMarkers(progenitors_unlabeled, ident.1 = 1, ident.2 = 9, min.pct = 0.1, logfc.threshold = 0.1)
    write.csv(cluster1_cluster9.markers, file = "cluster1_cluster9_markers.csv")
    cluster8_cluster_imn.markers <- FindMarkers(progenitors_unlabeled, ident.1 = 8, ident.2 = c(4,0,10,6,11,5), min.pct = 0.1, logfc.threshold = 0.1)
    write.csv(cluster8_cluster_imn.markers, file = "cluster8_cluster_imn_markers.csv")
    
    # subseting clusters of interest 
    NB_GMC_NWB_imn <- subset(progenitors_labeled, idents = c('GMCs - 1','New-born neurons - 8','Type I neuroblasts - 9',"Immature neurons - 5","Immature neurons - 11", "Immature neurons - 10","Immature neurons - 6", "Immature neurons - 4","Immature neurons - 0"))
    # grouping all imns
    new.cluster.ids <- c('Type I neuroblasts - 9','GMCs - 1','New-born neurons - 8',"Immature neurons","Immature neurons","Immature neurons","Immature neurons","Immature neurons","Immature neurons")
    names(new.cluster.ids) <- levels(NB_GMC_NWB_imn)
    NB_GMC_NWB_imn <- RenameIdents(NB_GMC_NWB_imn, new.cluster.ids)
    
    # GMCs vs T1NB
    GMCs.vs.T1NB.up.list <- c("Bacc","Rbfox1","lncRNA:CR34335","fax", "hdc","RpS27", "tap","RpL38","RpS29","spdo")
    DotPlot(NB_GMC_NWB_imn, features = GMCs.vs.T1NB.up.list, cols = c("grey","red"), col.min = 0, cluster.idents = FALSE, dot.min = 0, assay = "RNA", scale.by = "size") + RotatedAxis() 
    
    GMCs.vs.T1NB.down.list <- c("lncRNA:CR31386","Syp","Hsp27","fru","lncRNA:CR46003","CycB3","wdp","Rcc1", "IntS11","Prosalpha6")
    DotPlot(NB_GMC_NWB_imn, features = GMCs.vs.T1NB.down.list, cols = c("grey","red"), col.min = 0, cluster.idents = FALSE, dot.min = 0, assay = "RNA", scale.by = "size") + RotatedAxis() 
    
    # GMCs vs New-born
    GMCs.vs.NWB.up.list <- c("Pen","edl","Galphai","grh","SoxN","betaTub56D","feo","CNBP","lncRNA:CR45388","E(spl)mgamma-HLH")
    DotPlot(NB_GMC_NWB_imn, features = GMCs.vs.NWB.up.list, cols = c("grey","red"), col.min = 0, cluster.idents = FALSE, dot.min = 0, assay = "RNA", scale.by = "size") + RotatedAxis() 
    
    GMCs.vs.NWB.down.list <- c("hdc","fax","Liprin-gamma","E(spl)m6-BFM","robo2","fne","CG42540","26-29-p","Lim1","CG14989")
    DotPlot(NB_GMC_NWB_imn, features = GMCs.vs.NWB.down.list, cols = c("grey","red"), col.min = 0, cluster.idents = FALSE, dot.min = 0, assay = "RNA", scale.by = "size") + RotatedAxis() 
    
    # NWB vs imn
    NWB.vs.imn.up.list <- c("lncRNA:CR34335","lncRNA:cherub","gsb-n","Bacc","RpL38","lncRNA:CR45388","MFS3","CG10226","robo2","CG1552")
    DotPlot(NB_GMC_NWB_imn, features = NWB.vs.imn.up.list, cols = c("grey","red"), col.min = 0, cluster.idents = FALSE, dot.min = 0, assay = "RNA", scale.by = "size") + RotatedAxis() 
    
    NWB.vs.imn.down.list <- c("Alk", "fz2","trv","lncRNA:noe","Hsp27","Rbp6","mspo","Imp","CG17124","Fas2")
    DotPlot(NB_GMC_NWB_imn, features = NWB.vs.imn.down.list, cols = c("grey","red"), col.min = 0, cluster.idents = FALSE, dot.min = 0, assay = "RNA", scale.by = "size") + RotatedAxis() 
    
    # differentially expressed markers
    defining_makrers.list <- c("CycB3","Rcc1", "IntS11","wdp","Prosalpha6","E(spl)mgamma-HLH","fru","lncRNA:CR46003","Syp","lncRNA:CR31386","Hsp27","Pen","grh","betaTub56D",
                               "CNBP","feo","SoxN","lncRNA:cherub","lncRNA:CR45388","edl",
                               "Galphai","tap","spdo","gsb-n","Bacc","lncRNA:CR34335","RpL38","RpS27",
                               "CG10226","CG1552","MFS3","Lim1","robo2","E(spl)m6-BFM","RpS29",
                               "fax","Liprin-gamma","hdc","fne","CG14989","CG42540",
                               "Rbfox1","Imp","26-29-p","CG17124","mspo","Alk","trv","Fas2","Rbp6","fz2","lncRNA:noe")
    # used in figure 
    NB_GMC_NWB_imn_dotplot <- DotPlot(NB_GMC_NWB_imn, features = defining_makrers.list, cols = c("grey","red"), 
                                      col.min = 0, cluster.idents = FALSE, dot.min = 0, assay = "RNA", scale.by = "size") + RotatedAxis() 
    NB_GMC_NWB_imn_dotplot
    ggsave(filename = "Figure_6/NB_GMC_NWB_imn_dotplot.tiff", plot = NB_GMC_NWB_imn_dotplot, device = "tiff",
           scale = 1, width = 45, height = 11, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
}

# Figure 7 analysis - Mature neurons 
{
  # Sub set mature neurons and recluster
  {
    DotPlot(atlas.combined, features = c("elav","brp","nSyb"), cols = c("grey", "red"), col.min = 0, 
            cluster.idents = FALSE, dot.min = 0, scale.by = "size", dot.scale = 6) + RotatedAxis() 
    # brp and nSyb positive cluster not already identified
    mature.neurons.list <- c(0,2,3,4,5,8,10,12,13,18,21,22,23,28,27,30,34,35,39,46,48,50,51,52,53)
    mature_neurons_unlabeled <- subset(atlas.combined, idents = mature.neurons.list)
    
    # cluster mature neurons (51,596 cells)
    DefaultAssay(mature_neurons_unlabeled) <- "integrated"
    mature_neurons_unlabeled <- RunUMAP(mature_neurons_unlabeled, reduction = "pca", dims = 1:50)
    mature_neurons_unlabeled <- FindNeighbors(mature_neurons_unlabeled, reduction = "pca", dims = 1:50)
    mature_neurons_unlabeled <- FindClusters(mature_neurons_unlabeled, resolution = 0.37) 
    DefaultAssay(mature_neurons_unlabeled) <- "RNA"
    
    # save progress
    saveRDS(mature_neurons_unlabeled, file = "mature_neurons_unlabeled.rds")
    # read RDS file
    mature_neurons_unlabeled <- readRDS("mature_neurons_unlabeled.rds")
  }
  
  # UMAPs (plot used in figure)
  {
      MN_UMAP <- DimPlot(mature_neuron_unlabeled, reduction = "umap", label = FALSE, label.box = FALSE ) + NoLegend() + NoAxes()
      MN_UMAP
      ggsave(filename = "Figure_7/mature_neuron_unlabeled_UMAP.tiff", plot = MN_UMAP, device = "tiff",
             scale = 1, width = 25, height = 25, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
      
      # temporal changes 
      MN_1h <-subset(x = mature_neuron_unlabeled, subset = orig.ident == "atlas_1h_alh" ) 
      MN_24h <-subset(x = mature_neuron_unlabeled, subset = orig.ident == "atlas_24h_alh") 
      MN_48h <-subset(x = mature_neuron_unlabeled, subset = orig.ident == "atlas_48h_alh") 
      
      MN_1h_umap <- DimPlot(MN_1h, reduction = "umap", cols = c("grey","grey","grey","grey","#0CB702","#00BE67","grey",
                                                                "grey","grey","grey","grey","grey","grey","grey"), 
                            label = FALSE, label.box = FALSE) + NoLegend() + NoAxes()
      MN_24h_umap <- DimPlot(MN_24h, reduction = "umap", label = FALSE, cols = c("grey","grey","grey","grey","#0CB702","#00BE67","grey",
                                                                                "grey","grey","grey","grey","grey","grey","grey"),
                             label.box = FALSE) + NoLegend() + NoAxes()
      MN_48h_umap <- DimPlot(MN_48h, reduction = "umap", cols = c("grey","grey","grey","grey","#0CB702","#00BE67","grey",
                                                                  "grey","grey","grey","grey","grey","grey","grey"),
                             label = FALSE, label.box = FALSE) + NoLegend() + NoAxes()
      
      ggsave(filename = "Figure_7/MN_1h_umap.tiff", plot = MN_1h_umap, device = "tiff",
             scale = 1, width = 25, height = 25, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
      ggsave(filename = "Figure_7/MN_24h_umap.tiff", plot = MN_24h_umap, device = "tiff",
             scale = 1, width = 25, height = 25, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
      ggsave(filename = "Figure_7/MN_48h_umap.tiff", plot = MN_48h_umap, device = "tiff",
             scale = 1, width = 25, height = 25, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    }
    
  # Find all markers 
  {
    # Find cluster markers to identify cell identities 
    mature_neuron_unlabeled.markers <- FindAllMarkers(object = mature_neuron_unlabeled, assay = "RNA", 
                                                      logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                                      min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
    write.csv(mature_neuron_unlabeled.markers, file = "mature_neuron_markers.csv")
  }
  
  # Label cluster 
  {
      Mature.neuron.markers.list <- c("CadN","hdc","ChAT","Ace","Gad1","VGlut","Vmat","Ddc","Trh","CCAP","Burs","AstC",
                                      "Tbh","Tdc2","twit","Pka-R1","Pka-R2","Pka-C1","Rgk1","ITP","sNPF")
      DotPlot(mature_neurons_unlabeled, features = Mature.neuron.markers.list, 
              cols = c("grey", "red"), assay = "RNA", col.min = 0, cluster.idents = FALSE, dot.min = 0, scale.by = "size", dot.scale = 6) + RotatedAxis() 
      
      # label and order 
      mature_neurons_labeled <- mature_neurons_unlabeled
      DefaultAssay(mature_neurons_labeled) <- "RNA"
      new.cluster.ids <- c("Cholinergic - 0", "Unannotated - 1", "GABAergic - 2", "Glutamatergic - 3", "Undifferentiated - 4","Undifferentiated - 5", 
                           "Motor neurons - 6","Kenyon cells γ - 7", "Monoaminergic - 8", 
                           "Peptidergic - 9", "Unannotated - 10", "Unannotated - 11", "Octopaminergic - 12", "Neurosecretory cells - 13")
      names(new.cluster.ids) <- levels(mature_neurons_labeled)
      mature_neurons_labeled <- RenameIdents(mature_neurons_labeled, new.cluster.ids)
      order.list <- c( "Undifferentiated - 4","Undifferentiated - 5", "Cholinergic - 0","GABAergic - 2", "Glutamatergic - 3","Monoaminergic - 8","Peptidergic - 9",
                       "Octopaminergic - 12", "Motor neurons - 6", "Kenyon cells γ - 7","Neurosecretory cells - 13", "Unannotated - 1", "Unannotated - 10",
                       "Unannotated - 11")
      mature_neurons_labeled@active.ident <- factor(x = mature_neurons_labeled@active.ident, levels = order.list) 
      
      # save progress
      saveRDS(mature_neurons_labeled, file = "mature_neurons_labeled.rds")
      # read RDS file
      mature_neurons_labeled <- readRDS("mature_neurons_labeled.rds")
    }
    
  # Dot plot
  {
    # cell type defining 
    Mature.neuron.markers.ordered.list <- c("CadN","hdc","ChAT","Ace","Gad1","VGlut","Vmat","Ddc","Trh","CCAP","Burs","AstC",
                                            "Tbh","Tdc2","twit","Pka-R1","Pka-R2","Pka-C1","Rgk1","ITP","sNPF")
    # used in figure 
    cell_ident_dotplot <- DotPlot(mature_neurons_labeled, features = Mature.neuron.markers.ordered.list, 
                                  cols = c("grey", "red"), assay = "RNA", col.min = 0, cluster.idents = FALSE, 
                                  dot.min = 0, scale.by = "size", dot.scale = 6) + RotatedAxis() 
    cell_ident_dotplot
    ggsave(filename = "Figure_7/cell_ident_dotplot.tiff", plot = cell_ident_dotplot, device = "tiff",
             scale = 1, width = 32, height = 12, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
  }
  
  # Deferentially expressed genes between cluster 4 and 5 (undifferentiated)
  {
    # find deferentially expressed genes 
    {
      undiff_neurons_subset <- subset(mature_neurons_labeled, idents = c("Undifferentiated - 4","Undifferentiated - 5"))
      undiff_neurons_temporal.markers <- FindAllMarkers(object = undiff_neurons_subset, assay = "RNA", 
                                                        logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                                        min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
      write.csv(undiff_neurons_temporal.markers, file = "undiff_neurons_temporal_markers.csv")
    }
    
    # Dot plot 
    {
      # top 10 positive genes for each cluster
      cluster4.vs.cluster5.markers.list <- c("Hsp27","cib","Ldh","Liprin-gamma","Toll-7","Toll-6","beat-IIa","dpr9","dsh","Pal2",
                                             "lncRNA:roX2","28SrRNA-Psi:CR40596","Gs2","mt:srRNA","pros","lncRNA:roX1","jim","lncRNA:noe","Obp44a","mt:lrRNA")
      # used in figure 
      undiff_neurons <- DotPlot(undiff_neurons_subset, features = cluster4.vs.cluster5.markers.list, cols = c("grey","red"), col.min = 0, 
                                cluster.idents = FALSE, dot.min = 0, assay = "RNA", scale.by = "size") + RotatedAxis() 
      undiff_neurons
      ggsave(filename = "Figure_7/undiff_neurons_dotplot.tiff", plot = undiff_neurons, device = "tiff",
             scale = 1, width = 22, height = 12, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
    }
  }
  
  # Temporal expression patterns 
  {
      # cluster 0 
      {
        cluster_0 <- subset(mature_neurons_unlabeled, idents = 0)
        Idents(cluster_0) <- "orig.ident"
        cluster_0_subset <- subset(cluster_0, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh"  | orig.ident == "atlas_48h_alh" ))
        cluster0_temporal.markers <- FindAllMarkers(object = cluster_0_subset, assay = "RNA", 
                                                    logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                                    min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
        write.csv(cluster0_temporal.markers, file = "cluster0_temporal_markers.csv")
        
        # Top five positve, temporally deferentially expressed genes from each time point
        cluster0_sub_temp.list <- c("Pde6","CG14312","Rdh","Syn1","Tsp66E",
                                    "lncRNA:cherub","RpL10","RpS4","Hsp27","RpL8",
                                    "ATPsynO","kra", "COX6B","roh","RpS6")
        cluster_0_temp <- DotPlot(cluster_0_subset, features = cluster0_sub_temp.list, cols = c("grey","red"), col.min = 0, 
                                  cluster.idents = FALSE, dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis() 
        cluster_0_temp
      }
      
      # cluster 1
      {
        cluster_1 <- subset(mature_neurons_unlabeled, idents = 1)
        Idents(cluster_1) <- "orig.ident"
        cluster_1_subset <- subset(cluster_1, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))
        cluster0_temporal.markers <- FindAllMarkers(object = cluster_1_subset, assay = "RNA", 
                                                    logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                                    min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
        write.csv(cluster0_temporal.markers, file = "cluster1_temporal_markers.csv")
        
        # Top five positve, temporally deferentially expressed genes from each time point
        cluster1_sub_temp.list <- c("Hsp67Bc","Bacc","sif","CG17124","Tpi",
                                    "RpS17","RpL36A","RpL32","RpS29", "RpL36",
                                    "lncRNA:cherub","lncRNA:noe","Eip74EF","pros","hdc")
        cluster_1_temp <- DotPlot(cluster_1_subset, features = cluster1_sub_temp.list, cols = c("grey","red"), col.min = 0, 
                                  cluster.idents = FALSE, dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis() 
        cluster_1_temp
      }
      
      # cluster 2 
      {
        cluster_2 <- subset(mature_neurons_unlabeled, idents = 2)
        Idents(cluster_2) <- "orig.ident"
        cluster_2_subset <- subset(cluster_2, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))
        cluster2_temporal.markers <- FindAllMarkers(object = cluster_2_subset, assay = "RNA", 
                                                    logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                                    min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
        write.csv(cluster2_temporal.markers, file = "cluster2_temporal_markers.csv")
        
        # Top five positve, temporally deferentially expressed genes from each time point
        cluster2_sub_temp.list <- c("Cngl","Drep2","CG17124","CG30172","sif",
                                    "mt:ATPase6","mt:CoII","lncRNA:cherub","mt:CoI","RpS4",
                                    "lncRNA:noe","lncRNA:roX1","Top1")
        cluster_2_temp <- DotPlot(cluster_2_subset, features = cluster2_sub_temp.list, cols = c("grey","red"), col.min = 0, 
                                  cluster.idents = FALSE, dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis() 
        cluster_2_temp
      }
      
      # cluster 3
      {
        cluster_3 <- subset(mature_neurons_unlabeled, idents = 3)
        Idents(cluster_3) <- "orig.ident"
        cluster_3_subset <- subset(cluster_3, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))
        cluster3_temporal.markers <- FindAllMarkers(object = cluster_3_subset, assay = "RNA", 
                                                    logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                                    min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
        write.csv(cluster3_temporal.markers, file = "cluster3_temporal_markers.csv")
        
        # Top five positve, temporally deferentially expressed genes from each time point
        cluster3_sub_temp.list <- c("Hex-A","Cngl","Bacc","CG30172","Drep2",
                                    "Ldh","lncRNA:cherub","His3.3B","betaTub56D","mt:ATPase6",
                                    "lncRNA:noe","Eip74EF","VGlut","Msp300","hdc")
        cluster_3_temp <- DotPlot(cluster_3_subset, features = cluster3_sub_temp.list, cols = c("grey","red"), col.min = 0, 
                                  cluster.idents = FALSE, dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis() 
        cluster_3_temp
      }
      
      # cluster 6 
      {
        cluster_6 <- subset(mature_neurons_unlabeled, idents = 6)
        Idents(cluster_6) <- "orig.ident"
        cluster_6_subset <- subset(cluster_6, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))
        cluster6_temporal.markers <- FindAllMarkers(object = cluster_6_subset, assay = "RNA", 
                                                    logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                                    min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
        write.csv(cluster6_temporal.markers, file = "cluster6_temporal_markers.csv")
        
        # Top five positve, temporally deferentially expressed genes from each time point
        cluster6_sub_temp.list <- c("fax","Drep2","Bacc","Cngl","bel",
                                    "RpL36","mt:ND4","RpL37A","mt:ND5","RpS17",
                                    "Msp300","lncRNA:noe","Proc","Eip74EF","CG6329")
        cluster_6_temp <- DotPlot(cluster_6_subset, features = cluster6_sub_temp.list, cols = c("grey","red"), col.min = 0, 
                                  cluster.idents = FALSE, dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis() 
        cluster_6_temp
      }
      
      # cluster 7
      {
        cluster_7 <- subset(mature_neurons_unlabeled, idents = 7)
        Idents(cluster_7) <- "orig.ident"
        cluster_7_subset <- subset(cluster_7, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))
        cluster7_temporal.markers <- FindAllMarkers(object = cluster_7_subset, assay = "RNA", 
                                                    logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                                    min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
        write.csv(cluster7_temporal.markers, file = "cluster7_temporal_markers.csv")
        
        # Top five positve, temporally deferentially expressed genes from each time point
        cluster7_sub_temp.list <- c("rad","UQCR-11","VhaM9.7-a","CG14312","Nlg3",
                                    "His3.3B","cib", "Hsp26","Hsp27","lncRNA:CR40469")
        cluster_7_temp <- DotPlot(cluster_7_subset, features = cluster7_sub_temp.list, cols = c("grey","red"), col.min = 0, 
                                  cluster.idents = FALSE, dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis() 
        cluster_7_temp
      }
      
      # cluster 8
      {
        cluster_8 <- subset(mature_neurons_unlabeled, idents = 8)
        Idents(cluster_8) <- "orig.ident"
        cluster_8_subset <- subset(cluster_8, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))
        cluster8_temporal.markers <- FindAllMarkers(object = cluster_8_subset, assay = "RNA", 
                                                    logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                                    min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
        write.csv(cluster8_temporal.markers, file = "cluster8_temporal_markers.csv")
        
        # Top five positve, temporally deferentially expressed genes from each time point
        cluster8_sub_temp.list <- c("CG30172","Hex-A","Drep2","Tsp66E",
                                    "Ldh","cib","His3.3B","lncRNA:cherub","betaTub56D",
                                    "Eip74EF","tau","lncRNA:noe")
        cluster_8_temp <- DotPlot(cluster_8_subset, features = cluster8_sub_temp.list, cols = c("grey","red"), col.min = 0, 
                                  cluster.idents = FALSE, dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis() 
        cluster_8_temp
      }
      
      # cluster 9
      {
        cluster_9 <- subset(mature_neurons_unlabeled, idents = 9)
        Idents(cluster_9) <- "orig.ident"
        cluster_9_subset <- subset(cluster_9, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))
        cluster9_temporal.markers <- FindAllMarkers(object = cluster_9_subset, assay = "RNA", 
                                                    logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                                    min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
        write.csv(cluster9_temporal.markers, file = "cluster9_temporal_markers.csv")
        
        # Top five positve, temporally deferentially expressed genes from each time point 
        cluster9_sub_temp.list <- c("Bacc","pum","Hex-A","cdc14","bel",
                                    "Ldh", "lncRNA:CR40469","His3.3B","Hsp27","RpL3",
                                    "Hsp23","lncRNA:noe","lncRNA:cherub")
        cluster_9_temp <- DotPlot(cluster_9_subset, features = cluster9_sub_temp.list, cols = c("grey","red"), col.min = 0, 
                                  cluster.idents = FALSE, dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis() 
        cluster_9_temp
      }
      
      # cluster 10
      {
        cluster_10 <- subset(mature_neurons_unlabeled, idents = 10)
        Idents(cluster_10) <- "orig.ident"
        cluster_10_subset <- subset(cluster_10, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))
        cluster10_temporal.markers <- FindAllMarkers(object = cluster_10_subset, assay = "RNA", 
                                                     logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                                     min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
        write.csv(cluster10_temporal.markers, file = "cluster10_temporal_markers.csv")
        
        # Top five positve, temporally deferentially expressed genes from each time point 
        cluster10_sub_temp.list <- c("Drep2","CG30172","Cngl","Pde6","Syx1A",
                                     "cib", "Df31", "Hsp27","lncRNA:CR40469", "fabp",
                                     "lncRNA:noe","lncRNA:cherub","hdc")
        cluster_10_temp <- DotPlot(cluster_10_subset, features = cluster10_sub_temp.list, cols = c("grey","red"), col.min = 0, 
                                   cluster.idents = FALSE, dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis() 
        cluster_10_temp
      }
      
      # cluster 11
      {
        cluster_11 <- subset(mature_neurons_unlabeled, idents = 11)
        Idents(cluster_11) <- "orig.ident"
        cluster_11_subset <- subset(cluster_11, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))
        cluster11_temporal.markers <- FindAllMarkers(object = cluster_11_subset, assay = "RNA", 
                                                     logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                                     min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
        write.csv(cluster11_temporal.markers, file = "cluster11_temporal_markers.csv")
        
        # Top five positve, temporally deferentially expressed genes from each time point 
        cluster11_sub_temp.list <- c("Cngl","bel","Hex-A","nAChRalpha6","msi",
                                     "Ldh","Hsc70Cb", "RpL3", "Hsp27","Pal2",
                                     "lncRNA:noe")
        cluster_11_temp <- DotPlot(cluster_11_subset, features = cluster11_sub_temp.list, cols = c("grey","red"), col.min = 0, 
                                   cluster.idents = FALSE, dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis() 
        cluster_11_temp
      }
      
      # cluster 12
      {
        cluster_12 <- subset(mature_neurons_unlabeled, idents = 12)
        Idents(cluster_12) <- "orig.ident"
        cluster_12_subset <- subset(cluster_12, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))
        cluster12_temporal.markers <- FindAllMarkers(object = cluster_12_subset, assay = "RNA", 
                                                     logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                                     min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
        write.csv(cluster12_temporal.markers, file = "cluster12_temporal_markers.csv")
        
        # Top five positve, temporally deferentially expressed genes from each time point
        cluster12_sub_temp.list <- c("bel","Eip75B","luna","Bacc","Dsp1",
                                     "RpS23","fabp","RpS30","RpS18","mt:CoI",
                                     "Msp300", "hdc", "Arc1", "Hsp23","lncRNA:noe")
        cluster_12_temp <- DotPlot(cluster_12_subset, features = cluster12_sub_temp.list, cols = c("grey","red"), col.min = 0, 
                                   cluster.idents = FALSE, dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis() 
        cluster_12_temp
      }
      
      # cluster 13
      {
        cluster_13 <- subset(mature_neurons_unlabeled, idents = 13)
        Idents(cluster_13) <- "orig.ident"
        cluster_13_subset <- subset(cluster_13, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))
        cluster13_temporal.markers <- FindAllMarkers(object = cluster_13_subset, assay = "RNA", 
                                                     logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                                     min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
        write.csv(cluster13_temporal.markers, file = "cluster13_temporal_markers.csv")
        
        # Top five positve, temporally deferentially expressed genes from each time point
        cluster13_sub_temp.list <- c("CG34200","CG12384","lncRNA:CR44320","CG13585","Atg1")
        cluster_13_temp <- DotPlot(cluster_13_subset, features = cluster13_sub_temp.list, cols = c("grey","red"), col.min = 0, 
                                   cluster.idents = FALSE, dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis() 
        cluster_13_temp
      }
      
      # Differentiated neurons - temporal changes 
      {
        # sub set out labeled clusters 
        diff_neurons_subset <- subset(mature_neurons_unlabeled, idents = c(0,2,3,8,9,12,6,7,13))
        diff_neurons_subset_temp <- diff_neurons_subset
        Idents(diff_neurons_subset_temp) <- "orig.ident"
        diff_neurons_subset_temp <- subset(diff_neurons_subset_temp, subset = (orig.ident == "atlas_1h_alh" | orig.ident == "atlas_24h_alh" | orig.ident == "atlas_48h_alh" ))
        diff_neurons_temporal.markers <- FindAllMarkers(object = diff_neurons_subset_temp, assay = "RNA", 
                                                        logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.1,
                                                        min.diff.pct = 0.1,verbose = TRUE, return.thresh = 0.0501)
        write.csv(diff_neurons_temporal.markers, file = "diff_neurons_temporal_markers.csv")
        
        # Top temporally deferentially expressed genes from other lists with multiple hits (at least 2/9 as top five)
        diff_neurons_sub_temp.list <- c("CG14312","Tsp66E","Bacc", "CG30172","bel","Cngl", "Drep2", "Hex-A",
                                          "His3.3B","Ldh","cib","lncRNA:CR40469","Hsp27","betaTub56D","mt:ATPase6","mt:CoI","RpS4","lncRNA:cherub",
                                          "hdc", "Hsp23", "Msp300","lncRNA:noe","Eip74EF")
        # used in figure
        diff_neurons_temp <- DotPlot(diff_neurons_subset_temp, features = diff_neurons_sub_temp.list, cols = c("grey","red"), col.min = 0, 
                                       cluster.idents = FALSE, dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis() 
        diff_neurons_temp
        ggsave(filename = "Figure_7/diff_neurons_temp_dotplot.tiff", plot = diff_neurons_temp, device = "tiff",
               scale = 1, width = 32, height = 12, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
        
        # Top temporally deferentially expressed cell surface molecule genes from other lists with multiple hits (3/9 clusters) 
        diff_neurons_sub_temp_CSM.list <- c("5-HT1A","AstC-R2","DopEcR","GABA-B-R1","GluRIA","GluRIB","kek1", "nkd","Nlg2","Nlg3","side",
                                            "Dg","fz2","lbk","smog","Toll-7","trn","Alk","Eph")
        # used in figure
        diff_neurons_temp_CSM <- DotPlot(diff_neurons_subset_temp, features = diff_neurons_sub_temp_CSM.list, cols = c("grey","red"), col.min = 0, 
                                       cluster.idents = FALSE, dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis() 
        diff_neurons_temp_CSM
        ggsave(filename = "Figure_7/diff_neurons_temp_dotplot_CSM.tiff", plot = diff_neurons_temp_CSM, device = "tiff",
               scale = 1, width = 20, height = 12, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
        
        # Top temporally deferentially expressed cell surface molecule genes from other lists with multiple hits (3/9 clusters)  
        diff_neurons_sub_temp_TF.list <- c("salr","cg","maf-S","CG10543","Glut4EF","Eip78C", "psq","cic", "gro","crol","Eip75B","luna", "Pits","scrt",
                                          "hth","bi","chinmo", "Pdp1","SoxN","tara","Dp","BtbVII","crc","CtBP","kra", "lola","Nacalpha","Pfk", "sr","bic",
                                          "pros","Eip74EF")
        # used in figure
        diff_neurons_temp_TF <- DotPlot(diff_neurons_subset_temp, features = diff_neurons_sub_temp_TF.list, cols = c("grey","red"), col.min = 0, 
                                       cluster.idents = FALSE, dot.min = 0, assay = "RNA", group.by = "orig.ident", scale.by = "size") + RotatedAxis() 
        diff_neurons_temp_TF
        ggsave(filename = "Figure_7/diff_neurons_temp_dotplot_TTF.tiff", plot = diff_neurons_temp_TF, device = "tiff",
               scale = 1, width = 35, height = 12, units = c("cm"), dpi = 300, limitsize = TRUE, bg = "white")
      }
  }
  
}