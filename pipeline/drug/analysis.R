library(argparse)
library(future)
options(future.globals.maxSize = Inf)

parser <- ArgumentParser()

parser$add_argument("--count_matrix", 
                    required=TRUE,
                    help="Count matrix.")
parser$add_argument("--group_info", 
                    required=TRUE,
                    help="Group info file.")
parser$add_argument("--project", 
                    required=TRUE,
                    help="Project name.")
parser$add_argument("--group_by", 
                    default='plate',
                    help="Group by.")
parser$add_argument("--outdir", 
                    required=TRUE,
                    help="Output directory.")
parser$add_argument("--method", 
                    help="Normalize method for integration.",
                    choices=c('integrate', 'SCTransform'),
                    default='integrate')
parser$add_argument("--min_cells", 
                    default=4,
                    help="Include features detected in at least this many cells. 
                    Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.")
parser$add_argument("--min_features", 
                    default=200,
                    help="Include cells where at least this many features are detected.")
parser$add_argument("--dims", 
                    default=6,
                    help="PCA nums.")
parser$add_argument("--k_filter", 
                    default=80,
                    help="Mininum cell nums of sample.")
args <- parser$parse_args()


library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)

data_process <- function(matrix_file, group_info) {
  counts <- read.table(matrix_file,
                       sep='\t', header=T,
                       check.names = F)
  
  rownames(counts) <- counts$gene_name
  counts <- counts[-c(1,2)]
  
  group_info <- read.table(group_info,
                         sep='\t',
                         header = T)
  group_info$barcode <- paste(group_info$barcode, 
                              group_info$plate_tag, 
                              sep = '-')
  group_info <- group_info[group_info$barcode %in% colnames(counts), ]
  
  return(list(counts, group_info))
}


seurat_integrate <- function(counts, 
                             group_data,
                             project, 
                             outdir,
                             group_by='plate', 
                             min.cells=4,
                             min.features=200,
                             k.filter=80, 
                             dims=6) {
  # make seurat object
  pro <- CreateSeuratObject(counts = counts, project = project, 
                            min.cells = min.cells,
                            min.features = min.features)
  # count mito percent
  pro[["percent.mt"]] <- PercentageFeatureSet(pro, pattern = "^MT-")
  # filter
  
  pro <- subset(pro, subset = nCount_RNA > 65536 & nFeature_RNA > 6000 & percent.mt < 10)

  # match group
  group_data <- group_data[group_data$barcode %in% rownames(pro@meta.data), ]
  group_data <- group_data[match(rownames(pro@meta.data), group_data$barcode), ]
  
  pro@meta.data$plate <- group_data$plate_tag
  pro@meta.data$treatment <- group_data$treat
  n_features = nrow(pro[['RNA']]@counts)
  
  # split and integrate
  pro.list <- SplitObject(pro, split.by = group_by)
  # normalize and identify variable features for each dataset independently
  pro.list <- lapply(X = pro.list, FUN = function(x) {
    x <- NormalizeData(x, verbose=FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",
                              nfeatures = n_features, verbose=FALSE)
  })
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = pro.list, verbose=FALSE)
  pro.anchors <- FindIntegrationAnchors(object.list = pro.list, 
                                        anchor.features = features,
                                        k.filter = k.filter,
                                        verbose=FALSE)
  # this command creates an 'integrated' data assay
  pro.combined <- IntegrateData(anchorset = pro.anchors,
                                k.weight = k.filter,
                                verbose=FALSE)
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(pro.combined) <- "integrated"
  
  # normal run
  # Run the standard workflow for visualization and clustering
  pro.combined <- ScaleData(pro.combined, 
                            verbose=FALSE)
  
  pro.combined <- RunPCA(pro.combined, 
                         npcs = dims, 
                         verbose=FALSE)
  pro.combined <- RunUMAP(pro.combined, 
                          reduction = "pca", 
                          dims = 1:dims,
                          verbose=FALSE)
  pro.combined <- RunTSNE(pro.combined, 
                          reduction = "pca", 
                          dims = 1:dims,
                          verbose=FALSE)
  pro.combined <- FindNeighbors(pro.combined, 
                                reduction = "pca", 
                                dims = 1:dims,
                                verbose=FALSE)
  pro.combined <- FindClusters(pro.combined, 
                               resolution = 0.5,
                               verbose=FALSE)
  
  saveRDS(pro.combined,
          file = paste(outdir, '/integrate.seurat.rds', sep=''))
  return(pro.combined)
}


# SCT
seurat_sctransform <- function(counts, 
                               group_data,
                               project,
                               outdir,
                               group_by='plate',
                               min.cells=4,
                               min.features=200,
                               k.filter=80,
                               dims=6) {

  # make seurat object
  pro <- CreateSeuratObject(counts = counts, 
                            project = project, 
                            min.cells = min.cells,
                            min.features = min.features)
  # count mito percent
  pro[["percent.mt"]] <- PercentageFeatureSet(pro, pattern = "^MT-")
  # filter
  pro <- subset(pro, subset = nCount_RNA > 65536 & nFeature_RNA > 6000 & percent.mt < 10)

  # match group
  group_data <- group_data[group_data$barcode %in% rownames(pro@meta.data), ]
  group_data <- group_data[match(rownames(pro@meta.data), group_data$barcode), ]
  
  pro@meta.data$plate <- group_data$plate_tag
  pro@meta.data$treatment <- group_data$treat
  n_features = nrow(pro[['RNA']]@counts)

  # split and integrate
  pro.list <- SplitObject(pro, split.by = group_by)

  pro.list <- lapply(X = pro.list, FUN = SCTransform)

  features <- SelectIntegrationFeatures(object.list = pro.list, 
                                        nfeatures = n_features,
                                        verbose=FALSE)

  pro.list <- PrepSCTIntegration(object.list = pro.list, 
                                 anchor.features = features,
                                 verbose=FALSE)

  pro.anchors <- FindIntegrationAnchors(object.list = pro.list, 
                                        normalization.method = "SCT",
                                        anchor.features = features,
                                        k.filter = k.filter,
                                        verbose=FALSE)

  pro.combined.sct <- IntegrateData(anchorset = pro.anchors, 
                                    normalization.method = "SCT",
                                    k.weight = k.filter,
                                    verbose=FALSE)
  
  
  pro.combined.sct <- RunPCA(pro.combined.sct, 
                             npcs = dims,
                             verbose=FALSE)
  pro.combined.sct <- RunUMAP(pro.combined.sct, 
                              reduction = "pca", 
                              dims = 1:dims,
                              verbose=FALSE)
  pro.combined.sct <- RunTSNE(pro.combined.sct, 
                              reduction = "pca", 
                              dims = 1:dims,
                              verbose=FALSE)
  pro.combined.sct <- FindNeighbors(pro.combined.sct, 
                                reduction = "pca", 
                                dims = 1:dims,
                                verbose=FALSE)
  pro.combined.sct <- FindClusters(pro.combined.sct, 
                               resolution = 0.5,
                               verbose=FALSE)
  
  saveRDS(pro.combined.sct, 
          file = paste(outdir, '/sct.seurat.rds', sep=''))
  return(pro.combined.sct)
}


plt_module <- function(seurat_data, 
                       outdir,
                       height=25,
                       width=40,
                       units='cm') {
  seurat_obj <- seurat_data
  seurat_obj@meta.data <- rename(seurat_obj@meta.data, 
                                 gene_count=nFeature_RNA, 
                                 UMI_count=nCount_RNA)
  
  vln1 <- VlnPlot(seurat_obj,
          features = c("gene_count", "UMI_count", "percent.mt"),
          assay='RNA',
          group.by = 'orig.ident',
          ncol = 3,
          combine = TRUE,
          pt.size = 0)
  ggsave(filename = paste(outdir, '/vln_plot1.png', sep=''),
         plot = vln1,
         width = width,
         height = height,
         units = units)

  vln2 <- VlnPlot(seurat_obj,
                  features = c("gene_count", "UMI_count"),
                  assay='RNA',
                  group.by = 'plate',
                  ncol = 2,
                  pt.size = 0) +
    geom_boxplot(width=.2,col="black",fill="white") +
    geom_hline(yintercept = 65536, color='red')
  
  ggsave(filename = paste(outdir, '/vln_plot2.png', sep=''),
         plot = vln2,
         width = width,
         height = height,
         units = units)



  # pca percent
  percentage<-round(seurat_obj@reductions$pca@stdev / sum(seurat_obj@reductions$pca@stdev) * 100,2)
  percentage<-paste(colnames(seurat_obj@reductions$pca@cell.embeddings),"(", paste(as.character(percentage), "%", ")", sep=""))
  
  pca1 <- DimPlot(seurat_obj, reduction='pca', group.by = 'plate') + 
    xlab(percentage[1]) +
    ylab(percentage[2])
  pca2 <- DimPlot(seurat_obj, reduction='pca', group.by = 'treatment') + 
    xlab(percentage[1]) +
    ylab(percentage[2])
  pca3 <- DimPlot(seurat_obj, reduction='pca', group.by = 'seurat_clusters') + 
    xlab(percentage[1]) +
    ylab(percentage[2])
  pca <- pca1+pca2+pca3
  ggsave(filename = paste(outdir, '/pca_plot.png', sep=''),
         plot = pca,
         width = width,
         height = height,
         units = units)
  
  p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "plate")
  p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = 'seurat_clusters',
                label = TRUE, repel = TRUE)
  p3 <- DimPlot(seurat_obj, reduction = "umap", group.by = "treatment")
  umap_plt <- p1+p2+p3
  ggsave(filename = paste(outdir, '/umap_plt.png', sep=''),
         plot = umap_plt,
         width = width,
         height = height,
         units = units)
  
  p1 <- DimPlot(seurat_obj, reduction = "tsne", group.by = "plate")
  p2 <- DimPlot(seurat_obj, reduction = "tsne", group.by = 'seurat_clusters',
                label = TRUE, repel = TRUE)
  p3 <- DimPlot(seurat_obj, reduction = "tsne", group.by = "treatment")
  tsne_plt <- p1+p2+p3
  ggsave(filename = paste(outdir, '/tsne_plt.png', sep=''),
         plot = tsne_plt,
         width = width,
         height = height,
         units = units)

  
}


run <- function(args) {
  if (!dir.exists(args$outdir)) {
    dir.create(args$outdir)
  }

  data_list <- data_process(matrix_file = args$count_matrix,
                            group_info = args$group_info)
  if (args$method=='SCTransform') {
    seurat_obj <- seurat_sctransform(counts = data_list[[1]],
                       group_data = data_list[[2]],
                       outdir = args$outdir,
                       project = args$project,
                       group_by = args$group_by,
                       dims = as.integer(args$dims),
                       min.cells = as.integer(args$min_cells),
                       min.features = as.integer(args$min_features),
                       k.filter = as.integer(args$k_filter))
  } else if (args$method=='integrate') {
      seurat_obj <- seurat_integrate(counts = data_list[[1]],
                                     group_data = data_list[[2]],
                                     outdir = args$outdir,
                                     project = args$project,
                                     group_by = args$group_by,
                                     dims = as.integer(args$dims),
                                     min.cells = as.integer(args$min_cells),
                                     min.features = as.integer(args$min_features),
                                     k.filter = as.integer(args$k_filter))
  }


  plt_module(seurat_data = seurat_obj,
             outdir=args$outdir)
}

run(args)



