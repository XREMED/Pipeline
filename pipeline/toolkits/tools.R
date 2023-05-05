library(argparse)
library(future)
options(future.globals.maxSize = Inf)

library(ggplot2, quietly = TRUE)
library(patchwork, quietly = TRUE)
library(tidyverse, quietly = TRUE)

seurat_integrated <- function(counts, 
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
  
  pro@meta.data$tag <- group_data$tag
  pro@meta.data$treatment <- group_data$treatment
  # pro@meta.data$library <- sapply(pro@meta.data$tag,
  #                                 FUN = function(x) {
  #                                   s <- strsplit(x, '-', fixed = T)[[1]]
  #                                   return(s[2])
  #                                 })
  pro@meta.data$library <- group_data$library
  pro@meta.data$plate <- sapply(pro@meta.data$tag,
                                FUN = function(x) {
                                  s <- strsplit(x, '-', fixed = T)[[1]]
                                  return(s[2])
                                })
  pro@meta.data$detail_treatment <- group_data$detail_treatment
  # n_features = nrow(pro[['RNA']]@counts)
  n_features=5000
  # filter plate
  plate_df <- data.frame(table(pro@meta.data$tag))
  plate_df <- plate_df[plate_df$Freq>=30,]
  pro <- subset(pro, tag %in% plate_df$Var1)
  # split and integrated
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
  saveRDS(pro.combined, file=paste(outdir, '/seurat.integrated.rds', sep=''))
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
          file = paste(outdir, '/seurat.integrated.rds', sep=''))
  return(pro.combined)
}


score_compound <- function(seurat_obj,
                           outdir='./',
                           method='integrated',
                           NC='DMSO',
                           PC='PC') {
  # set data
  normalized_data <- data.frame(seurat_obj[[method]]@data, check.names = F)
  meta.data <- seurat_obj@meta.data
  # compute nc/pc median
  nc_samples <- rownames(meta.data[meta.data$treatment==NC,])
  nc_median <- apply(normalized_data[, nc_samples], 1, FUN = function(x) {
    return(median(x))
  })
  nc_median[nc_median==0] <- 1
  pc_samples <- rownames(meta.data[meta.data$treatment==PC,])
  pc_median <- apply(normalized_data[, pc_samples], 1, FUN = function(x) {
    return(median(x))
  })
  # compute gene weight
  weight <- data.frame(weight=abs(log2((pc_median/nc_median)+1)-1))
  
  # pc-nc diff
  pc_minus_nc <- pc_median-nc_median
  # compute score
  
  score <- apply(normalized_data, 2, FUN = function(x) {
    e <- (x-nc_median)/(pc_minus_nc)
    e[which(is.infinite(e))] <- 0
    e <- e*weight$weight
    e[is.na(e)] <- 0
    s <- sum(e)
    return(s)
  })
  # sort score and normalize score
  score_df <- data.frame(score=score)
  
  max_s <- max(score_df$score)
  min_s <- min(score_df$score)
  
  score_df$normalized_score <- sapply(score_df$score, FUN = function(x) {
    ifelse(x>=0, x/max_s, -1*x/min_s)
  })
  # add treat info
  score_df <- score_df[order(-score_df$normalized_score),]
  score_df$num <- c(1:nrow(score_df)) 
  score_df <- score_df[match(rownames(meta.data), rownames(score_df)),]
  score_df$treatment <- meta.data$treatment
  score_df$detail_treatment <- meta.data$detail_treatment
  score_df$seurat_clusters <- meta.data$seurat_clusters
  
  # plot
  score_plot <- ggplot(score_df, mapping = aes(num, 
                                               normalized_score, 
                                               color=detail_treatment)) +
    geom_point(size=1) + 
    ggtitle('Compounds score') +
    labs(x="Compounds tag number",y="Score") +
    theme(plot.title = element_text(hjust = 0.5))
  
  # save results
  write.csv(score_df, file = paste(outdir, '/score.csv', sep=''))
  ggsave(filename = paste(outdir, '/score_plot.png', sep=''),
         plot = score_plot,
         width = 40,
         height = 28,
         units = 'cm')
  ggsave(filename = paste(outdir, '/score_plot.pdf', sep=''),
         plot = score_plot,
         width = 40,
         height = 28,
         units = 'cm')
  
  return(score_df)
  
}

merge_df <- function(df1, df2) {
  df <- merge(df1, df2, by='Genes', all=TRUE)
  return(df)
}

check_dir <- function(outdir) {
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
}


plt_module <- function(seurat_data, 
                       outdir,
                       height=25,
                       width=40,
                       units='cm') {
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
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
  vln1[[2]] <- vln1[[2]] + 
    scale_y_continuous(breaks=seq(0, max(seurat_obj@meta.data$UMI_count), 100000))
  
  ggsave(filename = paste(outdir, '/vln_all.png', sep=''),
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
  vln2[[2]] <- vln2[[2]] + 
    scale_y_continuous(breaks=seq(0, max(seurat_obj@meta.data$UMI_count), 100000))
  
  ggsave(filename = paste(outdir, '/vln_plate.png', sep=''),
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
  pca2 <- DimPlot(seurat_obj, reduction='pca', group.by = 'detail_treatment') + 
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
  p3 <- DimPlot(seurat_obj, reduction = "umap", group.by = "detail_treatment")
  umap_plt <- p1+p2+p3
  ggsave(filename = paste(outdir, '/umap_plt.png', sep=''),
         plot = umap_plt,
         width = width,
         height = height,
         units = units)
  
  p1 <- DimPlot(seurat_obj, reduction = "tsne", group.by = "plate")
  p2 <- DimPlot(seurat_obj, reduction = "tsne", group.by = 'seurat_clusters',
                label = TRUE, repel = TRUE)
  p3 <- DimPlot(seurat_obj, reduction = "tsne", group.by = "detail_treatment")
  tsne_plt <- p1+p2+p3
  ggsave(filename = paste(outdir, '/tsne_plt.png', sep=''),
         plot = tsne_plt,
         width = width,
         height = height,
         units = units)
  
  
}

no_rep_deg <- function(counts, group, outdir, test, bcv=0.4) {
  y <- DGEList(counts = counts,genes = rownames(counts),group = group)
  keep <- rowSums(log2(cpm(y))>=2) >= 1
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y_bcv <- y
  et <- exactTest(y_bcv, dispersion = bcv ^ 2)
  degs <- decideTestsDGE(et, p.value = 0.05, lfc = 1)
  
  results <- cbind(y$genes,y$counts,et$table,degs)
  
  results$expression <- ifelse(abs(results$logFC)<2 | results$PValue>0.05, 'NoDiff', 
                               ifelse(results$logFC>0&results$PValue<=0.05, 'Up', 'Down'))
  
  volcano <- ggplot(results, mapping = aes(x=logFC,y=-log10(PValue), color=expression)) +
    theme_bw()+
    labs(x="LogFC",y="-Log10PValue")+
    geom_hline(yintercept = -log10(0.05), linetype='dashed')+
    geom_vline(xintercept = c(-2, 2), linetype='dashed')+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=c('Up'="red",'Down'="green",'NoDiff'="grey"),
                       name="Regulation",breaks=c("Up","Down","NoDiff"))+##！！！图例的说明需要适当修改
    ggtitle(paste(test, 'differential expressed genes', sep=' '))+
    geom_point()
  ggsave(filename = paste(outdir, '/',test, '_DEGs_volcano.png', sep=''),
         plot = volcano,
         width = 40,
         height = 28,
         units = 'cm')
  
  write.csv(results, file = paste(outdir, '/', test, '_DEGs.csv', sep=''),
            quote = F)
}

## plot enrichment analysis 
PlotBubble <- function(file,              # string, GO/KEGG results file.
                       x_lab='GeneRatio', # string, x-axis label.
                       y_lab='Category',  # string, y-axis label.
                       plot_title='Gene Enrichment', # string, plot title.
                       topn=10, # integer, number of category to show.
                       ...) {
  df <- read.csv(file, header = T)
  new_df <- subset(df, df$PValue<=0.05)
  print(dim(new_df))
  new_df <- new_df[order(new_df$PValue), ]
  if (nrow(new_df)>=topn) {
    plot_data <- new_df[c(1:topn), ]
  } else {
    plot_data <- new_df
  }
  bubble_plot <- ggplot(data=plot_data, aes(x=X., y=Term)) + 
    geom_point(aes(color=-log10(PValue), size=Count), alpha=0.8) + 
    theme_bw(base_size=15) + 
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank(),
      plot.title=element_text(hjust=0.5)
    ) + 
    ggtitle(plot_title) + 
    xlab(x_lab) +
    ylab(y_lab) +
    scale_y_discrete(labels=function(x) str_wrap(x, width=50)) +
    scale_colour_gradient(low="purple",high="red") +
    scale_size_area(max_size = 20)
  return(bubble_plot)
}
