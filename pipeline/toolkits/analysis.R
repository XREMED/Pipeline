library(Seurat, quietly = TRUE, verbose = FALSE)
library(DESeq2, quietly = TRUE, verbose = FALSE)
library(ComplexHeatmap, quietly = TRUE, verbose = FALSE)
library(SeuratDisk, quietly = TRUE, verbose = FALSE)
library(future, quietly = TRUE, verbose = FALSE)
library(argparse, quietly = TRUE, verbose = FALSE)
library(ggplot2, quietly = TRUE, verbose = FALSE)

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
                    default='tag',
                    help="Group by.")
parser$add_argument("--outdir",
                    required=TRUE,
                    help="Output directory.")
parser$add_argument("--normalization.method",
                    help="Normalize method for integration.",
                    choices=c('LogNormalize', 'SCT'),
                    default='LogNormalize')
parser$add_argument("--dims",
                    default=NULL,
                    help="PCA nums.")
parser$add_argument("--k_filter",
                    default=30,
                    help="Mininum cell nums of sample.")
args <- parser$parse_args()


check_dir <- function(outdir) {
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
}

check_file <- function(filename) {
  if (!file.exists(filename)) {
    stop(paste('ERROR: ', filename, ' not found!', sep=''))
  }
}

user.input <- function(prompt) {
  if (interactive()) {
    return(readline(prompt))
  } else {
    cat(prompt)
    return(readLines("stdin", n=1))
  }
}

read_countData <- function(count_file) {
  message('Reading in count file and process count data...')
  check_file(count_file)
  countData <- read.table(count_file,
                          sep='\t',
                          header = T,
                          check.names = F)
  rownames(countData) <- countData$gene_name
  countData <- countData[, -c(1:2)]
  colnames(countData) <- toupper(colnames(countData))
  return(countData)
}

read_colData <- function(group_info) {
  message('Reading in group info and process group data...')
  check_file(group_info)
  colData <- read.table(group_info, header = T,
                        sep='\t')
  colData$count_label <- paste(colData$barcode, '-', colData$tag, sep='')
  colData$plate_label <- paste(colData$well, '-', colData$tag, sep='')
  colData <- colData[colData$treatment!='BLANK',]
  colData <- data.frame(apply(colData, 2, function(x){toupper(x)}))
  return(colData)
}


make_seurat_list <- function(count_file,
                             group_info,
                             filename=NULL,
                             outdir='./',
                             min.good.wells=30,
                             split.by='tag') {
  countData <- read_countData(count_file)
  colData <- read_colData(group_info)
  message('Making seurat object list...')
  # check dir
  check_dir(outdir)
  # subset colData and countData
  both_samples <- intersect(colData$count_label, colnames(countData))
  colData <- colData[colData$count_label%in%both_samples,]
  colData <- colData[match(both_samples, colData$count_label),]
  countData <- countData[, both_samples]
  ### make seurat object
  obj <- CreateSeuratObject(counts = countData,
                            project=project,
                            min.cells = min.cells,
                            min.features = min.features,)
  # add colData info into meta.data
  obj@meta.data$count_label <- rownames(obj@meta.data)
  obj@meta.data <- merge(obj@meta.data, colData, by='count_label', all=TRUE)
  rownames(obj@meta.data) <- obj@meta.data$count_label
  # add library column
  obj@meta.data$library <- sapply(obj@meta.data$tag,
                                  function(x){
                                    return(strsplit(x, '-', fixed = T)[[1]][1])
                                  })
  # compute mito percent
  obj[['percent.mt']] <- PercentageFeatureSet(obj, pattern = "^MT-")
  # subset cells
  obj <- subset(obj, nCount_RNA>=65536 & nFeature_RNA>=6000 & percent.mt <=10)
  # plot vln
  vln_plot(obj, outdir)
  # remove plates with too little cells, N cells>=30
  plate_cells_df <- data.frame(table(obj@meta.data$tag))
  plates <- plate_cells_df[plate_cells_df$Freq>=min.good.wells,]$Var1
  obj <- subset(obj, tag%in%plates)
  # split by plates
  obj_list <- SplitObject(obj, split.by = split.by)
  
  if (is.null(filename)) {
    filename <- 'seurat.list.rds'
  }
  
  saveRDS(obj_list, 
          file = paste(outdir,'/', filename, sep=''))
}

### integrate seurat object list
integrate_seurat_list <- function(seurat_object_list,
                                  filename=NULL,
                                  outdir='./',
                                  normalization.method = 'LogNormalize',
                                  nfeatures=6000,
                                  k.filter=30) {
  message('Integrating seurat object list...')
  # check dir
  check_dir(outdir)
  check_file(seurat_object_list)
  # run
  obj_list <- readRDS(seurat_object_list)
  plan(multisession, workers=4)
  if (normalization.method == 'LogNormalize') {
    # preprocess
    obj_list <- lapply(X = obj_list, FUN = function(x) {
      x <- NormalizeData(x, verbose=FALSE)
      x <- FindVariableFeatures(x, 
                                selection.method = "dispersion",
                                nfeatures = nfeatures, 
                                verbose=FALSE)
    })
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = obj_list, 
                                          nfeatures = nfeatures,
                                          verbose=FALSE)
  } else {
    if (normalization.method == 'SCT') {
      obj_list <- lapply(X = obj_list, FUN = function(x){
        x <- SCTransform(x, 
                         vars.to.regress = "percent.mt", 
                         variable.features.n = nfeatures,
                         verbose = FALSE)
      })
      features <- SelectIntegrationFeatures(object.list = obj_list, 
                                            nfeatures = nfeatures)
      obj_list <- PrepSCTIntegration(object.list = obj_list, 
                                     anchor.features = features)
    }
  }
  obj.anchors <- FindIntegrationAnchors(object.list = obj_list, 
                                        normalization.method = normalization.method,
                                        anchor.features = features,
                                        k.filter = k.filter,
                                        verbose=FALSE)
  # this command creates an 'integrated' data assay
  obj.combined <- IntegrateData(anchorset = obj.anchors,
                                normalization.method = normalization.method,
                                k.weight = k.filter,
                                verbose=FALSE)
  if (normalization.method=='LogNormalize') {
    obj.combined <- ScaleData(obj.combined,
                              vars.to.regress = "percent.mt",
                              verbose=FALSE)
  }
  if (is.null(filename)){
    filename <- 'seurat.res'
  }
  
  SaveH5Seurat(obj.combined, 
               filename = paste(outdir, '/', filename,sep=''),
               overwrite = T,
               verbose=F)
}





## standard seurat analysis
standard_seurat_analysis <- function(seurat_data, 
                                     filename=NULL,
                                     outdir='./',
                                     dims=10,
                                     resolution=0.5) {
  message('Running standard seurat')
  # check dir
  check_dir(outdir)
  check_file(seurat_data)
  if (!file.exists(seurat_data)) {
    stop('ERROR: No seurat file found!')
  }
  obj <- LoadH5Seurat(seurat_data, verbose=F)
  ## RUN PCA
  obj <- RunPCA(obj, 
                verbose=FALSE)
  ## plot elbow to select the dims
  png(paste(outdir, '/ElbowPlot.png', sep=''),
      width = 1600,
      height = 800)
  elbow <- ElbowPlot(obj, reduction = 'pca')
  print(elbow)
  dev.off()
  
  ## run standard analysis
  obj <- RunUMAP(obj, 
                 reduction = "pca",
                 dims=1:dims,
                 verbose=FALSE)
  obj <- RunTSNE(obj,
                 reduction = "pca",
                 dims = 1:dims,
                 verbose=FALSE)
  obj <- FindNeighbors(obj,
                       reduction = "pca",
                       dims = 1:dims,
                       verbose=FALSE)
  obj <- FindClusters(obj,
                      resolution = resolution,
                      verbose=FALSE)
  ## save results
  if (is.null(filename)) {
    filename <- 'seurat.res'
  }
  SaveH5Seurat(obj, 
               filename = paste(outdir, '/', filename, sep=''),
               overwrite = T, verbose=F)
}

vln_plot <- function(obj, 
                     outdir,
                     width=30,
                     height=18,
                     units='cm') {
  obj@meta.data <- rename(obj@meta.data, 
                          'nFeature_RNA'='gene_count', 
                          'nCount_RNA'='UMI_count')
  # vln plot
  vln1 <- VlnPlot(obj,
                  features = c("gene_count", "UMI_count", "percent.mt"),
                  assay='RNA',
                  ncol = 3,
                  combine = TRUE,
                  pt.size = 0)
  for (i in 1:3) {
    vln1[[i]] <- vln1[[i]] +
      geom_boxplot(width=.2,col="black",fill="white")
  }
  vln1[[2]] <- vln1[[2]] + 
    scale_y_continuous(breaks=seq(0, max(obj@meta.data$UMI_count), 100000))
  
  ggsave(filename = paste(outdir, '/vln_all.png', sep=''),
         plot = vln1,
         width = width,
         height = height,
         units = units)
  
  vln2 <- VlnPlot(obj,
                  features = c("gene_count", "UMI_count"),
                  assay='RNA',
                  group.by = 'tag',
                  ncol = 2,
                  pt.size = 0)
  
  vln2[[1]] <- vln2[[1]] +
    geom_boxplot(width=.2,col="black",fill="white")
  vln2[[2]] <- vln2[[2]] + 
    geom_boxplot(width=.2,col="black",fill="white") +
    scale_y_continuous(breaks=seq(0, max(obj@meta.data$UMI_count), 100000)) +
    geom_hline(yintercept = 65536, color='red')
  
  ggsave(filename = paste(outdir, '/vln_plate.png', sep=''),
         plot = vln2,
         width = width,
         height = height,
         units = units)
}

seurat_plot <- function(seurat_data,
                        outdir='./',
                        width=30,
                        height=18,
                        units='cm') {
  message('Plotting seurat object...')
  check_dir(outdir)
  check_file(seurat_data)
  
  obj <- LoadH5Seurat(seurat_data, 
                      assays=c('integrated'),
                      verbose=F)
  
  # pca plot
  percentage<-round(obj@reductions$pca@stdev / 
                      sum(obj@reductions$pca@stdev) * 100,2)
  percentage<-paste(colnames(obj@reductions$pca@cell.embeddings),
                    "(", 
                    paste(as.character(percentage), "%", ")", sep=""))
  
  pca1 <- DimPlot(obj, reduction='pca', group.by = 'tag') + 
    xlab(percentage[1]) +
    ylab(percentage[2])
  pca2 <- DimPlot(obj, reduction='pca', group.by = 'seurat_clusters') + 
    xlab(percentage[1]) +
    ylab(percentage[2])
  pca3 <- DimPlot(obj, reduction='pca', group.by = 'detail_treatment') + 
    xlab(percentage[1]) +
    ylab(percentage[2])
  pca <- pca1+pca2+pca3
  ggsave(filename = paste(outdir, '/pca.png', sep=''),
         plot = pca,
         width = width,
         height = height,
         units = units)
  
  # umap plot
  p1 <- DimPlot(obj, reduction = "umap", group.by = "tag")
  p2 <- DimPlot(obj, reduction = "umap", group.by = 'seurat_clusters',
                label = TRUE, repel = TRUE)
  p3 <- DimPlot(obj, reduction = "umap", group.by = "detail_treatment")
  umap_plt <- p1+p2+p3
  ggsave(filename = paste(outdir, '/umap.png', sep=''),
         plot = umap_plt,
         width = width,
         height = height,
         units = units)
  
  # tsne plot
  p1 <- DimPlot(obj, reduction = "tsne", group.by = "tag")
  p2 <- DimPlot(obj, reduction = "tsne", group.by = 'seurat_clusters',
                label = TRUE, repel = TRUE)
  p3 <- DimPlot(obj, reduction = "tsne", group.by = "detail_treatment")
  tsne_plt <- p1+p2+p3
  ggsave(filename = paste(outdir, '/tsne.png', sep=''),
         plot = tsne_plt,
         width = width,
         height = height,
         units = units)
  
}

score_compound <- function(seurat_data,
                           outdir='./',
                           normalization.method = 'LogNormalize',
                           NC='DMSO',
                           PC='PC') {
  message('Scoring compounds...')
  check_dir(outdir)
  check_file(seurat_data)
  # set data
  seurat_obj <- LoadH5Seurat(seurat_data, verbose=F)
  if (normalization.method == 'LogNormalize') {
    normalized_data <- data.frame(seurat_obj[['RNA']]@data, check.names = F)
    normalized_data <- normalized_data[rownames(normalized_data)%in%rownames(seurat_obj[['integrated']]@scale.data),]
  } else {
    if (normalization.method == 'SCT') {
      normalized_data <- data.frame(seurat_obj[['SCT']]@data, check.names = F)
      normalized_data <- normalized_data[rownames(normalized_data)%in%rownames(seurat_obj[['SCT']]@scale.data),]
    }
  }
  
  
  
  meta.data <- seurat_obj@meta.data
  # compute nc/pc median
  nc_samples <- rownames(meta.data[meta.data$treatment=='NC',])
  nc_median <- apply(normalized_data[, nc_samples], 1, FUN = function(x) {
    return(median(x))
  })
  nc_median[nc_median==0] <- 1
  pc_samples <- rownames(meta.data[meta.data$treatment=='PC1',])
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
  
  max_s <- max(abs(score_df$score))
  
  max_s <- max(score_df$score)
  min_s <- min(score_df$score)
  
  score_df$normalized_score <- sapply(score_df$score, FUN = function(x) {
    x/max_s
    # ifelse(x>=0, x/max_s, -1*x/min_s)
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
  
}



### read in countdata and group data
# common parameters
min.cells = 4
min.features = 200
project = args$project
seurat_list_results = paste(args$project, '.seurat.rds', sep='')
seurat_results = paste(args$project, '.seurat.res', sep='')




if (!file.exists(paste(args$outdir,'/', seurat_results, '.h5seurat', sep=''))) {
  make_seurat_list(count_file = args$count_matrix,
                   group_info = args$group_info,
                   filename = seurat_list_results,
                   outdir = args$outdir)
  integrate_seurat_list(paste(args$outdir, seurat_list_results, sep='/'),
                        outdir=args$outdir,
                        normalization.method = args$normalization.method,
                        filename = seurat_results,
                        k.filter = as.integer(args$k_filter))
}

if (!is.null(args$dims)) {
  dims = as.integer(args$dims)
} else {
  dims <- user.input('Please select your dims (integer): ')
  dims <- as.integer(dims)
}

standard_seurat_analysis(paste(args$outdir,'/', seurat_results, '.h5seurat', sep=''), 
                         filename = seurat_results,
                         outdir = args$outdir,
                         dims=as.integer(dims))

seurat_plot(paste(args$outdir,'/', seurat_results, '.h5seurat', sep=''), 
            outdir = args$outdir)


score_compound(paste(args$outdir,'/', seurat_results, '.h5seurat', sep=''),
               outdir=args$outdir,
               normalization.method = args$normalization.method,
               NC = 'NC',
               PC='PC1')



