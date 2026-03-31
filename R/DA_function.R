


# object = readRDS('/home/rstudio/cell_type/GSE303823_AD_DLB_PDD/AD_obj.RDS')
#
# object = FindMetacells(object, reduction = 'pca', dims = 1:20)
#
#
# cluster.layers = get.layers(object = object)
# object$condition = as.character(object$group)
#
# condition_col = 'condition'
# sample_col = 'sample'



run_DA <- function(
    obj,
    condition_col,
    sample_col,
    batch_col=NULL,
    norm.method="TMM"
){

  ## Make design matrix
  design_df <- as_tibble(obj@meta.data[,c(sample_col, condition_col)])
  design_df <- distinct(design_df)
  # colnames(design_df)[1:2] = c('sample', 'condition')

  if (is.null(batch_col)) {
    design <- formula(paste('~', condition_col, collapse = ' '))
  } else {
    design <- formula(paste('~', batch_col, "+", condition_col, collapse = ' '))
  }

  condition_vec <- obj@meta.data[[condition_col]]
  sample_labels <- obj@meta.data[[sample_col]]

  clust.df <- data.frame("cell_id"=colnames(obj), "cluster"=obj$metacell)
  clust.df$Sample <- sample_labels
  clust.df$Condition <- condition_vec

  cluster.count <- table(clust.df$cluster, clust.df$Sample)
  attributes(cluster.count)$class <- "matrix"

  ## Test with same NB-GLM model as the Milo
  if(norm.method %in% c("TMM")){
    # message("Using TMM normalisation")
    dge <- DGEList(counts=cluster.count,
                   lib.size=colSums(cluster.count))
    dge <- calcNormFactors(dge, method="TMM")
  } else if(norm.method %in% c("logMS")){
    # message("Using logMS normalisation")
    dge <- DGEList(counts=cluster.count,
                   lib.size=colSums(cluster.count))
  }

  model <- model.matrix(design, data=design_df)
  rownames(model) <- design_df$synth_samples
  model <- model[colnames(cluster.count), ]

  dge <- estimateDisp(dge, model)
  fit <- glmQLFit(dge, model, robust=TRUE)

  n.coef <- ncol(model)
  DA.res <- as.data.frame(topTags(glmQLFTest(fit, coef=n.coef), sort.by='none', n=Inf))

  clust.df$logFC <- DA.res[match(clust.df$cluster,rownames(DA.res)), 'logFC']
  clust.df$FDR <- DA.res[match(clust.df$cluster,rownames(DA.res)), 'FDR']
  return(list(sc=clust.df, mc=DA.res))
}







run_hierarchical_DA <- function(
    obj,
    cluster.layers,
    condition_col,
    sample_col,
    batch_col=NULL,
    norm.method="TMM"
){

  ## Make design matrix
  design_df <- as_tibble(obj@meta.data[,c(sample_col, condition_col, batch_col)])
  design_df <- distinct(design_df)
  colnames(design_df)[1:2] = c('sample', 'condition')

  if (is.null(batch_col)) {
    design <- formula(paste('~', condition_col, collapse = ' '))
  } else {
    design <- formula(paste('~', batch_col, "+", condition_col, collapse = ' '))
  }

  condition_vec <- obj@meta.data[[condition_col]]
  sample_labels <- obj@meta.data[[sample_col]]

  logFC.matrix = matrix(NA, nrow=length(unique(cluster.layers[,dim(cluster.layers)[2]])), ncol=ncol(cluster.layers))
  dimnames(logFC.matrix) = list(sort(unique(cluster.layers[,dim(cluster.layers)[2]])),
                                colnames(cluster.layers))

  FDR.matrix = matrix(NA, nrow=length(unique(cluster.layers[,dim(cluster.layers)[2]])), ncol=ncol(cluster.layers))
  dimnames(FDR.matrix) = list(sort(unique(cluster.layers[,dim(cluster.layers)[2]])),
                              colnames(cluster.layers))

  clust.df <- data.frame("cell_id"=colnames(obj), "cluster"=cluster.layers[,dim(cluster.layers)[2]])
  clust.df$Sample <- sample_labels
  clust.df$Condition <- condition_vec

  for(kk in 1:dim(logFC.matrix)[2]){

    DA.layer = do.call(rbind, tapply(1:dim(cluster.layers)[1], cluster.layers[,kk], function(x){

      # print(cluster.layers[x[1],kk])

      clust.df.x = clust.df[x,]
      cluster.count <- table(clust.df.x$cluster, clust.df.x$Sample)
      attributes(cluster.count)$class <- "matrix"

      ## Test with same NB-GLM model as the Milo
      if(norm.method %in% c("TMM")){
        # message("Using TMM normalisation")
        dge <- DGEList(counts=cluster.count,
                       lib.size=colSums(cluster.count))
        dge <- calcNormFactors(dge, method="TMM")
      } else if(norm.method %in% c("logMS")){
       # message("Using logMS normalisation")
        dge <- DGEList(counts=cluster.count,
                       lib.size=colSums(cluster.count))
      }

      model <- model.matrix(design, data=design_df)
      rownames(model) <- design_df$sample
      model <- model[colnames(cluster.count), ]

      dge <- tryCatch({estimateDisp(dge, model)},
                      error = function(e) {
                        return(NULL)
                      })
      fit <- tryCatch({glmQLFit(dge, model, robust=TRUE)},
                      error = function(e) {
                        return(NULL)
                      })
      if(is.null(fit) | is.null(dge)){
        len = dim(cluster.count)[1]
        louvain.res = data.frame(logFC=rep(0, len),
                                 logCPM=rep(0,len),
                                 F=rep(0, len),
                                 PValue=rep(1,len),
                                 FDR=rep(1,len))
        rownames(louvain.res) = rownames(cluster.count)
        return(louvain.res)

      }else{
        n.coef <- ncol(model)
        louvain.res <- tryCatch({as.data.frame(topTags(glmQLFTest(fit, coef=n.coef), sort.by='none', n=Inf))},
                                error = function(e) {
                                  return(NULL)
                                })
        if(is.null(louvain.res)){
          len = dim(cluster.count)[1]
          louvain.res = data.frame(logFC=rep(0, len),
                                   logCPM=rep(0,len),
                                   F=rep(0, len),
                                   PValue=rep(1,len),
                                   FDR=rep(1,len))
          rownames(louvain.res) = rownames(cluster.count)
          return(louvain.res)
        }else{
          return(louvain.res)
        }
      }

    }))

    if(kk==dim(logFC.matrix)[2]){
      cluster.label = as.integer(rownames(DA.layer))
    }else{
      cluster.label = as.integer(sapply(strsplit(rownames(DA.layer), '.', fixed = T), '[', 2))
    }
    print(kk)
    print(rownames(DA.layer)[1:5])
    logFC.matrix[,kk] = DA.layer$logFC[order(cluster.label)]
    FDR.matrix[,kk] = DA.layer$FDR[order(cluster.label)]
  }
  FDR.matrix[is.na(FDR.matrix)] = 1

  return(list(logFC=logFC.matrix, FDR=FDR.matrix))
}


# aa = run_hierarchical_DA(obj = object, cluster.layers = cluster.layers, condition_col = 'condition', sample_col = 'sample')
#
#
# ump = data.frame(X=tapply(object@reductions$umap@cell.embeddings[,1], cluster.layers[,dim(cluster.layers)[2]], mean),
#                  Y=tapply(object@reductions$umap@cell.embeddings[,2], cluster.layers[,dim(cluster.layers)[2]], mean),
#                  size=tapply(1:dim(cluster.layers)[1], cluster.layers[,dim(cluster.layers)[2]], length))
#
# plot.layer.kk = function(kk){
#
#   ump.dat = ump
#   ump.dat$logFC = aa$logFC[,kk]
#   ump.dat$FDR = aa$FDR[,kk]
#   ump.dat$logFC[ump.dat$FDR > 0.1] = NA
#
#   ggplot(data=ump.dat) + geom_point(aes(x=X, y=Y, size=size, color=logFC)) + scale_color_gradient2(name='logFC', na.value = 'white')
# }
#
# meta = object@meta.data
# meta$X = object@reductions$umap@cell.embeddings[,1]
# meta$Y = object@reductions$umap@cell.embeddings[,2]
#
# ggplot(data=meta) + geom_point(aes(x=X, y=Y, color=factor(cell.type)))
#
#
#
# mc.type = tapply(object$cell.type, object$metacell, function(x){
#   cnt = table(x)
#   return(names(cnt)[which.max(cnt)])
# })
#
#
# which(aa$FDR[,10] < 0.1 & mc.type=='Inhibitory neuron')
#
# object$IN_grp1 = ifelse(object$metacell %in% c(588, 593, 594), 'grp1',
#                         ifelse(object$cell.type == 'Inhibitory neuron', 'grp2', 'oth'))
# Idents(object) = 'IN_grp1'
#
# DefaultAssay(object) = 'RNA'
# object = NormalizeData(object)
# mk = FindMarkers(object, ident.1 = 'grp2', ident.2 = 'grp1', slot='data')
# mk[mk$avg_log2FC < 0,]
# VlnPlot(object = object, features = c('MSR1','PDE1A','CDH12','SOX6','GRIK1','GRIK3','CDH8','ERBB4','CSMD3','GABRB2','GABRG3'),assay = 'RNA', stack=T, flip=T,slot='data')
#



