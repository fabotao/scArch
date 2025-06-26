#'
#' @title Identify metacells in single cell genomic data
#'
#' @description scArch aggregates cells into compact metacells based on cell condensation and recursive splitting with curated termination criteria.
#'
#' @usage FindMetacells(
#'   object,
#'   reduction = "pca",
#'   dims = 1:50,
#'   step = 20,
#'   min.cells = 10,
#'   seed = 1024,
#'   Q.final = 0.75,
#'   Q.sub = 0.6
#' )
#'
#' @param object A seurat object containing assays, reduction, meta data etc.
#' @param reduction The reduction slot for ICA analysis. Default is pca.
#' @param dims The target dimensions for ICA analysis. Default is 1:50.
#' @param steps The number of steps for diffusion condensation in the independent component (IC) space. The default value is set to 20.
#' @param min.cells The minimal number of cells in a cell cluster/subcluster. The default value is set to 10.
#' @param seed The seed for generating intial w.init for fastICA. Default is 1024.
#' @param Q.final The inner connectivity threshold for discarding metacells with less than 20 cells in the final refine step. The default value is set to 0.75.
#' @param Q.sub The inner connectivity threshold for evaluating over partition in the recursive step. The default value is set to 0.6.
#'
#'
#' @return Seurat object with metacell grouping information stored in the metacell of meta.data slot in object
#' @export
#' @examples
#' ## Starts from counts
#' library(scArch)
#' library(Seurat)
#' sc_object <- CreateSeuratObject(count=counts, project = "sc_object", min.cells = 3)
#' sc_object$percent.mt <- PercentageFeatureSet(sc_object, pattern = "^MT-")
#' sc_object <- subset(sc_object, percent.mt<20)
#' sc_object <- NormalizeData(sc_object)
#' sc_object <- FindVariableFeatures(sc_object, nfeatures=2000)
#' sc_object <- ScaleData(sc_object, do.center = F)
#' sc_object <- RunPCA(sc_object, npcs=50)
#' sc_object <- FindNeighbors(object = sc_object,
#'                            k.param = 20,
#'                            compute.SNN = F,
#'                            prune.SNN = 0,
#'                            reduction = "pca",
#'                            dims = 1:50,
#'                            force.recalc = F, return.neighbor = T)
#' sc_object <- FindMetacells(sc_object)
#'
FindMetacells <- function(object, reduction='pca', dims=1:50, steps=20, min.cells=10, seed=1024, Q.final=0.75, Q.sub=0.6){
  if(!reduction %in% names(object@reductions)){
    stop(paste0("The ", reduction," reduction slot does not exist."))
  }

  object <- RunICA(sob=object, dims=dims, reduction=reduction)

  k <- dim(object@neighbors$RNA.nn@nn.idx)[2]
  ## Build snn graph
  knn.100 <- list(object@neighbors$RNA.nn@nn.idx, object@neighbors$RNA.nn@nn.dist)
  names(knn.100) <- c('id', 'dist')
  knn.100[['k']] <- k
  knn.100[['sort']] <- T
  class(knn.100) <- 'kNN'
  snn.100 <- dbscan::sNN(knn.100, k=k, sort=T)
  snn.100$shared[,1] <- snn.100$shared[,2]
  object@graphs[['snn']] <- snn.100

  nmf.mat <- (object@reductions$ica@cell.embeddings)
  h.mat <- t(nmf.mat)
  for(i in 1:dim(h.mat)[1]){
    h.mat[i,] <- MeanShift1(V=nmf.mat[,i], M=object@graphs$snn$id[,1:k],N=(object@graphs$snn$shared[,1:k]), steps=steps)
  }

  object@reductions[['nmf_mod']] <- object@reductions$ica
  object@reductions$nmf_mod@cell.embeddings <- t(h.mat)

  cluster <- apply(h.mat, 2, function(x){
    pc.pos <- which.max(x)
    pc.limit <- abs(range(x))
    pc.neg <- which.min(x)
    if(pc.limit[1] > pc.limit[2]){
      return(-pc.neg)
    }else{
      return(pc.pos)
    }
  })

  #id.list <- tapply(1:dim(object)[2], cluster, function(x){return(x)})
  #cluster.adj <- .adjust.cluster(object=object, id.list) # Not significant with this step
  cluster.adj <- cluster
  uni.cluster <- sort(unique(cluster.adj))
  res.vec <- c()
  for(cl in uni.cluster){
    id <- which(cluster.adj == cl)
    id.cl <- paste0(id, '_', cl)
    res.vec <- c(res.vec, .div.mod(id.cl, sob=object, min.cells=min.cells, seed=seed, Q.sub = Q.sub))
  }

  id.cluster <- strsplit(res.vec, '_')
  id <- as.integer(sapply(id.cluster, '[', 1))
  subcluster.raw <- sapply(id.cluster, '[', 2)
  subcluster <- as.integer(factor(subcluster.raw[order(id)]))
  #sourceCpp('../Rpackage/Cage/src/knn_vote.cpp')
  id.list <- tapply(1:dim(object)[2], subcluster, function(x){x})
  cpt.subcluster <- sapply(id.list, function(x){.get.compact(sob = object, x)})

  count.mat <- matrix(0, dim(object)[2], length(id.list))
  for(m in 1:length(id.list)){
    m.id <- (id.list[[names(id.list)[m]]])
    cnt.m <- table(c(object@graphs$snn$id[m.id,1:min(c(length(m.id), 20)), drop=F]))
    count.mat[as.integer(names(cnt.m)),as.integer(m)] <- cnt.m
  }

  max.id <- apply(count.mat, 1, function(x){
    return(which.max(x))
  })

  max.rat <- tapply(max.id==subcluster, subcluster, function(x){sum(x)/length(x)})
  good.cl <- which(max.rat >= 0.5)
  good.id <- which(max.id==subcluster & subcluster %in% good.cl)
  subcluster0 <- subcluster

  change = T
  while(change){
    subcluster.adj <- knn_vote(V=subcluster0, M=object@graphs$snn$id[,], steps = 1)
    subcluster0.old <- subcluster0
    subcluster0 <- ifelse(1:dim(object)[2] %in% good.id, subcluster0, subcluster.adj)
    if(all(subcluster0 == subcluster0.old)){
      change=F
    }
  }

  id.list <- tapply(1:dim(object)[2], subcluster0, function(x){x})
  cnt <- sapply(id.list, length)
  cpt.subcluster <- sapply(id.list, function(x){.get.compact(sob = object, x)})
  small.cl <- union(names(cnt)[cnt < 5], intersect(names(cnt)[cnt < 20], names(cpt.subcluster)[cpt.subcluster < Q.final]))

  good.id <- which(!subcluster0 %in% as.integer(small.cl))
  change = T
  while(change){
    subcluster.adj <- knn_vote(V=subcluster0, M=object@graphs$snn$id[,], steps = 1)
    subcluster0.old <- subcluster0
    subcluster0 <- ifelse(1:dim(object)[2] %in% good.id, subcluster0, subcluster.adj)
    if(all(subcluster0 == subcluster0.old)){
      change=F
    }
  }

  object@meta.data['metacell'] <- subcluster0
  return(object)
}
