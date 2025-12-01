#'
#' @title Get cellular hierarchy based on the identified fine-grained cell clusters
#'
#' @description get.layers infers cellular hierarchy based on the simlarity of aggregated expression of identified fine-grained clusters.
#'
#' @usage get.layers(
#'   object
#' )
#'
#' @param object A seurat object processed by FindMetacells function.
#'
#'
#' @return A matrix containing multi-resolution clusters by columns, with each row representing a cell.
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
#' Layered.cluster = get.layers(object = sc_object)
#'
get.layers <- function(object){

  cell.clusters <- object$metacell
  cl.uniq <- sort(unique(cell.clusters))
  cell.num <- length(cell.clusters)
  pca.mat <- object@reductions$pca@cell.embeddings

  ave.mat.list <- tapply(1:(cell.num), cell.clusters, function(x){colMeans(pca.mat[x,,drop=F])})
  ave.mat <- (matrix(unlist(ave.mat.list), nrow=length(ave.mat.list), ncol=dim(pca.mat)[2], byrow = T))
  dimnames(ave.mat)[[1]] <- names(ave.mat.list)

  pmin.mat <- lsa::cosine(t(ave.mat))

  hc <- hclust(as.dist(1-(pmin.mat)), method='average')
  levels <- dim(pmin.mat)[1]

  sil.mat <- matrix(NA, levels, levels-1)
  cl.mat <- matrix(NA, levels, levels)
  for(k in 2:(levels-1)){
    ks <- cutree(hc, k=k)
    cl.mat[,k] = ks
    sil.mat[,k] <- cluster::silhouette(dist=1-pmin.mat, x=ks)[,3]
  }
  cl.mat[, levels] <- cutree(hc, k=levels)

  get.robust <- function(sil.mat){
    sil.ind <- matrix(0, nrow(sil.mat), ncol(sil.mat))
    for(kk in 1:dim(sil.mat)[1]){
      sil.kk <- sil.mat[kk,]
      dif.kk.3 <- which.max(sil.kk)
      ks <- cutree(hc, k=dif.kk.3)
      for(jj in dif.kk.3){
        sil.ind[which(ks==ks[kk]),jj] <- 1
      }
    }
    return(sil.ind)
  }

  sil.ind <- get.robust(sil.mat = sil.mat)

  sil.sum <- colSums(sil.mat)[2:dim(sil.mat)[2]]
  diff.sil <- unlist(lapply(2:(length(sil.sum)-1), function(x){
    return(sil.sum[x]*2 - sum(sil.sum[c(x-1, x+1)]))
  }))
  peak.ind <- union(which(diff.sil > sd(diff.sil)) + 2, which.max(sil.sum))

  cls.mat <- matrix(0, nrow = nrow(sil.ind), ncol=ncol(sil.ind))
  active.k = 1
  for(m in 1:dim(cls.mat)[2]){
    sil.m = sil.ind[,m]
    current.col = cls.mat[,active.k]
    if(all(sil.m==0) & !m %in% peak.ind){
      next
    }else{
      if(m %in% peak.ind){
        if(all(current.col==0)){
          cls.mat[,active.k] = cl.mat[,m]
          active.k = active.k + 1
        }else{
          cls.mat[,active.k+1] = cl.mat[,m]
          active.k = active.k + 2
        }
      }else{
        id.m <- which(sil.m==1)
        if(all(cls.mat[id.m,active.k]==0)){
          cls.mat[id.m,active.k] = m
        }else{
          active.k = active.k + 1
          cls.mat[id.m,active.k] = m
        }
      }
    }
  }
  N = dim(cls.mat)[1]
  for(zz in 1:dim(cls.mat)[2]){
    zero.zz <- which(cls.mat[,zz]==0)
    if(length(zero.zz)==N | length(zero.zz)==0) next
    levels = setdiff(unique(cls.mat[,zz]), 0)
    lens <- c()
    for(lev in levels){
      lens <- c(lens, length(table(cl.mat[zero.zz,lev])))
    }
    cand.level <- levels[which.max(lens)]
    cnt <- table(cl.mat[zero.zz,cand.level])
    for(vv in 1:length(cnt)){
      cls.mat[intersect(zero.zz, which(cl.mat[,cand.level]==as.integer(names(cnt)[vv]))),zz] <- max(levels) + vv
    }
  }

  cls.mat <- cls.mat[,colSums(cls.mat) > 0]
  cls.mat <- cbind(cls.mat, cl.mat[,dim(cl.mat)[2]])

  cluster.mat <- matrix(NA, length(cell.clusters), dim(cls.mat)[2])
  for(kk in 1:dim(cluster.mat)[2]){
    ks <- cls.mat[,kk]
    cluster.mat[,kk] <- ks[match(cell.clusters, as.integer(dimnames(pmin.mat)[[1]]))]
  }
  dimnames(cluster.mat)[[2]] = paste0('Lv.', 1:dim(cluster.mat)[2])
  return((cluster.mat))

}

