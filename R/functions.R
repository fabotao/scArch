
## Dimension reduction with NMF
RunICA <- function(sob,
                   n.comp=50, reduction='pca', seed=1024){
  # sob2 = sob
  # if(min(SeuratObject::GetAssayData(sob2,slot = 'scale.data'))<0){
  #   sob2 = Seurat::ScaleData(sob2,do.center = F)
  # }
  set.seed(seed)
  w.init <- matrix(rnorm(n.comp^2),n.comp,n.comp)
  vm <- sob@reductions[[reduction]]@cell.embeddings
  res <- fastICA::fastICA((vm), n.comp=n.comp, verbose = F, w.init = w.init)

  ## Seurat
  cell.embeddings = res$S
  feature.loadings = res$A
  assay = 'RNA'
  reduction.key <- "IC_"
  #rownames(x = feature.loadings) <- rownames(x = vm)
  #colnames(x = feature.loadings) <- paste0(reduction.key, 1:npcs)
  rownames(x = cell.embeddings) <- rownames(x = vm)
  colnames(x = cell.embeddings) <- paste0(reduction.key, 1:dim(cell.embeddings)[2])

  reduction.data <- SeuratObject::CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    stdev = 0,
    key = reduction.key,
    misc = list(tol = res$tol, iter=res$iter)
  )

  sob[['ica']] <- reduction.data

  # data("pbmc_small")
  # nmfReu = pbmc_small@reductions$pca
  # nmfReu@cell.embeddings <- t(res$h * res$d)
  # nmfReu@feature.loadings <-  res$w
  # nmfReu@key <- "Factor_"
  # nmfReu@stdev <- res$d
  # nmfReu@assay.used <- SeuratObject::DefaultAssay(sob)
  # sob@reductions$nmf <- sob@reductions$pca
  # sob@reductions$nmf@cell.embeddings <- t(res@h)
  # sob@reductions$nmf@feature.loadings <- res@w
  # sob@reductions$nmf@misc = res@misc
  # sob@reductions$nmf@key = "Factor_"
  # sob@reductions$nmf@stdev = res@d
  # sob@reductions[reduction.name]=nmfReu
  # dimnames(sob@reductions$nmf@cell.embeddings) <- list(dimnames(vm)[[2]], paste0('Factor_', 1:npcs))
  sob
}



## working horse for subcluster identification
.div.mod <- function(id.mod, sob, min.cells=10, steps=20, seed=1024, Q.sub=0.6){

  k <- sob@graphs$snn$k
  #exp.mat <- sob@assays$RNA@scale.data
  h.mat <- sob@reductions$nmf_mod@cell.embeddings
  input.list <- strsplit(id.mod, '_')
  id <- as.integer(sapply(input.list, '[', 1))
  mod <- unique(sapply(input.list, '[', 2))
  if(grepl('+', mod, fixed = T)){return(id.mod)}

  if(length(id) >= min.cells){
    npcs <- 10
    sub.mat <- h.mat[id,]
    set.seed(seed)
    nm <- try(
          fastICA::fastICA(sub.mat, n.comp = 10, verbose = F),
          silent=T)

    if(length(nm)==1){
      return(id.mod)
    }

    #nm <- fastICA::fastICA(sub.mat, n.comp = 10, verbose = F)
    subcluster <- apply(nm$S, 1, function(x){
      pc.pos <- which.max(x)
      pc.limit <- abs(range(x))
      pc.neg <- which.min(x)
      if(pc.limit[1] > pc.limit[2]){
        return(-pc.neg)
      }else{
        return(pc.pos)
      }
    })
    if(length(unique(subcluster)) == 1){
      return(id.mod)
    }else{
      subcluster0 <- subcluster
      id.list <- tapply(id, subcluster0, function(x){return(x)})
      cpt.sub0 <- sapply(id.list, function(x){.get.compact(sob, x)})
      cpt.less <- cpt.sub0[cpt.sub0 <= 0.5]
      if(length(cpt.less) < 2){
        subcluster.comb <- subcluster0
      }else{
        pair <- combn(names(cpt.less), 2)
        cpt.comb <- apply(pair, 2, function(x){.get.compact(sob, id[subcluster %in% x])})
        cpt1 <- cpt.sub0[pair[1,]]
        cpt2 <- cpt.sub0[pair[2,]]
        ind <- which(cpt.comb > cpt1 & cpt.comb > cpt2)
        if(length(ind)==0){
          subcluster.comb <- subcluster0
        }else{
          pair.comb <- pair[,ind,drop=F]
          check <- apply(pair.comb, 2, function(x){
            id.2 <- id[subcluster0==x[2]]
            id.1 <- id[subcluster0==x[1]]
            cpt.samp2 <- .get.compact(sob, c(id.1, sample(id.2, round(length(id.2)/2))))  ## 是否随机抽样，还是按照shared 数量排序
            cpt.samp1 <- .get.compact(sob, c(id.2, sample(id.1, round(length(id.1)/2))))
            if(cpt.samp1 > cpt.sub0[x[2]] & cpt.samp2 > cpt.sub0[x[1]]){
              return(T)
            }else{
              return(F)
            }
          })
          pair.valid <- pair.comb[,check,drop=F]
          if(sum(check)==0){
            subcluster.comb <- subcluster0
          }else{
            subcluster.comb <- subcluster0
            for(i in 1:dim(pair.valid)[2]){
              cl.i <- unique(subcluster.comb[subcluster0 %in% pair.valid[,i]])
              subcluster.comb[subcluster.comb %in% cl.i] <- paste0(unique(subcluster.comb[subcluster.comb %in% cl.i]), collapse = '+')
            }
          }
        }
      }
      id.list <- tapply(id, subcluster.comb, function(x){return(x)})
      cpt.sub <- sapply(id.list, function(x){.get.compact(sob, x)})
      len.id <- sapply(id.list, function(x){length(x)})
      if(all(cpt.sub==0) | length(unique(subcluster.comb))==1){return(id.mod)}
      bad.cl <- names(len.id)[len.id < 20 & cpt.sub <= Q.sub]
      if(length(bad.cl) == length(len.id)){
        return(id.mod)
      }
      cnt <- table(subcluster.comb)
      cnt.bad <- sort(cnt[bad.cl])
      for(cl in names(cnt.bad)){
        id.cl <- id[subcluster.comb==cl]
        id.cl.neib <- unique(unlist(sob@graphs$snn$id[id.cl, 1:min(c(20, length(id.cl)))]))
        share.len <- sapply(id.list, function(x){length(intersect(x, id.cl.neib))})
        id.neib <- names(which.max(share.len[names(share.len) != cl]))
        id.neib.now <- unique(subcluster.comb[id %in% id.list[[id.neib]][1]])
        subcluster.comb[which(id %in% id.cl)] <- id.neib.now
      }

      if(length(unique(subcluster.comb))==1){return(id.mod)}

      id.mod.sub <- paste0(id.mod, '.', subcluster.comb)
      id.mod.list <- tapply(id.mod.sub, subcluster.comb, function(x){return(x)})
      return(unlist(sapply(id.mod.list, function(x){
        .div.mod(x, sob)
      })))

    }
  }else{
    return(id.mod)
  }
}






.get.compact <- function(sob, id){
  if(length(id)<5){
    return(0)
  }else{
    len <- min(length(id), 20)
    snn.id <- c(sob@graphs$snn$id[id,1:len])
    return(length(snn.id[snn.id %in% id])/(len * length(id)))
  }
}





.find.mins <- function (x, thresh = 0)
{
  pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) > 0) + 2
  if (!missing(thresh)) {
    pks[x[pks - 1] - x[pks] > thresh]
  }
  else pks
}





metacell.expr <- function(sob){
  cell.clusters <- sob$subcluster
  cl.uniq <- sort(unique(cell.clusters))
  cell.num <- length(cell.clusters)

  exp.2000 <- sob@assays$RNA@counts

  ave.mat.list <- tapply(1:(cell.num), cell.clusters, function(x){Matrix::rowSums(exp.2000[,x,drop=F])})
  ave.mat <- (matrix(unlist(ave.mat.list), nrow=length(ave.mat.list), ncol=dim(exp.2000)[1], byrow = T))
  dimnames(ave.mat)[[1]] <- names(ave.mat.list)
  dimnames(ave.mat)[[2]] <- dimnames(exp.2000)[[1]]
  return(t(ave.mat))
}





