# scArch
scArch (R package) is an supervised and graph-based algorithm for metacell identification, that demonstrates superior performance to existing algorithms.


## Installation
```R
  # Install in R with devtools
  library(devtools)
  install_github('fabotao/scArch')
```

## Usage
Starts from count matrix that is proprocessed by Seurat package

```R
  library(scArch)
  library(Seurat) 
  
  # counts as input matrix for preprocessing using Seurat
  sc_object <- CreateSeuratObject(count=counts, project = "sc_object", min.cells = 3)
  sc_object$percent.mt <- PercentageFeatureSet(sc_object, pattern = "^MT-")
  sc_object <- subset(sc_object, percent.mt<20)
  sc_object <- NormalizeData(sc_object)
  sc_object <- FindVariableFeatures(sc_object, nfeatures=2000)
  sc_object <- ScaleData(sc_object, do.center = F)
  sc_object <- RunPCA(sc_object, npcs=50)
  sc_object <- FindNeighbors(object = sc_object,
                             k.param = 20, 
                             compute.SNN = F, 
                             prune.SNN = 0, 
                             reduction = "pca", 
                             dims = 1:50, 
                             force.recalc = F, 
                             return.neighbor = T)
  
  # Use scArch to derive metacells  
  sc_object <- FindMetacells(sc_object)
```

## Citation


