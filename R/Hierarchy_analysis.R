library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)

############################################################
## 1. Redistribution Mapping
############################################################

runRedistributionMapping <- function(
    object,
    da_results,
    cluster.layers
){
  fg.effect <-
    da_results$mc %>%
    tibble::rownames_to_column(
      "FG"
    ) %>%
    mutate(
      FG = as.character(FG),
      DA_score = sign(logFC) * (-log10(FDR + 1e-300))
    )

  meta.df <- data.frame(
    FG = as.character(object$metacell),
    cluster.layers,
    check.names = FALSE
  )

  calc.redistribution <- function(level.vec){
    fg.level <- data.frame(
      FG = unique(meta.df$FG),
      Parent = level.vec[match(unique(meta.df$FG), meta.df$FG)]
    )

    tmp <-
      left_join(
        fg.level,
        fg.effect,
        by="FG"
      )

    out <-
      lapply(
        unique(tmp$Parent),
        function(pp){
          tmp.pp <- tmp[tmp$Parent==pp,]
          effect.vec <- abs(tmp.pp$logFC) * (-log10(tmp.pp$FDR+1e-300))
          redistribution <-

            mean(
              abs(
                tmp.pp$logFC -
                  mean(
                    tmp.pp$logFC,
                    na.rm=TRUE
                  )
              ),
              na.rm=TRUE
            ) *

            sqrt(
              length(
                na.omit(
                  tmp.pp$logFC
                )
              )
            )

          data.frame(
            Parent = pp,
            nFG = nrow(tmp.pp),
            Redistribution = redistribution
          )
        }
      )

    bind_rows(out)
  }

  redistribution.df <-
    bind_rows(
      lapply(
        colnames(cluster.layers),
        function(L){
          x <- calc.redistribution(cluster.layers[,L])
          x$Level <- L
          x
        }
      )
    )


  parent.rank <-
    redistribution.df %>%
    arrange(
      desc(Redistribution)
    )


  ##########################################################
  ## Global IC scan
  ##########################################################

  IC.df <-
    as.data.frame(
      object@reductions$ica@cell.embeddings
    )

  IC.df$FG <-
    as.character(
      object$metacell
    )

  fg.program <-
    IC.df %>%
    group_by(FG) %>%
    summarise(
      across(
        everything(),
        mean
      )
    )


  fg.program <-
    left_join(
      fg.program,
      fg.effect,
      by="FG"
    )


  IC.names <-
    grep(
      "^IC",
      colnames(fg.program),
      value=TRUE
    )


  global.program.res <- list()

  for(ic in IC.names){
    fit <-
      lm(
        fg.program$logFC ~ fg.program[[ic]]
      )

    sm <- summary(fit)

    global.program.res[[ic]] <-
      data.frame(
        IC = ic,
        Beta = coef(fit)[2],
        P = sm$coefficients[2,4],
        R2 = sm$r.squared
      )
  }

  global.program.res <- bind_rows(global.program.res)

  global.program.res$FDR <-
    p.adjust(
      global.program.res$P,
      method="BH"
    )

  list(
    fg.effect = fg.effect,
    fg.program = fg.program,
    parent.rank = parent.rank,
    redistribution.df = redistribution.df,
    global.program.res = global.program.res
  )

}

############################################################
## 2. Global overview
############################################################

plotRedistributionOverview <- function(res){

  p1 <-
    ggplot(
      head(res$parent.rank, 30),
      aes(
        reorder(paste(Level,Parent),Redistribution),
        Redistribution
      )
    ) +
    geom_col() +
    coord_flip() +
    theme_bw()


  p2 <-
    ggplot(
      res$global.program.res,
      aes(Beta, -log10(FDR))
    ) +
    geom_point() +
    theme_bw()


  print(
    p1 + p2
  )

}

############################################################
## 3. Parent analysis
############################################################

analyzeParent <- function(
  object,
  res,
  cluster.layers,
  level,
  parent
){

  fg.effect <- res$fg.effect
  fg.program <- res$fg.program

  meta.df <- data.frame(
    FG = as.character(object$metacell),
    cluster.layers
  )

  fg.lookup <-
    meta.df %>%
    distinct(
      FG,
      .keep_all = TRUE
    )


  fg.top <- fg.lookup$FG[fg.lookup[[level]] == parent]

  fg.plot <-
    fg.effect %>%
    filter(FG %in% fg.top)

  fg.program.sub <-
    fg.program %>%
    filter(
      FG %in% fg.top
    )

  IC.names <- grep("^IC", colnames(fg.program.sub), value=TRUE)

  local.program <- list()

  for(ic in IC.names){
    fit <- lm(fg.program.sub$logFC ~ fg.program.sub[[ic]])

    sm <- summary(fit)
    local.program[[ic]] <-
      data.frame(
        IC = ic,
        Beta = coef(fit)[2],
        P = sm$coefficients[2,4],
        R2 = sm$r.squared
      )
  }

  local.program.res <-
    bind_rows(local.program)

  local.program.res$FDR <-
    p.adjust(
      local.program.res$P,
      method="BH"
    )

  local.program.res <-
    arrange(
      local.program.res,
      FDR
    )

  best.ic <- local.program.res$IC[1]

  ##########################################################
  ## UMAP/PCA landscape
  ##########################################################

  reds <- Reductions(object)

  if(any(grepl("umap",reds))){
    use.red <- grep("umap",reds, value=TRUE)[1]

  }else{
    if("pca.corrected" %in% reds){
      use.red <- "pca.corrected"
    }else{
      use.red <- "pca"
    }
  }

  coord <-
    as.data.frame(
      Embeddings(
        object,
        use.red
      )
    )

  coord <- coord[,1:2]
  colnames(coord) <- c("Dim1", "Dim2")

  coord$FG <- as.character(object$metacell)

  fg.coord <-
    coord %>%
    filter(
      FG %in% fg.top
    ) %>%
    group_by(FG) %>%
    summarise(
      Dim1 = mean(Dim1),
      Dim2 = mean(Dim2),
      CellNumber = n(),
      .groups="drop"

    )


  fg.coord <-
    left_join(
      fg.coord,
      fg.effect,
      by="FG"
    )

  fg.coord$ICscore <-
    fg.program.sub[[best.ic]][
      match(
        fg.coord$FG,
        fg.program.sub$FG
      )
    ]


  ##########################################################
  ## gene loading
  ##########################################################

  gene.loading <-
    object@reductions$pca_RNA@feature.loadings %*%
    object@reductions$ica@feature.loadings

  colnames(gene.loading) <- colnames(object@reductions$ica@cell.embeddings)


  gene.df <- data.frame(
      Gene = rownames(gene.loading),
      Loading = gene.loading[,best.ic]
    )

  gene.df <-
    gene.df %>%
    arrange(
      desc(
        abs(Loading)
      )
    )


  ##########################################################
  ## plots
  ##########################################################

  p1 <-
    ggplot(
      fg.plot,
      aes(
        reorder(FG,logFC),
        logFC,
        fill = logFC
      )
    ) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      na.value = "grey80"
    ) +
    labs(fill = "log2FC")


  p2 <-
    ggplot(
      local.program.res,
      aes(Beta, -log10(FDR))
    ) +
    geom_point() +
    theme_bw()


  p3 <-
    ggplot(
      fg.coord,
      aes(Dim1, Dim2,size=CellNumber,colour=logFC)
    ) +
    geom_point() +
    geom_text(aes(label=FG),size=3) +
    scale_colour_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "log2FC"
    ) + theme_bw()


  p4 <-
    ggplot(fg.coord,
      aes(Dim1, Dim2,size=CellNumber,colour=ICscore)
    ) +
    geom_point() +
    geom_text(aes(label=FG),size=3) +
    theme_bw()

  top.gene <-
    rbind(head(gene.df, 15), tail(gene.df,15))


  p5 <-
    ggplot(
      top.gene,
      aes(
        reorder(Gene, Loading),
        Loading,
        fill= Loading>0
      )
    ) +
    geom_col() +
    coord_flip() +
    theme_bw()


  ##########################################################
  ## enrichment
  ##########################################################

  if(requireNamespace("clusterProfiler", quietly=TRUE) &&
     requireNamespace("org.Hs.eg.db", quietly=TRUE)){

    suppressPackageStartupMessages(
      library(clusterProfiler)
    )

    suppressPackageStartupMessages(
      library(org.Hs.eg.db)
    )

    top.genes <- head(gene.df$Gene,200)

    ego <-
      enrichGO(
        gene = top.genes,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP"
      )
  }else{
    ego <- NULL
  }

  plotEnrichmentBubble <- function(
    ego,
    top.n = 15
  ){

    if(is.null(ego)){
      message("No enrichment result.")
      return(NULL)
    }

    ego.df <- as.data.frame(ego)

    if(nrow(ego.df)==0){
      message("No enriched pathways.")
      return(NULL)
    }

    ego.df <-
    ego.df %>%
      arrange(
        p.adjust
      ) %>%
      head(
        top.n
      )


    ego.df$Description <- factor(ego.df$Description, levels = rev(ego.df$Description))

    p <- ggplot(
      ego.df,
      aes(x = GeneRatio,
          y = Description,
          size =Count,
          colour = -log10(p.adjust))) +
      geom_point(alpha = 0.9) +
      scale_colour_gradient(
        low = "skyblue",
        high = "red"
      ) +
      theme_bw() +
      labs(
        x = "Gene Ratio",
        y =  NULL,
        colour = "-log10(FDR)",
        size = "Gene Count"
      )
    return(p)
  }


  p6 = plotEnrichmentBubble(ego = ego)

  print(
    (p1+p2)/
      (p3+p4)/
      (p5+p6)
  )

  return(
    list(
      best.ic = best.ic,
      local.program.res = local.program.res,
      fg.coord = fg.coord,
      gene.loading = gene.df,
      GO = ego
    )
  )
}



buildLeafOrder <- function(cluster.layers){

  current.order <-
    unique(cluster.layers[,1])

  for(level in 2:ncol(cluster.layers)){

    new.order <- c()
    for(parent in current.order){
      childs <-
        unique(
          cluster.layers[
            cluster.layers[,level-1] == parent,
            level
          ]
        )

      new.order <- c(new.order, childs)
    }
    current.order <- new.order
  }
  current.order
}


computeHierarchyLayout <- function(
    cluster.layers
){

  library(dplyr)

  ################################################
  ## leaf order
  ################################################

  leaf.order <-
    buildLeafOrder(
      cluster.layers
    )

  leaf.level <-
    ncol(cluster.layers)

  ################################################
  ## node table
  ################################################

  node.list <- list()

  for(i in seq_len(ncol(cluster.layers))){
    lv <- colnames(cluster.layers)[i]
    tmp <-
      data.frame(
        Node=
          unique(
            cluster.layers[,i]
          )
      )
    tmp$Level <- i
    tmp$LevelName <- lv
    node.list[[i]] <- tmp
  }

  node.df <-
    bind_rows(
      node.list
    )

  node.df$NodeID <-
    paste0(
      node.df$LevelName,
      "_",
      node.df$Node
    )

  ################################################
  ## edge table
  ################################################

  edge.list <- list()
  for(i in 2:ncol(cluster.layers)){
    parent.lv <- colnames(cluster.layers)[i-1]
    child.lv  <- colnames(cluster.layers)[i]
    tmp <-
      unique(
        data.frame(
          ParentID=
            paste0(
              parent.lv,
              "_",
              cluster.layers[,i-1]
            ),

          ChildID=
            paste0(
              child.lv,
              "_",
              cluster.layers[,i]
            )
        )
      )

    edge.list[[i]] <- tmp
  }

  edge.df <-
    bind_rows(
      edge.list
    )

  ################################################
  ## leaf x
  ################################################

  node.df$x <- NA_real_
  leaf.idx <-
    node.df$Level ==
    leaf.level

  node.df$x[leaf.idx] <-
    match(
      node.df$Node[leaf.idx],
      leaf.order
    )

  ################################################
  ## upward recursion
  ################################################

  for(level in seq(leaf.level-1,1,-1)){
    current.nodes <-
      node.df$NodeID[
        node.df$Level==level
      ]

    for(nn in current.nodes){
      child.ids <-
        edge.df$ChildID[
          edge.df$ParentID==nn
        ]

      child.x <-
        node.df$x[
          match(
            child.ids,
            node.df$NodeID
          )
        ]

      child.x <-
        child.x[
          !is.na(child.x)
        ]

      if(length(child.x)>0){
        node.df$x[
          node.df$NodeID==nn
        ] <-
          mean(child.x)
      }
    }
  }

  ################################################
  ## reorder inside each level
  ################################################

  node.df <-
    node.df %>%
    arrange(
      Level,
      x
    )

  ################################################
  ## y
  ################################################

  max.level <-
    max(
      node.df$Level
    )

  node.df$y <-
    max.level -
    node.df$Level + 1

  ################################################
  ## segments
  ################################################

  seg.df <-
    edge.df %>%
    left_join(
      node.df %>%
        dplyr::select(
          NodeID,
          x,
          y
        ),

      by=c(
        "ParentID"="NodeID"
      )
    ) %>%
    rename(
      x1=x,
      y1=y
    ) %>%
    left_join(
      node.df %>%
        dplyr::select(
          NodeID,
          x,
          y
        ),
      by=c(
        "ChildID"="NodeID"
      )
    ) %>%
    rename(
      x2=x,
      y2=y
    )

  ################################################
  ## sanity check
  ################################################

  cat(
    "Leaf order length:",
    length(leaf.order),
    "\n"
  )

  list(

    node.df =
      node.df,

    edge.df =
      edge.df,

    seg.df =
      seg.df,

    leaf.order =
      leaf.order

  )
}


plotRedistributionTree <- function(
    object,
    res,
    cluster.layers,
    group.by = NULL,
    max.level = NULL
){

  library(dplyr)
  library(ggplot2)
  library(scales)

  if(!is.null(max.level)){

    cluster.layers <-
      cluster.layers[
        ,
        seq_len(max.level),
        drop = FALSE
      ]
  }

  ################################################
  ## hierarchy layout
  ################################################

  layout.res <-
    computeHierarchyLayout(
      cluster.layers
    )

  node.df <- layout.res$node.df
  seg.df  <- layout.res$seg.df

  ################################################
  ## cell number
  ################################################

  node.df$CellNumber <- 0

  for(i in seq_len(ncol(cluster.layers))){

    tmp <-
      table(
        cluster.layers[,i]
      )

    idx <-
      node.df$Level == i

    node.df$CellNumber[idx] <-
      as.numeric(
        tmp[
          match(
            node.df$Node[idx],
            names(tmp)
          )
        ]
      )
  }

  ################################################
  ## Mean LogFC
  ################################################

  node.df$MeanLogFC <- 0

  if(!is.null(res$fg.effect)){

    fg.map <-
      data.frame(
        FG =
          unique(
            as.character(
              object$metacell
            )
          )
      )

    for(i in seq_len(ncol(cluster.layers))){

      fg.map[[colnames(cluster.layers)[i]]] <-

        cluster.layers[
          match(
            fg.map$FG,
            as.character(
              object$metacell
            )
          ),
          i
        ]
    }

    for(i in seq_len(ncol(cluster.layers))){

      lv <- colnames(cluster.layers)[i]

      tmp <-

        fg.map %>%

        left_join(
          res$fg.effect,
          by = "FG"
        ) %>%

        group_by(
          .data[[lv]]
        ) %>%

        summarise(
          MeanLogFC =
            mean(
              logFC,
              na.rm = TRUE
            ),
          .groups = "drop"
        )

      idx <-
        node.df$Level == i

      node.df$MeanLogFC[idx] <-

        tmp$MeanLogFC[
          match(
            node.df$Node[idx],
            tmp[[1]]
          )
        ]
    }
  }

  ################################################
  ## point size
  ################################################

  node.df$PointSize <-

    scales::rescale(
      log10(
        node.df$CellNumber + 1
      ),
      to = c(3,15)
    )

  ################################################
  ## CellType annotation
  ################################################

  leaf.df <- NULL

  if(
    !is.null(group.by) &&
    group.by %in% colnames(object@meta.data)
  ){

    leaf.level <-
      ncol(cluster.layers)

    leaf.nodes <-

      node.df %>%

      filter(
        Level == leaf.level
      )

    fg.type <-

      data.frame(

        FG =
          as.character(
            object$metacell
          ),

        CellType =
          object@meta.data[[group.by]]

      )

    fg.type <-

      fg.type %>%

      distinct(
        FG,
        .keep_all = TRUE
      )

    tmp <-

      data.frame(
        FG =
          unique(
            fg.type$FG
          )
      )

    tmp$LeafNode <-

      cluster.layers[
        match(
          tmp$FG,
          as.character(
            object$metacell
          )
        ),
        leaf.level
      ]

    tmp <-

      left_join(
        tmp,
        fg.type,
        by = "FG"
      )

    leaf.type <-

      tmp %>%

      group_by(
        LeafNode,
        CellType
      ) %>%

      summarise(
        n = n(),
        .groups = "drop"
      ) %>%

      group_by(
        LeafNode
      ) %>%

      slice_max(
        n,
        n = 1,
        with_ties = FALSE
      ) %>%

      ungroup()

    leaf.df <-

      leaf.nodes %>%

      mutate(

        CellType =

          leaf.type$CellType[
            match(
              Node,
              leaf.type$LeafNode
            )
          ]

      )
  }

  ################################################
  ## plot
  ################################################

  p <-

    ggplot() +

    ################################################
  ## edges
  ################################################

  geom_segment(

    data = seg.df,

    aes(
      x = x1,
      y = y1,
      xend = x2,
      yend = y2
    ),

    colour = "grey70",
    linewidth = 0.4

  ) +

  ################################################
  ## nodes
  ################################################

  geom_point(
    data = node.df,
    aes(
      x = x,
      y = y,
      size = log10(CellNumber+1),
      fill = MeanLogFC
    ),
    shape = 21,
    colour = "black",
    stroke = 0.4
  ) +
    scale_size_continuous(
      name = "Cell Number",
      # trans = "log10",
      range = c(2,8),
      breaks = log10(
        c(
          10,
          100,
          1000,
          10000
        ) + 1
      ),

      labels = c(
        "10",
        "100",
        "1000",
        "10000"
      )
    ) +

    scale_fill_gradient2(

      low  = "#2166AC",
      mid  = "white",
      high = "#B2182B",

      midpoint = 0,

      name = "Mean logFC"

    ) +

    ################################################
  ## labels
  ################################################

  geom_text(

    data = node.df,

    aes(
      x = x,
      y = y,
      label = Node
    ),

    size = 2.6

  ) +

    ################################################
  ## axes
  ################################################

  scale_y_continuous(

    breaks =

      rev(
        seq_len(
          ncol(cluster.layers)
        )
      ),

    labels =
      colnames(cluster.layers)

  ) +

    theme_bw() +

    theme(

      panel.grid =
        element_blank(),

      axis.title =
        element_blank()

    )

  ################################################
  ## CellType legend
  ################################################

  if(!is.null(leaf.df)){

    p <-

      p +

      ggnewscale::new_scale_fill() +

      geom_tile(

        data = leaf.df,

        aes(

          x = x,

          y =
            min(node.df$y) - 0.8,

          fill =
            CellType

        ),

        width = 0.8,

        height = 0.35

      ) +

      scale_fill_brewer(

        palette = "Set3",

        name = group.by

      )
  }

  print(p)

  invisible(

    list(

      plot =
        p,

      node.df =
        node.df,

      seg.df =
        seg.df,

      leaf.df =
        leaf.df

    )

  )
}

############################################################
## Usage
############################################################

# res <-
#   runRedistributionMapping(
#     object,
#     da_results,
#     cluster.layers
#   )
#
# plotRedistributionOverview(res)
# parent.res <-
#   analyzeParent(
#     object = object,
#     res = res,
#     cluster.layers = cluster.layers,
#     level = res$parent.rank$Level[26],
#     parent = res$parent.rank$Parent[26]
#   )
#
#
# tree.res <-
#   plotRedistributionTree(
#     object =object,
#     res =  res,
#     cluster.layers =cluster.layers[,1:dim(cluster.layers)[2]],
#     group.by = "CellType",
#     max.level =  12
#   )
#
# parent.res <-
#   analyzeParent(
#     object = object,
#     res = res,
#     cluster.layers = cluster.layers,
#     level = 'Lv.3',
#     parent = 4
#   )



