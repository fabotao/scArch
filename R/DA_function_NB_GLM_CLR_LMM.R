run_hierarchical_DA <- function(
    obj,
    cluster.layers,
    condition_col,
    sample_col,
    batch_col = NULL,
    norm.method = "TMM",
    ref.level = NULL,
    min_cells = 10,
    min_samples = 3,
    pseudo.count = 0.5
){

  ## Make design matrix
  design_df <- as_tibble(obj@meta.data[, c(sample_col, condition_col, batch_col)])
  design_df <- distinct(design_df)
  colnames(design_df)[1:2] <- c('sample', 'condition')

  ## 设置参照水平
  design_df$condition <- factor(design_df$condition)
  if(!is.null(ref.level)){
    design_df$condition <- relevel(design_df$condition, ref = ref.level)
  }

  ## 定义公式
  if(is.null(batch_col)){
    design      <- formula('~ condition')
    lmm_formula <- formula('clr_value ~ condition')
  } else {
    design      <- formula(paste('~ condition +', batch_col))
    lmm_formula <- formula(paste('clr_value ~ condition +', batch_col))
  }

  condition_vec <- obj@meta.data[[condition_col]]
  sample_labels <- obj@meta.data[[sample_col]]

  ## 缺陷三：验证样本对齐
  meta_samples <- unique(design_df$sample)
  obj_samples  <- unique(sample_labels)
  missing_in_meta <- setdiff(obj_samples, meta_samples)
  missing_in_obj  <- setdiff(meta_samples, obj_samples)
  if(length(missing_in_meta) > 0){
    stop("Samples in obj not found in meta.data: ",
         paste(missing_in_meta, collapse = ", "))
  }
  if(length(missing_in_obj) > 0){
    warning("Samples in meta.data not found in obj, will be ignored: ",
            paste(missing_in_obj, collapse = ", "))
  }

  ## 初始化输出矩阵
  n_leaf <- length(unique(cluster.layers[, dim(cluster.layers)[2]]))
  n_layer <- ncol(cluster.layers)
  mat_dimnames <- list(
    sort(unique(cluster.layers[, dim(cluster.layers)[2]])),
    colnames(cluster.layers)
  )

  logFC.matrix  <- matrix(NA, nrow = n_leaf, ncol = n_layer, dimnames = mat_dimnames)
  FDR.matrix    <- matrix(NA, nrow = n_leaf, ncol = n_layer, dimnames = mat_dimnames)
  PValue.matrix <- matrix(NA, nrow = n_leaf, ncol = n_layer, dimnames = mat_dimnames)
  method.matrix <- matrix(NA, nrow = n_leaf, ncol = n_layer, dimnames = mat_dimnames)

  clust.df <- data.frame(
    "cell_id" = colnames(obj),
    "cluster" = cluster.layers[, dim(cluster.layers)[2]]
  )
  clust.df$Sample    <- sample_labels
  clust.df$Condition <- condition_vec

  for(kk in 1:n_layer){

    DA.layer <- do.call(rbind, tapply(
      1:dim(cluster.layers)[1],
      cluster.layers[, kk],
      function(x){

        clust.df.x    <- clust.df[x, ]
        cluster.count <- table(clust.df.x$cluster, clust.df.x$Sample)
        attributes(cluster.count)$class <- "matrix"

        ## 缺陷三：验证cluster.count列名和design_df样本的对齐
        unmatched <- setdiff(colnames(cluster.count), design_df$sample)
        if(length(unmatched) > 0){
          warning("Sample name mismatch, removing unmatched samples: ",
                  paste(unmatched, collapse = ", "))
          cluster.count <- cluster.count[
            , !colnames(cluster.count) %in% unmatched, drop = FALSE
          ]
        }

        ## 过滤全零样本和全零子节点
        valid_samples  <- colSums(cluster.count) > 0
        valid_clusters <- rowSums(cluster.count) > 0
        if(sum(valid_samples) < ncol(cluster.count)){
          warning("Removing ", sum(!valid_samples),
                  " samples with zero counts in this node")
          cluster.count <- cluster.count[, valid_samples, drop = FALSE]
        }
        cluster.count <- cluster.count[valid_clusters, , drop = FALSE]

        ## 节点跳过检查
        skip_node   <- FALSE
        skip_reason <- NULL

        if(sum(cluster.count) < min_cells){
          skip_node   <- TRUE
          skip_reason <- "too_few_cells"
        } else if(sum(colSums(cluster.count) > 0) < min_samples){
          skip_node   <- TRUE
          skip_reason <- "insufficient_sample_coverage"
        } else if(nrow(cluster.count) < 2){
          skip_node   <- TRUE
          skip_reason <- "single_child_node"
        }

        if(skip_node){
          warning("Skipping node: ", skip_reason)
          len <- nrow(cluster.count)
          res <- data.frame(
            logFC  = rep(NA, len),
            logCPM = rep(NA, len),
            F      = rep(NA, len),
            PValue = rep(NA, len),
            FDR    = rep(NA, len),
            method = rep(paste0("skipped:", skip_reason), len)
          )
          rownames(res) <- rownames(cluster.count)
          return(res)
        }

        ## 检查自由度
        n_coef <- ncol(model.matrix(design, data = design_df))
        if(ncol(cluster.count) <= n_coef){
          warning("Insufficient degrees of freedom for node, skipping")
          len <- nrow(cluster.count)
          res <- data.frame(
            logFC  = rep(NA, len),
            logCPM = rep(NA, len),
            F      = rep(NA, len),
            PValue = rep(NA, len),
            FDR    = rep(NA, len),
            method = rep("skipped:insufficient_df", len)
          )
          rownames(res) <- rownames(cluster.count)
          return(res)
        }

        # --------------------------------------------------------
        # 分支：NB-GLM 或 CLR+LMM
        # --------------------------------------------------------
        if(norm.method %in% c("TMM", "logMS")){

          ## NB-GLM路线
          actual_method <- norm.method

          if(norm.method == "TMM"){
            dge <- DGEList(
              counts   = cluster.count,
              lib.size = colSums(cluster.count)
            )
            if(nrow(cluster.count) < 8){
              warning("Too few clusters (n=", nrow(cluster.count),
                      "), switching to logMS")
              actual_method <- "NB_logMS"
            } else {
              dge <- calcNormFactors(dge, method = "TMM")
              actual_method <- "NB_TMM"
            }
          } else {
            dge <- DGEList(
              counts   = cluster.count,
              lib.size = colSums(cluster.count)
            )
            actual_method <- "NB_logMS"
          }

          model <- model.matrix(design, data = design_df)
          rownames(model) <- design_df$sample
          model <- model[colnames(cluster.count), , drop = FALSE]

          dge <- tryCatch(
            estimateDisp(dge, model),
            error = function(e) NULL
          )
          fit <- tryCatch(
            glmQLFit(dge, model, robust = TRUE),
            error = function(e) NULL
          )

          if(is.null(fit) | is.null(dge)){
            len <- nrow(cluster.count)
            res <- data.frame(
              logFC  = rep(NA, len),
              logCPM = rep(NA, len),
              F      = rep(NA, len),
              PValue = rep(NA, len),
              FDR    = rep(NA, len),
              method = rep("failed:NB_GLM", len)
            )
            rownames(res) <- rownames(cluster.count)
            return(res)
          }

          n.coef <- ncol(model)
          louvain.res <- tryCatch(
            as.data.frame(topTags(
              glmQLFTest(fit, coef = n.coef),
              sort.by = 'none', n = Inf
            )),
            error = function(e) NULL
          )

          if(is.null(louvain.res)){
            len <- nrow(cluster.count)
            res <- data.frame(
              logFC  = rep(NA, len),
              logCPM = rep(NA, len),
              F      = rep(NA, len),
              PValue = rep(NA, len),
              FDR    = rep(NA, len),
              method = rep("failed:NB_GLM", len)
            )
            rownames(res) <- rownames(cluster.count)
            return(res)
          }

          ## 缺陷一：健壮提取PValue列名
          pval_col <- grep("^[Pp].*[Vv]alue$",
                           colnames(louvain.res), value = TRUE)
          if(length(pval_col) == 0){
            warning("PValue column not found in topTags output, ",
                    "columns are: ", paste(colnames(louvain.res), collapse = ", "))
            louvain.res$PValue <- NA
          } else if(pval_col != "PValue"){
            ## 统一重命名为PValue
            louvain.res$PValue <- louvain.res[[pval_col]]
          }

          louvain.res$method <- actual_method
          return(louvain.res)

        } else if(norm.method == "CLR_LMM"){

          ## CLR+LMM路线
          res <- run_clr_lmm_node(
            cluster.count = cluster.count,
            design_df     = design_df,
            lmm_formula   = lmm_formula,
            pseudo.count  = pseudo.count
          )

          ## 层内BH校正，和NB-GLM的topTags保持一致
          valid <- !is.na(res$PValue)
          if(sum(valid) > 0){
            res$FDR[valid] <- p.adjust(
              res$PValue[valid],
              method = "BH"
            )
          }

          return(res)
        }
      }
    ))

    ## 缺陷二：cluster.label解析加验证
    if(kk == n_layer){
      cluster.label <- suppressWarnings(as.integer(rownames(DA.layer)))
    } else {
      cluster.label <- suppressWarnings(
        as.integer(
          sapply(strsplit(rownames(DA.layer), '.', fixed = TRUE), '[', 2)
        )
      )
    }

    if(any(is.na(cluster.label))){
      warning("Layer ", kk, ": cluster.label parsing produced NA values. ",
              "Check cluster naming format. ",
              "Problematic rownames: ",
              paste(rownames(DA.layer)[is.na(cluster.label)], collapse = ", "))
      ## 用行序号兜底，避免赋值错位
      cluster.label[is.na(cluster.label)] <- which(is.na(cluster.label))
    }

    print(kk)
    print(rownames(DA.layer)[1:5])

    ## 缺陷一：健壮提取PValue
    pval_col <- grep("^[Pp].*[Vv]alue$", colnames(DA.layer), value = TRUE)
    pval_col <- if(length(pval_col) > 0) pval_col[1] else NULL

    logFC.matrix[, kk]  <- DA.layer$logFC[order(cluster.label)]
    FDR.matrix[, kk]    <- DA.layer$FDR[order(cluster.label)]
    PValue.matrix[, kk] <- if(!is.null(pval_col)){
      DA.layer[[pval_col]][order(cluster.label)]
    } else {
      NA
    }
    method.matrix[, kk] <- DA.layer$method[order(cluster.label)]
  }

  ## 问题二：logFC单位标注
  logfc_note <- switch(
    norm.method,
    "TMM"     = "log2 fold change (NB-GLM with TMM normalisation)",
    "logMS"   = "log2 fold change (NB-GLM with library size normalisation)",
    "CLR_LMM" = "CLR regression estimate (not log2 fold change, unit is natural log ratio)"
  )

  return(list(
    logFC    = logFC.matrix,
    FDR      = FDR.matrix,
    PValue   = PValue.matrix,
    method   = method.matrix,
    logFC_note = logfc_note
  ))
}


run_clr_lmm_node <- function(cluster.count,
                             design_df,
                             lmm_formula,
                             pseudo.count = 0.5){

  library(lme4)
  library(lmerTest)

  ## 子节点数不足
  if(nrow(cluster.count) < 2){
    len <- nrow(cluster.count)
    res <- data.frame(
      logFC  = rep(NA, len),
      logCPM = rep(NA, len),
      F      = rep(NA, len),
      PValue = rep(NA, len),
      FDR    = rep(NA, len),
      method = rep("skipped:single_child_node", len)
    )
    rownames(res) <- rownames(cluster.count)
    return(res)
  }

  ## CLR变换
  count_pseudo <- cluster.count + pseudo.count
  geo_mean     <- apply(count_pseudo, 2,
                        function(x) exp(mean(log(x))))
  clr_matrix   <- sweep(log(count_pseudo), 2, log(geo_mean), "-")

  ## 对每个子节点做LMM
  results <- lapply(rownames(clr_matrix), function(node){

    node_data <- data.frame(
      clr_value = as.numeric(clr_matrix[node, ]),
      design_df[match(colnames(clr_matrix), design_df$sample), ]
    )

    ## 确保condition因子水平和主函数一致
    node_data$condition <- factor(
      node_data$condition,
      levels = levels(design_df$condition)
    )

    ## 拟合LMM
    actual_method <- "CLR_LMM"
    fit <- tryCatch({
      lmer(lmm_formula, data = node_data,
           REML = FALSE,
           control = lmerControl(
             optimizer = "bobyqa",
             optCtrl   = list(maxfun = 2e5)
           ))
    }, error = function(e){
      warning("LMM failed for node ", node,
              ", falling back to lm: ", e$message)
      NULL
    })

    ## 问题四：LMM退化为LM时记录
    if(is.null(fit)){
      actual_method <- "CLR_LM_fallback"
      fixed_formula <- as.formula(
        paste("clr_value ~",
              paste(attr(terms(lmm_formula), "term.labels"),
                    collapse = " + "))
      )
      fit <- tryCatch(
        lm(fixed_formula, data = node_data),
        error = function(e2) NULL
      )
    }

    if(is.null(fit)){
      return(data.frame(
        logFC  = NA,
        logCPM = NA,
        F      = NA,
        PValue = NA,
        FDR    = NA,
        method = "failed:CLR_LMM",
        row.names = node
      ))
    }

    ## 提取系数表
    coef_table <- tryCatch(
      coef(summary(fit)),
      error = function(e) NULL
    )

    if(is.null(coef_table)){
      return(data.frame(
        logFC  = NA,
        logCPM = NA,
        F      = NA,
        PValue = NA,
        FDR    = NA,
        method = "failed:CLR_LMM",
        row.names = node
      ))
    }

    ## 稳健匹配condition系数行
    condition_row <- grep("^condition", rownames(coef_table), value = TRUE)

    if(length(condition_row) == 0){
      warning("condition coefficient not found in node ", node,
              ". Available rows: ",
              paste(rownames(coef_table), collapse = ", "))
      return(data.frame(
        logFC  = NA,
        logCPM = NA,
        F      = NA,
        PValue = NA,
        FDR    = NA,
        method = "failed:condition_not_found",
        row.names = node
      ))
    }

    ## 缺陷一：健壮提取PValue列名
    pval_col <- grep("^[Pp].*[Vv]alue\\]?$|^Pr",
                     colnames(coef_table), value = TRUE)
    pval_col <- if(length(pval_col) > 0) pval_col[1] else NULL

    if(is.null(pval_col)){
      warning("PValue column not found in coef table for node ", node,
              ". Available columns: ",
              paste(colnames(coef_table), collapse = ", "))
      return(data.frame(
        logFC  = NA,
        logCPM = NA,
        F      = NA,
        PValue = NA,
        FDR    = NA,
        method = "failed:pvalue_not_found",
        row.names = node
      ))
    }

    data.frame(
      logFC  = coef_table[condition_row, "Estimate"],
      logCPM = NA,
      F      = coef_table[condition_row, "t value"]^2,
      PValue = coef_table[condition_row, pval_col],
      FDR    = NA,
      method = actual_method,
      row.names = node
    )
  })

  do.call(rbind, results)
}
