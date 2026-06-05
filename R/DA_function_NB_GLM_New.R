run_hierarchical_DA <- function(
    obj,
    cluster.layers,
    condition_col,
    sample_col,
    batch_col       = NULL,
    norm.method     = c("TMM", "RLE", "logMS"),
    ref.level       = NULL,
    min_cells       = 10,
    min_samples     = 3,
    max_layer       = NULL,
    model.contrasts = NULL,
    robust          = TRUE,
    use.glmm        = FALSE,
    REML            = TRUE,
    glmm.solver     = "Fisher",
    max.iters       = 50,
    max.tol         = 1e-5,
    force           = FALSE
){

  ## -------------------------------------------------------
  ## 第一步：验证cluster.layers和obj的对齐关系
  ## -------------------------------------------------------
  if(nrow(cluster.layers) != ncol(obj)){
    stop("cluster.layers has ", nrow(cluster.layers),
         " rows but obj has ", ncol(obj),
         " cells. They must match exactly.")
  }

  if(!is.null(rownames(cluster.layers))){
    cell_ids  <- colnames(obj)
    layer_ids <- rownames(cluster.layers)

    if(!identical(cell_ids, layer_ids)){
      common_ids <- intersect(cell_ids, layer_ids)

      if(length(common_ids) < length(cell_ids)){
        missing <- setdiff(cell_ids, layer_ids)
        stop("cluster.layers missing ", length(missing),
             " cells. First few: ",
             paste(head(missing, 5), collapse = ", "))
      }

      warning("cluster.layers row order differs from obj. Reordering.")
      cluster.layers <- cluster.layers[cell_ids, , drop = FALSE]
    }
  } else {
    warning("cluster.layers has no rownames. Assuming row order matches ",
            "colnames(obj). Set rownames(cluster.layers) <- colnames(obj) ",
            "to ensure correctness.")
    rownames(cluster.layers) <- colnames(obj)
  }

  ## -------------------------------------------------------
  ## 确保cluster.layers有有意义的列名
  ## -------------------------------------------------------
  if(is.null(colnames(cluster.layers))){
    colnames(cluster.layers) <- paste0("layer_", seq_len(ncol(cluster.layers)))
    warning("cluster.layers has no colnames. Auto-assigned: ",
            paste(colnames(cluster.layers), collapse = ", "))
  }

  ## -------------------------------------------------------
  ## 行列定义分离
  ## 行：原始最高分辨率全部叶节点（用cluster标签作为行名）
  ## 列：受max_layer限制的检验层次（用层次标签作为列名）
  ## -------------------------------------------------------
  n_layers_full <- ncol(cluster.layers)

  leaf_labels_full <- cluster.layers[, n_layers_full]

  ## 关键修复：行名使用最高分辨率层的cluster标签
  ## sort保证行名顺序一致，as.character确保是字符串
  all_leaf_nodes <- sort(unique(as.integer(leaf_labels_full)))
  all_leaf_names <- paste0("cluster_", all_leaf_nodes)
  n_leaf_full    <- length(all_leaf_nodes)

  ## 建立cluster整数值到行名的映射
  leaf_name_map <- setNames(all_leaf_names, as.character(all_leaf_nodes))

  if(!is.null(max_layer)){
    n_layers_test <- min(max_layer, n_layers_full)
    message("DA testing limited to layer 1:", n_layers_test,
            " (out of ", n_layers_full, " total layers)")
  } else {
    n_layers_test <- n_layers_full
  }

  cluster.layers.test <- cluster.layers[, 1:n_layers_test, drop = FALSE]

  ## 每层的cluster数量（用于message）
  n_clusters_per_layer <- apply(
    cluster.layers.test, 2,
    function(x) length(unique(x))
  )
  message("cluster.layers validated: ",
          nrow(cluster.layers), " cells x ",
          n_layers_full, " layers | clusters per layer: ",
          paste(apply(cluster.layers, 2,
                      function(x) length(unique(x))),
                collapse = " -> "))

  ## normalization方法验证
  norm.method <- norm.method[1]
  if(!norm.method %in% c("TMM", "RLE", "logMS")){
    stop("norm.method must be one of: TMM, RLE, logMS")
  }

  ## Make design matrix
  design_df <- as_tibble(
    obj@meta.data[, c(sample_col, condition_col, batch_col)]
  )
  design_df <- distinct(design_df)
  colnames(design_df)[1:2] <- c('sample', 'condition')

  design_df$condition <- factor(design_df$condition)
  if(!is.null(ref.level)){
    design_df$condition <- relevel(design_df$condition, ref = ref.level)
  }

  if(is.null(batch_col)){
    design <- formula('~ condition')
  } else {
    if(use.glmm){
      design <- formula(paste('~ condition + (1|', batch_col, ')'))
    } else {
      design <- formula(paste('~ condition +', batch_col))
    }
  }

  ## GLMM样本量检查
  n_samples_design <- nrow(design_df)
  if(use.glmm && n_samples_design < 60 && !force){
    stop("Attempting to use GLMM with N=", n_samples_design,
         ". Strongly discouraged for N<60. Set force=TRUE to override.")
  } else if(use.glmm && n_samples_design < 60 && force){
    warning("Running GLMM with small sample size N=", n_samples_design,
            ". Estimates may not be reliable.")
  }

  condition_vec <- obj@meta.data[[condition_col]]
  sample_labels <- obj@meta.data[[sample_col]]

  ## 样本对齐验证
  meta_samples    <- unique(design_df$sample)
  obj_samples     <- unique(sample_labels)
  missing_in_meta <- setdiff(obj_samples, meta_samples)
  missing_in_obj  <- setdiff(meta_samples, obj_samples)
  if(length(missing_in_meta) > 0){
    stop("Samples in obj not found in meta.data: ",
         paste(missing_in_meta, collapse = ", "))
  }
  if(length(missing_in_obj) > 0){
    warning("Samples in meta.data not found in obj: ",
            paste(missing_in_obj, collapse = ", "))
  }

  ## -------------------------------------------------------
  ## 初始化输出矩阵
  ## 行名：cluster_X（最高分辨率层的cluster标签）
  ## 列名：层次标签（来自cluster.layers的列名）
  ## -------------------------------------------------------
  mat_dimnames <- list(
    all_leaf_names,                          # 行名：cluster_0, cluster_1, ...
    colnames(cluster.layers.test)            # 列名：layer_1, layer_2, ...
  )

  logFC.matrix  <- matrix(NA, nrow = n_leaf_full, ncol = n_layers_test,
                          dimnames = mat_dimnames)
  FDR.matrix    <- matrix(NA, nrow = n_leaf_full, ncol = n_layers_test,
                          dimnames = mat_dimnames)
  PValue.matrix <- matrix(NA, nrow = n_leaf_full, ncol = n_layers_test,
                          dimnames = mat_dimnames)
  method.matrix <- matrix(NA, nrow = n_leaf_full, ncol = n_layers_test,
                          dimnames = mat_dimnames)

  ## clust.df使用带前缀的cluster名称
  clust.df <- data.frame(
    "cell_id"    = colnames(obj),
    "cluster"    = leaf_name_map[as.character(leaf_labels_full)],
    stringsAsFactors = FALSE
  )
  clust.df$Sample    <- sample_labels
  clust.df$Condition <- condition_vec

  ## -------------------------------------------------------
  ## 逐层检验
  ## -------------------------------------------------------
  for(kk in 1:n_layers_test){

    parent_labels_kk <- cluster.layers.test[, kk]

    DA.layer.list <- tapply(
      seq_len(nrow(cluster.layers.test)),
      parent_labels_kk,
      function(leaf_idx){

        ## 该父节点下的叶节点（用带前缀的名称）
        leaf_nodes_raw       <- leaf_labels_full[leaf_idx]
        leaf_nodes_in_parent <- leaf_name_map[as.character(leaf_nodes_raw)]

        clust.df.x <- clust.df[
          clust.df$cluster %in% leaf_nodes_in_parent, ]

        cluster.count <- table(clust.df.x$cluster, clust.df.x$Sample)
        attributes(cluster.count)$class <- "matrix"

        ## 样本名对齐验证
        unmatched <- setdiff(colnames(cluster.count), design_df$sample)
        if(length(unmatched) > 0){
          warning("Sample name mismatch, removing: ",
                  paste(unmatched, collapse = ", "))
          cluster.count <- cluster.count[
            , !colnames(cluster.count) %in% unmatched, drop = FALSE]
        }

        ## 过滤全零样本和全零叶节点
        valid_samples  <- colSums(cluster.count) > 0
        valid_clusters <- rowSums(cluster.count) > 0
        if(sum(valid_samples) < ncol(cluster.count)){
          cluster.count <- cluster.count[, valid_samples, drop = FALSE]
        }
        cluster.count <- cluster.count[valid_clusters, , drop = FALSE]

        all_leaves_in_parent  <- as.character(
          sort(unique(leaf_nodes_in_parent))
        )
        valid_leaves_in_count <- rownames(cluster.count)
        filtered_leaves       <- setdiff(all_leaves_in_parent,
                                         valid_leaves_in_count)

        ## 辅助函数：构建跳过结果
        make_skip_res <- function(nodes, reason){
          len <- length(nodes)
          res <- data.frame(
            logFC  = rep(NA, len),
            logCPM = rep(NA, len),
            F      = rep(NA, len),
            PValue = rep(NA, len),
            FDR    = rep(NA, len),
            method = rep(paste0("skipped:", reason), len)
          )
          rownames(res) <- as.character(nodes)
          res
        }

        ## 节点跳过检查
        if(sum(cluster.count) < min_cells){
          return(make_skip_res(all_leaves_in_parent, "too_few_cells"))
        }
        if(sum(colSums(cluster.count) > 0) < min_samples){
          return(make_skip_res(all_leaves_in_parent,
                               "insufficient_sample_coverage"))
        }
        if(nrow(cluster.count) < 2){
          return(make_skip_res(all_leaves_in_parent, "single_child_node"))
        }

        ## 构建model matrix
        model <- tryCatch({
          if(use.glmm){
            fe_formula <- formula(paste(
              '~', ifelse(is.null(batch_col), 'condition',
                          paste('condition +', batch_col))
            ))
            m <- model.matrix(fe_formula, data = design_df)
          } else {
            m <- model.matrix(design, data = design_df)
          }
          rownames(m) <- design_df$sample
          m[colnames(cluster.count), , drop = FALSE]
        }, error = function(e){
          warning("model.matrix failed: ", e$message)
          NULL
        })

        if(is.null(model)){
          return(make_skip_res(all_leaves_in_parent, "failed:model_matrix"))
        }

        if(ncol(cluster.count) <= ncol(model)){
          return(make_skip_res(all_leaves_in_parent, "insufficient_df"))
        }

        lib_size <- colSums(cluster.count)

        ## Normalization
        if(norm.method == "TMM"){
          dge <- DGEList(counts = cluster.count, lib.size = lib_size)
          if(nrow(cluster.count) < 8){
            warning("Too few clusters (n=", nrow(cluster.count),
                    ") for stable TMM, using logMS")
            actual_norm <- "logMS"
          } else {
            dge         <- calcNormFactors(dge, method = "TMM")
            actual_norm <- "TMM"
          }
        } else if(norm.method == "RLE"){
          dge         <- DGEList(counts = cluster.count, lib.size = lib_size)
          dge         <- calcNormFactors(dge, method = "RLE")
          actual_norm <- "RLE"
        } else {
          dge         <- DGEList(counts = cluster.count, lib.size = lib_size)
          actual_norm <- "logMS"
        }

        ## NB-GLM路线
        if(!use.glmm){

          dge <- tryCatch(estimateDisp(dge, model), error = function(e) NULL)
          fit <- tryCatch(
            glmQLFit(dge, model, robust = robust, legacy = TRUE),
            error = function(e) NULL
          )

          if(is.null(fit) | is.null(dge)){
            return(make_skip_res(all_leaves_in_parent, "failed:NB_GLM"))
          }

          qlf_res <- tryCatch({
            if(!is.null(model.contrasts)){
              mod.contrast <- makeContrasts(
                contrasts = model.contrasts, levels = model
              )
              glmQLFTest(fit, contrast = mod.contrast)
            } else {
              n.coef <- ncol(model)
              glmQLFTest(fit, coef = n.coef)
            }
          }, error = function(e) NULL)

          if(is.null(qlf_res)){
            return(make_skip_res(all_leaves_in_parent, "failed:NB_GLM_test"))
          }

          test_res <- tryCatch(
            as.data.frame(topTags(
              qlf_res,
              sort.by       = 'none',
              n             = Inf,
              adjust.method = 'none'
            )),
            error = function(e) NULL
          )

          if(is.null(test_res)){
            return(make_skip_res(all_leaves_in_parent, "failed:topTags"))
          }

          ## 健壮提取PValue列名
          pval_col <- grep("^[Pp].*[Vv]alue$",
                           colnames(test_res), value = TRUE)
          if(length(pval_col) > 0 && !"PValue" %in% colnames(test_res)){
            test_res$PValue <- test_res[[pval_col[1]]]
          }

          test_res$FDR    <- NA
          test_res$method <- paste0("NB_GLM_", actual_norm)

          ## 补回被过滤的叶节点
          if(length(filtered_leaves) > 0){
            filler <- make_skip_res(filtered_leaves, "filtered:zero_counts")
            missing_cols <- setdiff(colnames(test_res), colnames(filler))
            for(mc in missing_cols) filler[[mc]] <- NA
            filler   <- filler[, colnames(test_res), drop = FALSE]
            test_res <- rbind(test_res, filler)
          }

          return(test_res)

        } else {

          ## GLMM路线
          if(is.null(batch_col)){
            warning("use.glmm=TRUE but no batch_col, falling back to NB-GLM")
            dge <- tryCatch(estimateDisp(dge, model), error = function(e) NULL)
            fit <- tryCatch(
              glmQLFit(dge, model, robust = robust, legacy = TRUE),
              error = function(e) NULL
            )
            if(is.null(fit) | is.null(dge)){
              return(make_skip_res(all_leaves_in_parent,
                                   "failed:NB_GLM_fallback"))
            }
            n.coef   <- ncol(model)
            qlf_res  <- tryCatch(
              glmQLFTest(fit, coef = n.coef),
              error = function(e) NULL
            )
            if(is.null(qlf_res)){
              return(make_skip_res(all_leaves_in_parent,
                                   "failed:NB_GLM_fallback"))
            }
            test_res <- tryCatch(
              as.data.frame(topTags(
                qlf_res,
                sort.by       = 'none',
                n             = Inf,
                adjust.method = 'none'
              )),
              error = function(e) NULL
            )
            if(is.null(test_res)){
              return(make_skip_res(all_leaves_in_parent,
                                   "failed:NB_GLM_fallback"))
            }
            test_res$FDR    <- NA
            test_res$method <- paste0("NB_GLM_fallback_", actual_norm)
            return(test_res)
          }

          z_formula <- formula(paste('~ 0 +', batch_col))
          z_model   <- tryCatch({
            zm <- model.matrix(
              z_formula,
              data = design_df[match(colnames(cluster.count),
                                     design_df$sample), ]
            )
            rownames(zm) <- colnames(cluster.count)
            zm
          }, error = function(e) NULL)

          if(is.null(z_model)){
            return(make_skip_res(all_leaves_in_parent, "failed:GLMM_zmodel"))
          }

          dge <- tryCatch(estimateDisp(dge, model), error = function(e) NULL)
          if(is.null(dge)){
            return(make_skip_res(all_leaves_in_parent, "failed:GLMM_disp"))
          }

          dispersion  <- dge$tagwise.dispersion
          offsets     <- log(dge$samples$norm.factors)
          rand.levels <- lapply(seq_len(ncol(z_model)), function(j){
            unique(z_model[, j])
          })
          names(rand.levels) <- colnames(z_model)

          glmm.ctrl <- list(
            theta.tol = max.tol,
            max.iter  = max.iters,
            solver    = glmm.solver
          )

          fit_list <- lapply(seq_len(nrow(dge$counts)), function(i){
            tryCatch(
              fitGLMM(
                X              = model,
                Z              = z_model,
                y              = dge$counts[i, ],
                offsets        = offsets,
                random.levels  = rand.levels,
                REML           = REML,
                dispersion     = 1 / dispersion[i],
                glmm.control   = glmm.ctrl,
                intercept.type = "fixed"
              ),
              error = function(e) NULL
            )
          })

          ret.beta <- ncol(model)

          glmm_res <- do.call(rbind, lapply(seq_along(fit_list), function(i){
            fit_i     <- fit_list[[i]]
            node_name <- rownames(cluster.count)[i]
            if(is.null(fit_i) || !isTRUE(fit_i$converged)){
              return(data.frame(
                logFC  = NA, logCPM = NA, F = NA,
                PValue = NA, FDR    = NA,
                method = "failed:GLMM_convergence",
                row.names = node_name
              ))
            }
            data.frame(
              logFC  = fit_i$FE[ret.beta],
              logCPM = log2(mean(dge$counts[i, ] /
                                   colSums(dge$counts)) * 1e6),
              F      = fit_i$t[ret.beta]^2,
              PValue = fit_i$PVALS[ret.beta],
              FDR    = NA,
              method = paste0("NB_GLMM_", actual_norm),
              row.names = node_name
            )
          }))

          if(length(filtered_leaves) > 0){
            filler <- make_skip_res(filtered_leaves, "filtered:zero_counts")
            missing_cols <- setdiff(colnames(glmm_res), colnames(filler))
            for(mc in missing_cols) filler[[mc]] <- NA
            filler   <- filler[, colnames(glmm_res), drop = FALSE]
            glmm_res <- rbind(glmm_res, filler)
          }

          return(glmm_res)
        }
      }
    )

    ## -------------------------------------------------------
    ## 按叶节点标签写入矩阵
    ## 行名现在是cluster_X格式，直接用字符串索引
    ## -------------------------------------------------------
    for(parent_id in names(DA.layer.list)){
      res_parent <- DA.layer.list[[parent_id]]
      if(is.null(res_parent)) next

      valid_leaves <- intersect(
        rownames(res_parent),
        rownames(logFC.matrix)   # 直接用矩阵行名匹配，更安全
      )
      if(length(valid_leaves) == 0) next

      logFC.matrix[valid_leaves, kk]  <- res_parent[valid_leaves, "logFC"]
      method.matrix[valid_leaves, kk] <- res_parent[valid_leaves, "method"]

      pval_col <- grep("^[Pp].*[Vv]alue$",
                       colnames(res_parent), value = TRUE)
      if(length(pval_col) > 0){
        PValue.matrix[valid_leaves, kk] <- res_parent[
          valid_leaves, pval_col[1]
        ]
      }
    }

    ## 逐层BH校正
    valid_pvals <- !is.na(PValue.matrix[, kk])
    if(sum(valid_pvals) > 0){
      FDR.matrix[valid_pvals, kk] <- p.adjust(
        PValue.matrix[valid_pvals, kk],
        method = "BH"
      )
    }

    message("Layer ", kk, " (", colnames(cluster.layers.test)[kk], ")",
            " done | tested: ", sum(valid_pvals),
            " nodes | n_clusters: ", n_clusters_per_layer[kk])
  }

  return(list(
    logFC           = logFC.matrix,
    FDR             = FDR.matrix,
    PValue          = PValue.matrix,
    method          = method.matrix,
    norm.method     = norm.method,
    model           = ifelse(use.glmm, "NB-GLMM", "NB-GLM"),
    n_layers_tested = n_layers_test,
    n_layers_full   = n_layers_full,
    leaf_nodes      = all_leaf_nodes,      # 整数cluster编号
    leaf_node_names = all_leaf_names       # cluster_X格式的行名
  ))
}
