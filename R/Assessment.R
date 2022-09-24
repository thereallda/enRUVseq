#' Assessments of normalization performance
#'
#' @param data.ls List containing normalized counts and adjust factors for 
#' adjusting unwanted variation. Output of \code{ApplyNormalization}. 
#' @param bio_group Vector of index indicating the column index of samples of 
#' each biological groups in the raw/normalized count data matrix. 
#' @param assay_group Vector of index indicating the column index of 
#' enrichment and input samples in the raw/normalized count data matrix. 
#' @param batch_group Vector of index indicating the column index of 
#' each batch groups in the raw/normalized count data matrix. 
#' @param pam_krange Integer or vector of integers indicates the number of 
#' clusters for PAM clustering, default: 2:6. 
#' @param pc_k Integer indicates the metrics will be calculated in the first kth PCs, default: 3.
#' @param log Whether to perform log2-transformation with 1 offset on data matrix, 
#' default: TRUE. 
#' @param pos.eval.set Vector of genes id
#' @param neg.eval.set Vector of genes id
#'
#' @return List containing the score of metrics matrix and the ranking matrix, 
#' both sorted by the performance of methods from top to bottom. 
#' @export
#'
#' @importFrom stats dist cor lm na.omit var
#' @importFrom cluster silhouette 
#' @importFrom MatrixGenerics rowMedians colMedians colIQRs
#' @importFrom fpc pamk
AssessNormalization <- function(data.ls, bio_group=NULL, assay_group=NULL, batch_group=NULL,
                                pam_krange=2:6, pc_k=3, log=TRUE, 
                                pos.eval.set=NULL, neg.eval.set=NULL) {
  
  score.ls <- lapply(data.ls, function(x) {
    data <- as.matrix(x$dataNorm)
    # Clustering properties
    if (log) {
      data.log <- log2(data + 1)
    } else {
      data.log <- data
    }
    # PCA on expression matrix
    pca.expr <- prcomp(scale(t(data.log)))
    # compute right singular value by svd
    expr_sv <- svd(scale(t(data.log), center = TRUE, scale = TRUE),
                   nu = pc_k, nv = 0)$u
    # calculate euclidean distance in the space of first k PCs (default: 3)
    dist.pca.expr <- dist(scale(pca.expr$x[, 1:pc_k]), method = "euclidean")
    # dist.pca.expr <- dist(expr_sv, method = "euclidean")
    # silhouette width
    if (length(bio_group) == ncol(data)) {
      bio_sil <- mean(cluster::silhouette(bio_group, dist.pca.expr)[,"sil_width"])
    } else {
      bio_sil <- 0
    }
    if (length(assay_group) == ncol(data)) {
      assay_sil <- mean(cluster::silhouette(assay_group, dist.pca.expr)[,"sil_width"])
    } else {
      assay_sil <- 0
    }
    if (length(batch_group) == ncol(data)) {
      batch_sil <- mean(cluster::silhouette(batch_group, dist.pca.expr)[,"sil_width"])
    } else {
      batch_sil <- 0
    }
    
    prk <- fpc::pamk(pca.expr$x[,1:pc_k], krange=pam_krange) # PAM clustering with user specified k
    # prk <- cluster::pam(expr_sv, k=pam_k) # PAM clustering with user specified k
    pam_sil <- prk$pamobject$silinfo$avg.width
    
    # Global distribution properties
    data.log.rle <- data.log - rowMedians(data.log)
    # Mean squared Median RLE
    rle_med <- mean(colMedians(data.log.rle)^2)
    # Variance of IQR of RLE
    rle_iqr <- var(colIQRs(data.log.rle))
    
    # Association with control genes
    # wanted factors from positive set 
    if (!is.null(pos.eval.set)) {
      wv_factors <- svd(scale(t(data.log[pos.eval.set,]), center = TRUE, scale = TRUE),
                        nu = pc_k, nv = 0)$u
      # weighted coefficient of determination
      wv_cor <- 1 - sum(unlist(apply(expr_sv, 2, function(y) {
        lm(y ~ wv_factors)$residual
      })) ^ 2) / sum(scale(expr_sv, scale = FALSE) ^ 2)
    }
    
    # unwanted factors from negative set
    if (!is.null(neg.eval.set)) {
      uv_factors <- svd(scale(t(data.log[neg.eval.set,]), center = TRUE, scale = TRUE),
                        nu = pc_k, nv = 0)$u
      # weighted coefficient of determination
      uv_cor <- 1 - sum(unlist(apply(expr_sv, 2, function(y) {
        lm(y ~ uv_factors)$residual
      })) ^ 2) / sum(scale(expr_sv, scale = FALSE) ^ 2)
    }
  
    score <- c(
      BIO_SIL = bio_sil,
      ASSAY_SIL = assay_sil,
      BATCH_SIL = batch_sil,
      PAM_SIL = pam_sil,
      RLE_MED = rle_med,
      RLE_IQR = rle_iqr,
      EXP_WV_COR = wv_cor,
      EXP_UV_COR = uv_cor
      )
  })
  
  # reduce list of scores into table
  # with methods in row and measures in column
  score <- data.frame(do.call(rbind, score.ls))
  
  # multiplying by +/- 1 so that large values correspond to good performance
  performance <- t(t(score) * c(1,1,-1,1,-1,-1,1,-1))  # BIO_SIL,ASSAY_SIL,BATCH_SIL,PAM_SIL,RLE_MED,RLE_IQR,EXP_WV_COR,EXP_UV_COR
  # rank performance
  ranked_performance <- apply(na.omit(performance), 2, rank, ties.method = "min")
  # mean performance rank
  if (is.null(dim(ranked_performance))) {
    mean_performance_rank <- ranked_performance
  } else {
    # if performance all 1, remove it before scoring
    metrics.keep <- colSums(ranked_performance==1) != nrow(ranked_performance)
    mean_performance_rank <- rowMeans(ranked_performance[,metrics.keep])
  }
  
  ranked_performance <- as.data.frame(ranked_performance)
  ranked_performance$PERF_SCORE <- mean_performance_rank
  ranked_performance <- ranked_performance[order(mean_performance_rank, decreasing = TRUE),]
  
  score <- score[order(mean_performance_rank, decreasing = TRUE), ]
    
  return(list(
    score = score,
    performance = ranked_performance
  ))
}


