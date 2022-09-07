#' Remove unwanted variation using control genes
#'
#' @param object A counts matrix.
#' @param log Whether to perform log2-transformation with 1 offset on data matrix, 
#' default: TRUE. 
#' @param k The number of factors of unwanted variation to be estimated from the data.
#' @param tolerance Tolerance in the selection of the number of positive singular 
#' values, i.e., a singular value must be larger than tolerance to be considered positive.
#' @param control.idx ID of control genes.
#' @param drop The number of singular values to drop in the estimation of 
#' unwanted variation, default drop the first singular value that represent the 
#' difference between enrichment and input
#'
#' @return list contain a matrix of unwanted factors (W) and corrected counts matrix (normalizedCounts).
#' @export
#'
enRUVg <- function(object, log=TRUE, k=1, tolerance=1e-8, control.idx, drop=1) {
  if(log) {
    Y <- t(log2(object+1))
  } else {
    Y <- t(object)
  }
  Ycenter <- apply(Y, 2, function(x) scale(x, center = TRUE, scale = FALSE))
  svdWa <- svd(Ycenter[, control.idx])
  first <- 1+drop
  # d	a vector containing the singular values of x, of length min(n, p), sorted decreasingly.
  k <- min(k, max(which(svdWa$d > tolerance)))
  # u	a matrix whose columns contain the left singular vectors of x, present if nu > 0. Dimension c(n, nu).
  W <- svdWa$u[, (first:k), drop = FALSE]
  colnames(W) <- paste("W", seq(1, ncol(W)), sep="_")
  
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
  correctedY <- Y - W %*% alpha
  return(list(W = W, normalizedCounts = t(correctedY), alpha = alpha))
}
