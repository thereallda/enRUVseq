
#' Applies normalization on sequencing data
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is the number of samples and p is the number of features. 
#' @param method Vector of normalization methods that are applied to the data.
#'   Available methods are: \code{c("CPM", "UQ", "TMM", "DESeq", "RUV")}. 
#'   Select one or multiple methods. By default all normalization methods will be applied.
#' @param control.idx Vector of control genes' id. 
#' @param sc.idx A numeric matrix specifying the replicate samples for which to 
#' compute the count differences used to estimate the factors of unwanted variation.
#'
#' @return List of objects containing normalized data and associated normalization factors. 
#' @export
#'
ApplyNormalization <- function(data, 
                               method = c("CPM", "UQ", "TMM", "DESeq", "RLE", "RUV"),
                               control.idx = NULL,
                               sc.idx = NULL) {
  method <- match.arg(method,
                      choices = c("CPM", "UQ", "TMM", "DESeq", "RLE", "RUV"),
                      several.ok = TRUE)
  
  data.norm <- list()
  
  data.norm[["Raw"]] <- list(dataNorm = data, normFactor = rep(1, ncol(data)))
  
  if ("CPM" %in% method) {
    data.norm[["CPM"]] <- normCPM(data)
  }
  if ("UQ" %in% method) {
    data.norm[["UQ"]] <- normUQ(data)
  }
  if ("TMM" %in% method) {
    data.norm[["TMM"]] <- normTMM(data)
  }
  if ("DESeq" %in% method) {
    data.norm[["DESeq"]] <- normDESeq(data)
  }
  if ("RLE" %in% method) {
    data.norm[["RLE"]] <- normRLE(data)
  }
  if ("RUV" %in% method) {
    data.norm[["RUVg"]] <- normRUV(data, control.idx, method = "RUVg")
    data.norm[["RUVs"]] <- normRUV(data, control.idx, sc.idx, method = "RUVs")
  }
  
  return(data.norm)
}

#' Perform CPM normalization
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is 
#' the number of samples and p is the number of features. 
#'
#' @return List containing normalized counts and normalized factors for library size. 
#' @export
#'
normCPM <- function(data) {
  normFactor <- rep(1,ncol(data))
  sizeFactor <- normFactor*colSums(data)/1e6
  dataNorm <- t(t(data)/sizeFactor)
  
  return(list(
    dataNorm = dataNorm,
    normFactor = normFactor
  ))
}

#' Perform upper-quartile normalization
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is 
#' the number of samples and p is the number of features. 
#'
#' @return List containing normalized counts and normalized factors for library size. 
#' @export
#'
#' @importFrom edgeR calcNormFactors
normUQ <- function(data) {
  normFactor <- edgeR::calcNormFactors(data, method = "upperquartile")
  sizeFactor <- normFactor*colSums(data)/1e6
  dataNorm <- t(t(data)/sizeFactor)
  
  return(list(
    dataNorm = dataNorm,
    normFactor = normFactor
  ))
}

#' Perform TMM normalization
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is 
#' the number of samples and p is the number of features.
#'
#' @return List containing normalized counts and normalized factors for library size. 
#' @export
#'
#' @importFrom edgeR calcNormFactors
normTMM <- function(data) {
  normFactor <- edgeR::calcNormFactors(data, method = "TMM")
  sizeFactor <- normFactor*colSums(data)/1e6
  dataNorm <- t(t(data)/sizeFactor)
  
  return(list(
    dataNorm = dataNorm,
    normFactor = normFactor
  ))
}

#' Perform DESeq2 normalization
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is 
#' the number of samples and p is the number of features.
#'
#' @return List containing normalized counts and normalized factors for library size. 
#' @export
#'
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
normDESeq <- function(data) {
  sizeFactor <- DESeq2::estimateSizeFactorsForMatrix(data)
  normFactor <- 1e7*sizeFactor/colSums(data)
  dataNorm <- t(t(data)/sizeFactor)
  return(list(
    dataNorm = dataNorm,
    normFactor = normFactor
  ))  
}

#' Perform RLE normalization
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is 
#' the number of samples and p is the number of features.
#'
#' @return List containing normalized counts and normalized factors for library size. 
#' @export
#'
#' @importFrom edgeR calcNormFactors
normRLE <- function(data) {
  normFactor <- edgeR::calcNormFactors(data, method = "RLE")
  sizeFactor <- normFactor*colSums(data)/1e6
  dataNorm <- t(t(data)/sizeFactor)
  
  return(list(
    dataNorm = dataNorm,
    normFactor = normFactor
  ))
}

#' Perform RUV normalization 
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is 
#' the number of samples and p is the number of features.
#' @param control.idx Vector of control genes' id. 
#' @param sc.idx A numeric matrix specifying the replicate samples for which to 
#' compute the count differences used to estimate the factors of unwanted variation.
#' @param method Perform RUVg or RUVs normalization. 
#'
#' @return List containing normalized counts and adjust factors for adjusting unwanted variation. 
#' @export
#'
#' @import RUVSeq
normRUV <- function(data, 
                    control.idx = NULL, 
                    sc.idx = NULL, 
                    method = c("RUVg", "RUVs")
                    # adjustFactor = NULL,
                    # alpha = NULL
) {
  
  if (is.null(control.idx)) {
    control.idx <- rownames(data)
  }
  
  # control
  data.control <- data[control.idx, ]
  normFactor.control <- edgeR::calcNormFactors(data.control, method = "RLE")
  sizeFactor.control <- normFactor.control*colSums(data.control) / 1e6
  dataNorm.control <- t(t(data.control) / sizeFactor.control)
  
  # non-control
  if (nrow(data.control) < nrow(data)) {
    data.noncontrol <- data[rownames(data) %in% control.idx, ]
    normFactor.noncontrol <- edgeR::calcNormFactors(data.noncontrol, method = "RLE")
    sizeFactor.noncontrol <- normFactor.noncontrol*colSums(data.noncontrol) / 1e6
    dataNorm.noncontrol <- t(t(data.noncontrol) / sizeFactor.noncontrol)
    dataNorm <- rbind(dataNorm.noncontrol, dataNorm.control)
    dataNorm <- log2(dataNorm + 1)
  } 
  else {
    dataNorm <- log2(dataNorm.control + 1)
  }
  
  if (method == "RUVg") {
    # ruv.set <- RUVSeq::RUVg(dataNorm, cIdx = control.idx, k = 2, drop = 1, isLog = TRUE)
    ruv.set <- enRUVg(dataNorm, control.idx = control.idx, k = 2, drop = 1, isLog = TRUE)
    dataNorm <- 2^(ruv.set$normalizedCounts)-1
    
    # if (all(!is.null(adjustFactor), !is.null(alpha))) {
    #   dataCorrected <- t(dataNorm) - adjustFactor %*% alpha 
    #   dataNorm <- 2^(t(dataCorrected)) - 1
    #   ruv.set <- list(W = adjustFactor, alpha = alpha)
    # }
    
    return(list(
      dataNorm = dataNorm,
      adjustFactor = ruv.set$W,
      alpha = ruv.set$alpha
    ))
  }
  
  if (method == "RUVs") {
    # ruv.set <- RUVSeq::RUVs(dataNorm, cIdx = control.idx, k = 1, scIdx = sc.idx, isLog = TRUE)
    ruv.set <- enRUVs(dataNorm, control.idx = control.idx, k = 1, sc.idx = sc.idx, isLog = TRUE)
    dataNorm <- 2^(ruv.set$normalizedCounts)-1
    
    # if (all(!is.null(adjustFactor), !is.null(alpha))) {
    #   dataCorrected <- t(dataNorm) - adjustFactor %*% alpha 
    #   dataNorm <- 2^(t(dataCorrected)) - 1
    #   ruv.set <- list(W = adjustFactor, alpha = alpha)
    # }
    
    return(list(
      dataNorm = dataNorm,
      adjustFactor = ruv.set$W,
      alpha = ruv.set$alpha
    ))
  }
}
