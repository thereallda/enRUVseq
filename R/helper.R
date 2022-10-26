#' Normalization and assessment in one function
#'
#' @param object Enone object.
#' @param auto Whether to automatically select negative control, positive 
#' evaluation and negative evaluation genes, default: TRUE. 
#' @param n.neg.control Number of negative control genes for RUV normalization, default: 1000. 
#' @param n.pos.eval Number of positive evaluation genes for wanted variation assessment, default: 500.
#' @param n.neg.eval Number of negative evaluation genes for unwanted variation assessment, default: 500.
#' @param neg.control Vector of negative control genes' id for RUV normalization, default: NULL. 
#' @param pos.eval Vector of positive evaluation genes' id for wanted variation assessment, default: NULL.
#' @param neg.eval Vector of negative evaluation genes' id for unwanted variation assessment, default: NULL.
#' @param scaling.method Vector of normalization methods that are applied to the data.
#'   Available methods are: \code{c("TC", "UQ", "TMM", "DESeq", "PossionSeq")}. 
#'   Select one or multiple methods. By default all normalization methods will be applied.
#' @param ruv.norm Whether to perform RUV normalization. 
#' @param ruv.k The number of factors of unwanted variation to be estimated from the data.
#' @param ruv.drop The number of singular values to drop in the estimation of 
#' unwanted variation, default: 0.  
#' @param pam.krange Integer or vector of integers indicates the number of 
#' clusters for PAM clustering, default: 2:6. 
#' @param pc.k Integer indicates the metrics will be calculated in the first kth PCs, default: 3.
#'
#' @return Enone object.
#' @export
#'
#' @importFrom utils head
#' @importFrom stringr str_extract
#' @importFrom stats as.formula model.matrix setNames
enONE <- function(object,
                  auto = TRUE, 
                  n.neg.control = 1000, n.pos.eval = 500, n.neg.eval = 500,
                  neg.control = NULL, pos.eval = NULL, neg.eval = NULL,
                  scaling.method = c("TC", "UQ", "TMM", "DESeq", "PossionSeq"),
                  ruv.norm = TRUE, ruv.k = 1, ruv.drop = 0,
                  pam.krange = 2:6, pc.k = 3) {
  
  # retrieve parameters from Enone object
  bio.group <- object$condition
  enrich.group <- object$enrich
  
  if (!any(is.na(object$batch))) {
    batch.group <- object$batch
  } else {
    batch.group <- NULL
  }
  
  spike.in.prefix <- object@parameter$spike.in.prefix
  input.id <- object@parameter$input.id
  enrich.id <- object@parameter$enrich.id
  synthetic.id <- object@parameter$synthetic.id
  
  # create group matrix
  sc_mat <-  CreateGroupMatrix(bio.group)
  enrich_mat <- CreateGroupMatrix(enrich.group)
  
  # get counts
  data <- SummarizedExperiment::assay(object)
  counts_nsp <- data[grep(spike.in.prefix, rownames(data), invert = TRUE),]
  counts_sp <- data[grep(spike.in.prefix, rownames(data)),]
  
  ## gene selection 
  if (auto) {
    ### 1. negative control genes for RUV
    cat(paste("The number of negative control genes for RUV:",n.neg.control,"\n"))
    designMat <- model.matrix(~0+enrich.group)
    deg.en <- edgeRDE(counts_sp,
                      group = enrich.group,
                      design.formula = as.formula("~0+condition"),
                      contrast.df = data.frame(Group1=enrich.id, Group2=input.id)
    )
    # top 1000 (default) non-sig de 
    res_tab <- deg.en$res.ls[[paste(enrich.id, input.id, sep='_')]]
    # res_tab <- subset(res_tab, FDR > 0.05)
    neg.control.set <- head(res_tab[order(res_tab$FDR, decreasing = TRUE),]$GeneID, n=n.neg.control)
    
    ### 2. positive evaluation genes (default 500)
    # if provided, preclude synthetic RNA from evaluation set 
    cat(paste("The number of positive evaluation genes:",n.pos.eval,"\n"))
    deg.en <- edgeRDE(counts_nsp[!rownames(counts_nsp) %in% synthetic.id,],
                      group = enrich.group,
                      design.formula = as.formula("~0+condition"),
                      contrast.df = data.frame(Group1=enrich.id, Group2=input.id)
    )
    res_tab <- deg.en$res.ls[[paste(enrich.id, input.id, sep='_')]]
    pos.eval.set <- head(res_tab[order(res_tab$FDR),]$GeneID, n=n.pos.eval)
    
    ### 3. negative evaluation genes (default 500)
    # if provided, preclude synthetic RNA from evaluation set 
    cat(paste("The number of negative evaluation genes:",n.neg.eval,"\n"))
    de.all <- edgeRDE(counts_nsp[!rownames(counts_nsp) %in% synthetic.id,],
                      group = bio.group,
                      design.formula = as.formula("~condition"),
                      coef = 2:length(unique(bio.group))
    )
    res_tab <- de.all$res.ls[[1]]
    neg.eval.set <- head(res_tab[order(res_tab$FDR, decreasing = TRUE),]$GeneID, n=n.neg.eval)
  } else {
    neg.control.set <- neg.control
    pos.eval.set <- pos.eval
    neg.eval.set <- neg.eval
  }
  
  # save gene set to object
  SummarizedExperiment::rowData(object)$NegControl <- SummarizedExperiment::rowData(object)$GeneID %in% neg.control.set
  SummarizedExperiment::rowData(object)$NegEvaluation <- SummarizedExperiment::rowData(object)$GeneID %in% neg.eval.set
  SummarizedExperiment::rowData(object)$PosEvaluation <- SummarizedExperiment::rowData(object)$GeneID %in% pos.eval.set
  
  ## apply normalization 
  cat("Apply normalization...\n")
  norm.ls <- ApplyNormalization(data,
                                scaling.method = scaling.method, 
                                ruv.norm = ruv.norm, ruv.k = ruv.k, ruv.drop = ruv.drop,
                                spike.in.prefix = spike.in.prefix,
                                # below parameters are generated inside function
                                control.idx = neg.control.set, 
                                sc.idx = sc_mat, 
                                enrich.idx = enrich_mat)
  
  ## assessment 
  bio_group_index <- as.numeric(factor(bio.group, levels=unique(bio.group)))
  assay_group_index <- as.numeric(factor(enrich.group, levels=unique(enrich.group)))
  if (!is.null(batch.group)) {
    batch_group_index <- as.numeric(factor(batch.group, levels=unique(batch.group)))
  } else {
    batch_group_index <- NULL
  }
  cat("Perform assessment...\n")
  norm.eval <- AssessNormalization(norm.ls,
                                   pam.krange = pam.krange,
                                   pc.k = pc.k,
                                   batch.group = batch_group_index,
                                   # below parameters are created inside function
                                   bio.group = bio_group_index, 
                                   assay.group = assay_group_index, 
                                   pos.eval.set = pos.eval.set,
                                   neg.eval.set = neg.eval.set)
  
  ## save metrics to Enone object
  object@enone_metrics <- norm.eval$metrics
  ## save score to Enone object
  object@enone_score <- norm.eval$score
  ## add run parameter in Enone object
  parameter.run <- list(
    n.neg.control=n.neg.control,
    n.pos.eval=n.pos.eval,
    n.neg.eval=n.neg.eval,
    scaling.method=scaling.method,
    ruv.norm=ruv.norm,
    ruv.k=ruv.k,
    ruv.drop=ruv.drop,
    pam.krange=pam.krange,
    pc.k=pc.k
  )
  object@parameter <- c(object@parameter, parameter.run)

  # only store normalization method names in object for reducing memory cost
  norm.methods <- names(norm.ls)
  object@counts$sample <- stats::setNames(vector('list', length(norm.methods)), nm=norm.methods)
  object@enone_factor$sample <- stats::setNames(vector('list', length(norm.methods)), nm=norm.methods)
  # except 'Raw' matrix
  Counts(object, slot='sample', method='Raw') <- counts_nsp
  Counts(object, slot='spike_in', method='Raw') <- counts_sp
  
  return(object)
}

#' Find enriched genes between enrich and input samples
#'
#' @param object An object.
#' @param slot Which slot, one of \code{sample} or \code{spike_in}.  
#' @param method Which normalization methods, must be one of the methods presented in the selected slot.  
#' @param logfc.cutoff Filter genes by at least X-fold difference (log2-scale) 
#' between the two groups of samples, default: 1. 
#' @param p.cutoff Filter genes by no more than Y adjusted p-value, default: 0.05. 
#' @param ... Additional parameters can be passed to \code{edgeRDE}. 
#' 
#' @name FindEnrichment 
#' @aliases FindEnrichment  FindEnrichment,Enone,character,character,numeric,numeric-method
#' 
#' @return data.frame
#' @export
#'
FindEnrichment <-  function(object, slot=c("sample","spike_in"), method, 
                            logfc.cutoff=1, p.cutoff=0.05, ...) {
  
  contrast_df <- data.frame(Group1 = unique(grep(object@parameter$enrich.id, object$condition, value = TRUE)),
                            Group2 = unique(grep(object@parameter$input.id, object$condition, value = TRUE)))
  # extract sample or spike-in counts
  slot <- match.arg(slot, choices = c("sample","spike_in"))
  
  # test if chosen method in object 
  if (is.null(names(object@enone_factor[[slot]]))) {
    stop("Normalizations for ", slot, " not found. At least one normalization should be performed.")
  }
  method <- match.arg(method, choices = names(object@enone_factor[[slot]]))
  
  if (slot == "spike_in") {
    counts_df <- SummarizedExperiment::assay(object)[SummarizedExperiment::rowData(object)$SpikeIn,]
  } 
  else if (slot == "sample" & any(SummarizedExperiment::rowData(object)$Synthetic)) {
    counts_df <- SummarizedExperiment::assay(object)[!SummarizedExperiment::rowData(object)$SpikeIn & !SummarizedExperiment::rowData(object)$Synthetic,]
  } 
  else {
    counts_df <- SummarizedExperiment::assay(object)[!SummarizedExperiment::rowData(object)$SpikeIn,]
  }
  
  # get list of factors 
  factor.ls <- getFactor(object, slot=slot, method=method)
  if ("normFactor" %in% names(factor.ls)) {
    if ("adjustFactor" %in% names(factor.ls)) {
      # if norm factors and adjust factors were both provided
      de <- edgeRDE(counts = counts_df,
                    group = object$condition,
                    contrast.df = contrast_df,
                    norm.factors = factor.ls$normFactor,
                    adjust.factors = factor.ls$adjustFactor,
                    logfc.cutoff=logfc.cutoff, p.cutoff=p.cutoff, ...)
    } else {
      # if only norm factors were provided
      de <- edgeRDE(counts = counts_df,
                    group = object$condition,
                    contrast.df = contrast_df,
                    norm.factors = factor.ls$normFactor,
                    design.formula = as.formula("~0+condition"),
                    logfc.cutoff=logfc.cutoff, p.cutoff=p.cutoff, ...)
    }
  } else {
    stop("One or both of 'normFactor' and 'adjustFacotr' should be provided.")
  }
  # save enrichment in Enone object
  object@enrichment[[slot]] <- de$res.ls
  object@enrichment_filtered[[slot]] <- de$res.sig.ls
  
  return(object)
}

#' Apply specific normalization method
#'
#' @param object Enone object. 
#' @param slot Which slot, one of \code{sample} or \code{spike_in}.  
#' @param method Which normalization methods to perform. 
#'
#' @return object
#' @export
#' 
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom stringr str_extract
UseNormalization <- function(object, slot=c("sample","spike_in"), method) {
  
  # sample or spike-in 
  slot <- match.arg(slot, choices = c("sample","spike_in"))
  
  # get raw counts
  if (slot == "spike_in") {
    # only spike-in counts
    counts_df <- SummarizedExperiment::assay(object)[SummarizedExperiment::rowData(object)$SpikeIn,]
  } 
  else {
    # All counts
    counts_df <- SummarizedExperiment::assay(object)
  }
  
  # method.curr[1]: scaling; method.curr[2]: RUV; method.curr[3]: number of k
  method.curr <- unlist(strsplit(method, split = "_"))
  
  # check whether selected method can be provided
  if (!method.curr[1] %in% c('Raw','TC','UQ','DESeq','TMM','PossionSeq')) {
    stop('Scaling method: ', method.curr[1], ' is not provided. It should be one of ', c('Raw','TC','UQ','DESeq','TMM','PossionSeq'))
  }
  if (length(method.curr) > 1 & !method.curr[2] %in% c('RUVg','RUVs','RUVse')) {
    stop('RUV method: ', method.curr[2], ' is not provided. It should be one of ', c('RUVg','RUVs','RUVse'))
  }
  if (length(method.curr) > 2 & is.na(stringr::str_extract(method.curr[3],'k'))) {
    stop('Use x number of factors to estimate unwanted variation as "kx", but not ', method.curr[3])
  }  
  if (length(method.curr) > 2 & as.numeric(stringr::str_extract(method.curr[3],'\\d')) > ncol(object)) {
    stop('Number of required factors exceed, try least.')
  } 
  
  # scaling
  neg.control <- getGeneSet(object, "NegControl")
  if (method.curr[1] == "Raw") {
    counts_scale <- list(dataNorm=counts_df, normFactor=rep(1, ncol(counts_df)))
  } 
  else {
    normScaling <- get(paste0("norm", method.curr[1]))
    
    if (slot == "spike_in") {
      counts_scale <- normScaling(counts_df)
    } 
    else {
      counts_sp_scale <- normScaling(counts_df[SummarizedExperiment::rowData(object)$SpikeIn,])
      counts_nsp_scale <- normScaling(counts_df[!SummarizedExperiment::rowData(object)$SpikeIn,])
      dataNorm <- rbind(counts_nsp_scale$dataNorm, counts_sp_scale$dataNorm[neg.control,])
      # counts_scale contain both normalized data and normalization factors
      counts_scale <- list(dataNorm = dataNorm, 
                           normFactor = counts_nsp_scale$normFactor)
    }
  }
  
  # RUV
  sc_idx <- CreateGroupMatrix(object$condition)
  enrich_idx <- CreateGroupMatrix(object$enrich)
  if (method.curr[2] %in% c("RUVg", "RUVs", "RUVse")) {
    
    sc.idx <- switch(method.curr[2],
                     "RUVg" = NULL,
                     "RUVs" = sc_idx,
                     "RUVse" = enrich_idx)
    
    if (slot == "spike_in") {
      counts_norm <- normRUV(counts_scale$dataNorm,
                             control.idx = neg.control,
                             sc.idx = sc.idx,
                             method = method.curr[2],
                             k = as.numeric(gsub("k", "", method.curr[3])))
    } 
    else {
      counts_norm <- normRUV(counts_scale$dataNorm,
                             control.idx = neg.control,
                             sc.idx = sc.idx,
                             method = method.curr[2],
                             k = as.numeric(gsub("k", "", method.curr[3])))
      # return only non-spike-in counts
      counts_norm$dataNorm <- counts_norm$dataNorm[!rownames(counts_norm$dataNorm) %in% neg.control,]
    }
    
  } 
  else {
    counts_norm <- counts_scale
  }
  
  Counts(object, slot=slot, method=method) <- counts_norm$dataNorm
  
  object@enone_factor[[slot]][[method]] <- list(normFactor=counts_scale$normFactor,
                                                adjustFactor=counts_norm$adjustFactor)
  return(object)
}

#' Enrichment level of synthetic RNA 
#'
#' @param object Enone object. 
#' @param method Which normalization methods to perform. 
#'
#' @return enrichment level of synthetic RNA
#' @export
#'
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom stringr str_extract
synEnrichment <- function(object, method="TC") {
  
  if (!any(SummarizedExperiment::rowData(object)$Synthetic)) {
    stop("Synthetic RNA id not provided in the object")
  }
  
  # get raw counts
  counts_df <- SummarizedExperiment::assay(object)
  
  # method.curr[1]: scaling; method.curr[2]: RUV; method.curr[3]: number of k
  method.curr <- unlist(strsplit(method, split = "_"))
  
  # check whether selected method can be provided
  if (!method.curr[1] %in% c('Raw','TC','UQ','DESeq','TMM')) {
    stop('Scaling method: ', method.curr[1], ' is not provided. It should be one of ', c('Raw','TC','UQ','DESeq','TMM'))
  }
  if (length(method.curr) > 1 & !method.curr[2] %in% c('RUVg','RUVs','RUVse')) {
    stop('RUV method: ', method.curr[2], ' is not provided. It should be one of ', c('RUVg','RUVs','RUVse'))
  }
  if (length(method.curr) > 2 & is.na(stringr::str_extract(method.curr[3],'k'))) {
    stop('Use x number of factors to estimate unwanted variation as "kx", but not ', method.curr[3])
  }  
  if (length(method.curr) > 2 & as.numeric(stringr::str_extract(method.curr[3],'\\d')) > ncol(object)) {
    stop('Number of required factors exceed, try least.')
  } 
  
  # scaling
  if (method.curr[1] == "Raw") {
    counts_scale <- list(dataNorm=counts_df, normFactor=rep(1,ncol(counts_df)))
  } else {
    normScaling <- get(paste0("norm", method.curr[1]))
    counts_scale <- normScaling(counts_df)
  }
  
  # RUV
  sc_idx <- CreateGroupMatrix(object$condition)
  enrich_idx <- CreateGroupMatrix(object$enrich)
  if (method.curr[2] %in% c("RUVg", "RUVs", "RUVse")) {
    
    sc.idx <- switch(method.curr[2],
                     "RUVg" = NULL,
                     "RUVs" = sc_idx,
                     "RUVse" = enrich_idx)
    counts_norm <- normRUV(counts_scale$dataNorm,
                           control.idx = getGeneSet(object, "NegControl"),
                           sc.idx = sc.idx,
                           method = method.curr[2],
                           k = as.numeric(gsub("k", "", method.curr[3])))
    
  } else {
    counts_norm <- counts_scale
  }
  
  # calculate enrichment of synthetic RNA
  counts_norm <- counts_norm$dataNorm
  syn_id <- SummarizedExperiment::rowData(object)$Synthetic
  syn_en <- counts_norm[syn_id,enrich_idx[1,]]/counts_norm[syn_id,enrich_idx[2,]]
  
  return(syn_en)
}


#' Create a matrix for RUVSeq
#'
#' @param group.vec A vector indicating membership in a group.
#'
#' @return A matrix. 
#' @export
#'
#' @examples 
#' CreateGroupMatrix(c('a','b','b','c','c','c','a','d','d'))
CreateGroupMatrix <- function(group.vec) {
  group.vec <- factor(group.vec)
  group.mat <- matrix(-1, nrow = length(levels(group.vec)), ncol = max(table(group.vec)))
  for (i in 1:length(levels(group.vec))) {
    idxs <- which(group.vec == levels(group.vec)[i])
    group.mat[i, 1:length(idxs)] <- idxs
  }
  group.mat
}

#' Count numbers of each members
#'
#' @param group.vec Vector of members. 
#'
#' @return vector
#' @export
#'
#' @examples 
#' countReplicate(c('a','b','b','c','c','c','a','d','d'))
countReplicate <- function(group.vec) {
  group.vec <- factor(group.vec, levels = unique(group.vec)) # keep factor levels input order
  rep.vec <- vector("double")
  for (i in levels(group.vec)) {
    curr.group <- grep(i, group.vec, value = TRUE)
    curr.reps <- which(curr.group == i)
    rep.vec <- c(rep.vec,curr.reps)
  }
  rep.vec
}

#' Calculate size factors for scaling raw library size
#'
#' @param data A counts matrix.
#' @param method Normalization methods to be used. Should be one of c("TMM","TMMwsp","RLE","upperquartile","none","spikein").
#' @param control.prefix Prefix of control gene id.
#'
#' @return values of size factors
#' @export 
#'
#' @importFrom edgeR calcNormFactors
calcSizeFactors <- function(data, method=NULL, control.prefix=NULL) {
  if (!is.null(control.prefix)) {
    counts <- data[grep(control.prefix, rownames(data)),]
    counts_nsp <- data[grep(control.prefix, rownames(data), invert=T),]
    counts_all <- data
  } else {
    counts <- data
  }
  
  if (method %in% c("TMM","TMMwsp","RLE","upperquartile","none")) {
    size_factor <- edgeR::calcNormFactors(counts, method=method)
  } else if (method == 'spikein'){
    spike_in_factor <- colSums(counts)/colSums(counts_all)
    size_factor <- spike_in_factor/colSums(counts_nsp)
    size_factor <- size_factor/prod(size_factor)^(1/length(size_factor))
  }
  return(size_factor)
}

#' Calculate pairwise set similarity using jaccard index
#'
#' @param x vector of set 1.
#' @param y vector of set 2.
#'
#' @return pairwise similarity
#'
JacIdx <- function(x, y) {
  int <- length(intersect(x, y))
  uni <- length(union(x, y))
  int/uni
}

#' Calculate pairwise set similarity using jaccard index 
#'
#' @param set.ls list of sets
#' @param pair If TRUE return pairwise similarity, otherwise return overall similarity (average). 
#'
#' @return Vector of similarity
#' @export
#'
setSimilarity <- function(set.ls, pair=FALSE) {
  sim.mat <- matrix(1, nrow = length(set.ls), ncol = length(set.ls))
  for (i in 1:(length(set.ls) - 1)) {
    for (j in (i+1):length(set.ls)) {
      jidx <- JacIdx(set.ls[[i]]$GeneID, set.ls[[j]]$GeneID)
      sim.mat[i, j] <- jidx
      sim.mat[j, i] <- jidx
    }
  }
  if (pair) {
    sim.mat
  } else {
    sum(sim.mat[upper.tri(sim.mat)]) / sum(upper.tri(sim.mat))  
  }
}


#' Visualization of average jaccard index in different k
#'
#' @param k.cor.vec vector of overall similarity in different k.
#' @param ref.cor value of reference similarity.
#'
#' @return ggplot2 object
#' @export
#' 
#' @import dplyr
#' @import ggplot2
JCPlot <- function(k.cor.vec, ref.cor=NULL) {
  data.frame(
    id=seq_along(k.cor.vec),
    cor=k.cor.vec) %>% 
    ggplot(aes(id, cor)) +
    geom_point() +
    geom_hline(yintercept = ref.cor, color='red') +
    theme_classic() +
    theme(axis.text = element_text(color='black')) +
    scale_x_continuous(breaks=seq_along(k.cor.vec), labels=seq_along(k.cor.vec)) +
    labs(x='k', y="Average Jaccard Similarity")
}

#' PCA plot from counts matrix
#'
#' @param object A count matrix.
#' @param use.pc Which two PCs to be used, default PC1 in x-axis and PC2 in y-axis.
#' @param color Vector indicates the color mapping of samples, default NULL.
#' @param label Vector of sample names or labels, default NULL.
#' @param shape Vector indicates the shape mapping of samples, default NULL.
#' @param vst.norm Whether to perform \code{vst} transformation, default FALSE.
#' @param palette The color palette for different groups.
#' @param repel Whether to use \code{ggrepel} to avoid overlapping text labels or not, default TRUE.
#'
#' @return ggplot2 object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom DESeq2 vst
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats prcomp
#' @importFrom paintingr paint_palette
ggPCA <- function(object, use.pc=c(1,2),
                  color=NULL, label=NULL, shape=NULL,
                  vst.norm=FALSE, palette=NULL, repel=TRUE) {
  if (vst.norm) {
    counts_norm <- DESeq2::vst(as.matrix(object))
  } else {
    counts_norm <- object
  }
  
  # perform PCA
  pca <- prcomp(t(counts_norm))
  pc.var <- round(summary(pca)$importance[2,], 3)
  pca_dat <- as.data.frame(pca$x)
  
  # check if use.pc exceed the range of pcs
  use.pc <- paste0('PC', use.pc)
  if (!all(use.pc %in% colnames(pca_dat))) {
    stop(use.pc, "exceed the range of PCs.")
  }
  # mapping data
  var.ls <- list(color = color,
                 shape = shape
  )
  var.length <- unlist(lapply(var.ls, length))
  var.ls <- var.ls[var.length == max(var.length)]
  map_df <- as.data.frame(Reduce(cbind, var.ls))
  colnames(map_df) <- names(var.ls)
  # combine with pca_dat if not empty
  if (!any(dim(map_df) == 0)) {
    pca_dat <- cbind(pca_dat, map_df)
  }
  
  # generate color palette
  if (is.null(palette)) {
    palette <- paintingr::paint_palette("Spring", length(unique(pca_dat$color)), 'continuous')
  }
  
  # create aes mapping
  map_ls <- list(x = use.pc[1],
                 y = use.pc[2],
                 color = "color",
                 shape = "shape")
  mapping <- do.call(ggplot2::aes_string, map_ls)
  
  p <- ggplot(pca_dat, mapping) +
    geom_point(size=3) +
    geom_vline(xintercept=0, color='grey80', lty=2) +
    geom_hline(yintercept=0, color='grey80', lty=2) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = 'top',
          axis.text = element_text(color='black')) +
    scale_color_manual(values = palette) +
    labs(x=paste0(use.pc[1], ': ', pc.var[1]*100, '%'),
         y=paste0(use.pc[2], ': ', pc.var[2]*100, '%'))
  
  # add text label
  if (!is.null(label)) {
    if (repel) {
      p <- p + ggrepel::geom_text_repel(label=label, max.overlaps = 20, color='black')
    } else {
      p <- p + geom_text(label=label, color='black')
    }
  }
  return(p)
}


#' Biplot of individuals and variables
#'
#' @param object PCA object returned by \code{prcomp}. 
#' @param score Vector of performance scores from \code{Enone} evaluation score. 
#' @param pt.label Whether to plot the point labels, default: TRUE. 
#' @param interactive Whether to demonstrate the plot interactively, default: FALSE. 
#'
#' @return plot
#' @export
#'
#' @import ggplot2
#' @importFrom paintingr paint_palette
#' @importFrom ggrepel geom_text_repel
#' @importFrom plotly ggplotly
ggPCA_Biplot <- function(object, score, pt.label=TRUE, interactive=FALSE) {
  
  # get data matrix 
  X <- object$x %*% solve(object$rotation)
  
  pc.var <- round(summary(object)$importance[2,], 3)
  pc.score <- as.data.frame(object$x)
  pc.score$method.id <- rownames(pc.score)
  pc.score$Performance <- score
  
  p <- ggplot(pc.score, aes(PC1, PC2, color=Performance, text=method.id)) +
    geom_point(size=3) +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.text = element_text(color='black')) +
    geom_hline(yintercept=0, lty='dashed') +
    geom_vline(xintercept=0, lty='dashed') +
    scale_color_gradientn(colors = paint_palette('Vesuvius', 100, 'continuous')) +
    labs(x=paste0('Dim1: ', pc.var[1]*100, '%'),
         y=paste0('Dim2: ', pc.var[2]*100, '%'))
  
  # add arrow, inspired by: 
  # https://stats.stackexchange.com/questions/276645/arrows-of-underlying-variables-in-pca-biplot-in-r
  for (i in 1:ncol(X)) {
    x.cord <- cor(X[,i], object$x[,1]) * sqrt(nrow(X)-1) * 0.4
    y.cord <- cor(X[,i], object$x[,2]) * sqrt(nrow(X)-1) * 0.4
    p <- p +   
      annotate("segment", x=0, y=0,
               xend=x.cord*0.9, 
               yend=y.cord*0.9, 
               arrow=arrow(), color="#4F99B4") +
      annotate("text", x=x.cord, y=y.cord, 
               label=colnames(X)[i], color="#4F99B4") 
  }
  
  if (pt.label) {
    p <- p + ggrepel::geom_text_repel(aes(label=method.id), max.overlaps = 20, color='black')
  }
  
  if (interactive) {
    p <- plotly::ggplotly(p)
  }
  
  return(p)
}

#' Combine list of DE results
#'
#' @param res.ls Named list of differential analysis results tables. 
#' Each elements in the list correspond to a table of differential analysis 
#' results between two groups of samples. 
#' @param fc.col Name of the fold change column. The fold change column name 
#' must be consistent with your results tables. 
#' @param levels Factor levels of the groups, default order by the element order of \code{res.ls}. 
#'
#' @return data.frame
#' @export
#'
#' @import dplyr
reduceRes <- function(res.ls, fc.col, levels=names(res.ls)) {
  df <- data.frame()
  for (id in names(res.ls)) {
    curr <- res.ls[[grep(id, names(res.ls), value=TRUE)]] 
    # curr$GeneID <- rownames(curr)
    df1 <- curr %>% 
      dplyr::mutate(Group = factor(rep(id, nrow(curr)), levels = levels)) %>% 
      dplyr::select(GeneID, !!sym(fc.col), Group)
    df <- rbind(df, df1)
  }
  return(df)
}

#' Box-violin plot comparing values between groups
#'
#' @param data A data frame (or a tibble).
#' @param x The grouping variable from the \code{data}.
#' @param y The value variable from the \code{data}.
#' @param color The color variable from the \code{data}.
#' @param palette The color palette for different groups.
#' @param test Perform "wilcox.test" or "t.test" or no test.
#' @param step.increase numeric vector with the increase in fraction of total height for every additional comparison to minimize overlap.
#' @param comparisons	A list of length-2 vectors specifying the groups of interest to be compared. For example to compare groups "A" vs "B" and "B" vs "C", the argument is as follow: comparisons = list(c("A", "B"), c("B", "C"))
#' @return ggplot2 object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import rstatix
#' @importFrom paintingr paint_palette
#' @importFrom ggpubr stat_pvalue_manual
#' @importFrom stats as.formula
#' @importFrom rlang .data
BetweenStatPlot <- function(data, x, y, color, palette = NULL,
                            test = c('wilcox.test', 't.test', 'none'),
                            comparisons = NULL,
                            step.increase=0.3) {
  stat.formula <- as.formula(paste(y, "~", x))
  
  test <- match.arg(test, choices = c('wilcox.test', 't.test', 'none'))
  if (test != 'none') {
    if (test == 'wilcox.test') {
      stat_dat <- data %>%
        wilcox_test(stat.formula, comparisons = comparisons)
    }
    if (test == 't.test') {
      stat_dat <- data %>%
        t_test(stat.formula, comparisons = comparisons)
    }
    stat_dat <- stat_dat %>%
      adjust_pvalue() %>%
      p_format(.data$p.adj, digits = 2, leading.zero = FALSE,
               trailing.zero = TRUE, add.p = TRUE, accuracy = 2e-16) %>%
      add_xy_position(x = x, dodge=0.8, step.increase=step.increase)
  }
  
  x.labs <- paste0(unique(data[,x]), "\n(n=", tabulate(as.factor(data[,x])),")")
  x.num <- length(unique(data[,color])) # number of x types
  if (is.null(palette)) palette <- paint_palette("Spring", x.num, 'continuous')
  
  p <- data %>%
    ggplot(aes_string(x, y, color = color)) +
    geom_violin(width = 0.8) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.text = element_text(color='black')) +
    scale_color_manual(values = palette) +
    scale_x_discrete(labels = x.labs) +
    labs(x='')
  
  if (exists('stat_dat')) {
    p <- p + stat_pvalue_manual(data = stat_dat, label = "p.adj", tip.length = 0.01, size = 3)
  }
  
  return(p)
}


#' Assess differential analysis performance
#'
#' @param de.raw.res Raw counts DE results
#' @param de.norm.res Normalized counts DE results
#'
#' @return List of similarity
#' @export
#'
AssessDiffAnalysis <- function(de.raw.res, de.norm.res) {
  
  # raw.sim <- setSimilarity(de.raw.res, pair = TRUE)
  # norm.sim <- setSimilarity(de.norm.res, pair = TRUE)
  # sim.diff <- calcCorDiff(raw.cor = raw.sim, norm.cor = norm.sim)
  raw.sim <- mean(setSimilarity(de.raw.res, pair = TRUE)[-1,1])
  norm.sim <- mean(setSimilarity(de.norm.res, pair = TRUE)[-1,1])
  sim.diff <- (norm.sim-raw.sim)/raw.sim  
  
  return(sim.diff)
}

#' Calculate difference of correlation
#'
#' @param raw.cor Correlation matrix (p x p) computed from raw count matrix, where p is the number of samples. 
#' @param norm.cor Correlation matrix (p x p) computed from normalized count matrix, where p is the number of samples. 
#'
#' @return Vector of correlation difference
#' @export
#'
calcCorDiff <- function(raw.cor, norm.cor) {
  mean.raw.cor <- sum(raw.cor[upper.tri(raw.cor)])/sum(upper.tri(raw.cor))
  mean.norm.cor <- sum(norm.cor[upper.tri(norm.cor)])/sum(upper.tri(norm.cor))
  return((mean.norm.cor - mean.raw.cor) / mean.raw.cor)
}

#' Applying Differential analysis
#'
#' @param data.raw A un-normalized count data matrix of shape n x p, where n is 
#' the number of samples and p is the number of features.
#' @param data.norm.ls List containing normalized counts and adjust factors for 
#' adjusting unwanted variation. 
#' @param group Vector of length p mapping the columns of \code{data.raw} to 
#' corresponding samples group. 
#'
#' @return List containing differential analysis object, result table and filtered result table.  
#' @export
#' @importFrom stringr str_extract
DiffAnalysis <- function(data.raw, 
                         data.norm.ls,
                         group) {
  de.ls <- list()
  
  contrast_df <- data.frame(Group1 = unique(grep("Enrich", group, value = TRUE)),
                            Group2 = unique(grep("Input", group, value = TRUE)))
  # compare only enrich and input
  contrast_ei <- data.frame(Group1 = "Enrich", Group2 = "Input")
  group_ei <- stringr::str_extract(group, '(Input)|(Enrich)')
  
  normalization.methods <- names(data.norm.ls)
  
  if ("CPM" %in% normalization.methods) {
    design.formula <- as.formula("~0+condition")
    de.ls[["CPM"]] <- edgeRDE(counts = data.raw,
                              group = group,
                              design.formula = design.formula,
                              contrast.df = contrast_df,
                              norm.factors = data.norm.ls[["CPM"]]$normFactor)
    
    de.tmp <- edgeRDE(counts = data.raw,
                      group = group_ei,
                      design.formula = design.formula,
                      contrast.df = contrast_ei,
                      norm.factors = data.norm.ls[["CPM"]]$normFactor)
    de.ls[["CPM"]]$res.sig.ls <- c(de.tmp$res.sig.ls, de.ls[["CPM"]]$res.sig.ls)
    
  }
  
  if ("UQ" %in% normalization.methods) {
    design.formula <- as.formula("~0+condition")
    de.ls[["UQ"]] <- edgeRDE(counts = data.raw,
                              group = group,
                              design.formula = design.formula,
                              contrast.df = contrast_df,
                              norm.factors = data.norm.ls[["UQ"]]$normFactor)
    
    de.tmp <- edgeRDE(counts = data.raw,
                      group = group_ei,
                      design.formula = design.formula,
                      contrast.df = contrast_ei,
                      norm.factors = data.norm.ls[["UQ"]]$normFactor)
    de.ls[["UQ"]]$res.sig.ls <- c(de.tmp$res.sig.ls, de.ls[["UQ"]]$res.sig.ls)
  }
  
  if ("TMM" %in% normalization.methods) {
    design.formula <- as.formula("~0+condition")
    de.ls[["TMM"]] <- edgeRDE(counts = data.raw,
                             group = group,
                             design.formula = design.formula,
                             contrast.df = contrast_df,
                             norm.factors = data.norm.ls[["TMM"]]$normFactor)
    
    de.tmp <- edgeRDE(counts = data.raw,
                      group = group_ei,
                      design.formula = design.formula,
                      contrast.df = contrast_ei,
                      norm.factors = data.norm.ls[["TMM"]]$normFactor)
    de.ls[["TMM"]]$res.sig.ls <- c(de.tmp$res.sig.ls, de.ls[["TMM"]]$res.sig.ls)
  }
  
  # if ("DESeq" %in% normalization.methods) {
  #   design.formula <- as.formula("~condition")
  #   de.ls[["DESeq"]] <- DESeq2DE(counts = data.raw, 
  #                                group = group,
  #                                design.formula = design.formula, 
  #                                contrast.df = contrast_df)
  #   
  #   de.tmp <- DESeq2DE(counts = data.raw,
  #                     group = group_ei,
  #                     design.formula = design.formula,
  #                     contrast.df = contrast_ei)
  #   de.ls[["DESeq"]]$res.sig.ls <- c(de.tmp$res.sig.ls, de.ls[["DESeq"]]$res.sig.ls)
  # }
  
  if ("DESeq" %in% normalization.methods) {
    design.formula <- as.formula("~0+condition")
    de.ls[["DESeq"]] <- edgeRDE(counts = data.raw,
                              group = group,
                              design.formula = design.formula,
                              contrast.df = contrast_df,
                              norm.factors = data.norm.ls[["DESeq"]]$normFactor)
    
    de.tmp <- edgeRDE(counts = data.raw,
                      group = group_ei,
                      design.formula = design.formula,
                      contrast.df = contrast_ei,
                      norm.factors = data.norm.ls[["DESeq"]]$normFactor)
    de.ls[["DESeq"]]$res.sig.ls <- c(de.tmp$res.sig.ls, de.ls[["DESeq"]]$res.sig.ls)
  }
  
  if ("RUVg" %in% normalization.methods) {
    de.ls[["RUVg"]] <- edgeRDE(counts = data.raw,
                             group = group,
                             contrast.df = contrast_df,
                             adjust.factors = data.norm.ls[["RUVg"]]$adjustFactor)
    
    de.tmp <- edgeRDE(counts = data.raw,
                      group = group_ei,
                      contrast.df = contrast_ei,
                      adjust.factors = data.norm.ls[["RUVg"]]$adjustFactor)
    de.ls[["RUVg"]]$res.sig.ls <- c(de.tmp$res.sig.ls, de.ls[["RUVg"]]$res.sig.ls)
  }
  
  if ("RUVs" %in% normalization.methods) {
    de.ls[["RUVs"]] <- edgeRDE(counts = data.raw,
                               group = group,
                               contrast.df = contrast_df,
                               adjust.factors = data.norm.ls[["RUVs"]]$adjustFactor)
    
    de.tmp <- edgeRDE(counts = data.raw,
                      group = group_ei,
                      contrast.df = contrast_ei,
                      adjust.factors = data.norm.ls[["RUVs"]]$adjustFactor)
    de.ls[["RUVs"]]$res.sig.ls <- c(de.tmp$res.sig.ls, de.ls[["RUVs"]]$res.sig.ls)
  }
  if ("RUVse" %in% normalization.methods) {
    de.ls[["RUVse"]] <- edgeRDE(counts = data.raw,
                               group = group,
                               contrast.df = contrast_df,
                               adjust.factors = data.norm.ls[["RUVse"]]$adjustFactor)
    
    de.tmp <- edgeRDE(counts = data.raw,
                      group = group_ei,
                      contrast.df = contrast_ei,
                      adjust.factors = data.norm.ls[["RUVse"]]$adjustFactor)
    de.ls[["RUVse"]]$res.sig.ls <- c(de.tmp$res.sig.ls, de.ls[["RUVse"]]$res.sig.ls)
  }
  return(de.ls)
}

#' Wrapper of DESeq2 procedure
#'
#' @param counts A un-normalized counts data matrix. 
#' @param group Vector of length p mapping the columns of \code{counts} to 
#' corresponding samples group. 
#' @param design.formula Formula
#' @param contrast.df Data frame of contrast, where extracting results as 
#' first column vs. second column. 
#'
#' @return List containing differential analysis object, result table and filtered result table.  
#' @export
#'
#' @import DESeq2
#' @import dplyr
DESeq2DE <- function(counts, 
                     group, 
                     design.formula, 
                     contrast.df) {
  
  col.data <- data.frame(condition = group, row.names = colnames(counts))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                colData = col.data,
                                design = design.formula)
  de <- DESeq2::DESeq(dds)
  # extract DE results
  res.ls <- list()
  for (i in 1:nrow(contrast.df)) {
    namei <- paste(contrast.df[i,], collapse = '-')
    res.ls[[namei]] <- DESeq2::results(de, 
                               contrast = c('condition', 
                                            contrast.df[i,1], 
                                            contrast.df[i,2]), 
                               tidy=TRUE) %>% 
      dplyr::rename(GeneID = row)
  }
  # cutoff for significant DEGs
  res.sig.ls <- lapply(res.ls, function(x) { x[x$log2FoldChange >=1 & x$padj < 0.05,] })
  
  return(list(de.obj = de, res.ls = res.ls, res.sig.ls = res.sig.ls))
}

#' Wrapper of edgeR procedure
#'
#' @param counts A un-normalized counts data matrix. 
#' @param group Vector of length p mapping the columns of \code{counts} to 
#' corresponding samples group. 
#' @param norm.factors Vector of normalization factors with p length. 
#' @param adjust.factors Matrix with each column indicates the adjusting factors 
#' that estimated from RUV. 
#' @param design.formula Formula
#' @param contrast.df Data frame of contrast, where extracting results as 
#' first column vs. second column. 
#' @param coef Integer or character vector indicating which coefficients of the 
#' linear model are to be tested equal to zero. Values must be columns or column names of design. 
#' @param logfc.cutoff Filter genes by at least X-fold difference (log2-scale) 
#' between the two groups of samples, default: 1. 
#' @param p.cutoff Filter genes by no more than Y adjusted p-value, default: 0.05. 
#' @param only.pos Only return positive genes in filtered results \code{res.sig.ls}, default: TRUE.
#' 
#' @return List containing differential analysis object, result table and filtered result table.  
#' @export
#'
#' @import edgeR
#' @import dplyr
#' @importFrom stats model.matrix
#' @importFrom limma makeContrasts
#' @importFrom rlang .data
#' @importFrom tibble rownames_to_column
edgeRDE <- function(counts, 
                    group,  
                    norm.factors = NULL, 
                    adjust.factors = NULL, 
                    design.formula = NULL, 
                    contrast.df = NULL,
                    coef = NULL,
                    logfc.cutoff = 1, p.cutoff = 0.05, only.pos = TRUE) {
  
  degs <- edgeR::DGEList(counts, group = group)
  
  if (is.null(norm.factors)) {
    degs <- edgeR::calcNormFactors(degs, method = "RLE") # Default: perform RLE normalization
  } else {
    degs$samples$norm.factors <- norm.factors
  }
  
  if (is.null(adjust.factors)) {
    design.df <- data.frame(condition = group, row.names = colnames(counts))
    design.mat <- stats::model.matrix(design.formula, data = design.df)
  } else {
    design.df <- data.frame(condition = group, adjust.factors, row.names = colnames(counts))
    design.mat <- stats::model.matrix(as.formula(paste("~0+", paste(colnames(design.df), collapse = '+'))), data = design.df)
  }
  
  degs <- edgeR::estimateDisp(degs, design = design.mat)
  fit.glm <- edgeR::glmFit(degs, design = design.mat)
  
  # extract DE results
  if (is.null(coef) & !is.null(contrast.df)) {
    contrast.vec <- apply(contrast.df, 1, function(x) { paste(paste0('condition',x), collapse='-') })
    contrast.mat <- limma::makeContrasts(contrasts = contrast.vec, levels = design.mat)
    lrt.ls <- apply(contrast.mat, 2, function(x) { edgeR::glmLRT(fit.glm, contrast = x) })  
    names(lrt.ls) <- gsub('-','_',gsub('condition','',contrast.vec))
  } 
  if (!is.null(coef) & is.null(contrast.df)) { 
    lrt.ls <- list(edgeR::glmLRT(fit.glm, coef = coef)) 
    # names(lrt.ls) <- gsub('condition','',paste(colnames(design.mat)[coef], collapse = '_'))
  }
  
  res.ls <- lapply(lrt.ls, function(x) {
    res1 <- edgeR::topTags(x, n = Inf, adjust.method = "BH")
    res.tab <- res1$table %>% tibble::rownames_to_column() %>% dplyr::rename(GeneID = rowname)
    return(res.tab)
    })
  
  # names(res.ls) <- gsub('condition', '', contrast.vec)
  # significant DEGs
  if (only.pos) {
    # return only positive regulated genes
    res.sig.ls <- lapply(res.ls, function(x) { x[x$logFC >= logfc.cutoff & x$FDR < p.cutoff,] })
  } else {
    # return both positive and negative regulated genes
    res.sig.ls <- lapply(res.ls, function(x) { x[abs(x$logFC) >= logfc.cutoff & x$FDR < p.cutoff,] })
  }
  
  return(list(de.obj = degs, res.ls = res.ls, res.sig.ls = res.sig.ls))
}

#' Calculate fold-change of synthetic RNA 
#'
#' @param dat.norm.ls List containing normalized counts and adjust factors for 
#' adjusting unwanted variation. 
#' @param syn.id Vector of synthetic RNA ids.
#' @param enrich.idx Matrix with two rows indicating the column index of 
#' enrichment and input samples in the raw/normalized count data matrix. 
#' The first row is the column index of input and the second row is the 
#' column index of enrichment samples. 
#'
#' @return Data frame with fold-change of synthetic RNA
#' @export
#' 
#' @importFrom stringr str_split
SynFC <- function(dat.norm.ls, syn.id, enrich.idx) {
  
  syn.fc.df <- data.frame(row.names = syn.id)
  for (i in names(dat.norm.ls)) {
    data.normi <- as.matrix(dat.norm.ls[[i]]$dataNorm)
    fci <- data.normi[syn.id, enrich.idx[2,]] / data.normi[syn.id, enrich.idx[1,]]
    colnames(fci) <- paste(colnames(fci), i, sep = '.')
    syn.fc.df <- cbind(syn.fc.df, fci)
  }
  syn.fc.df <- as.data.frame(t(syn.fc.df))
  syn.fc.df$id <- stringr::str_split(rownames(syn.fc.df), pattern = '\\.', simplify = TRUE)[,1]
  syn.fc.df$method <- stringr::str_split(rownames(syn.fc.df), pattern = '\\.', simplify = TRUE)[,2]
  syn.fc.df$method <- factor(syn.fc.df$method, levels = unique(syn.fc.df$method))
  return(syn.fc.df)
}

#' Dot-plot with mean_sd bar
#'
#' @param data A data.frame (or a tibble).
#' @param x The grouping variable from the \code{data}.
#' @param y The value variable from the \code{data}.
#' @param fill The fill variable from the \code{data}.  
#' @param palette The fill palette for different groups.
#'
#' @return ggplot2 object
#' @export
#' 
#' @import ggplot2
#' @importFrom paintingr paint_palette
ggDotPlot <- function(data, x, y, fill, palette = NULL) {
  
  if (is.null(palette)) {
    palette <- paintingr::paint_palette("Splash",length(unique(data[,fill])),'continuous')
  }
  
  ggplot(data, aes_string(x, y, fill = fill)) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', color = NA, 
                 dotsize = 0.8, position = 'dodge') +
    stat_summary(fun.data = .mean_sd, size = 0.5, shape = 19, 
                 position = position_dodge(width=0.9), show.legend = FALSE) +
    theme_minimal() +
    scale_fill_manual(values = palette) +
    labs(x='', y='Fold Change')
}

#' Statistics summary (mean and +/- sd)
#'
#' @param x Value
#'
#' @return Vector of mean and mean +/- sd 
#'
#' @importFrom stats sd
.mean_sd <- function(x) {
  m <- mean(x)
  ymin <- m - stats::sd(x)
  ymax <- m + stats::sd(x)
  return(c(y=m, ymin=ymin, ymax=ymax))
}

# For adjusting no visible binding
## reduceRes
utils::globalVariables(c("GeneID", "Group"))
## ggPCA
utils::globalVariables(c("PC1", "PC2", "group"))
## ggPCA_Biplot
utils::globalVariables(c("Performance", "method.id"))
## edgeRDE
utils::globalVariables(c('rowname'))