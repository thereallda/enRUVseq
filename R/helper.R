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
#' @return vector of similarity
#' @export
#'
setSimilarity <- function(set.ls, pair=FALSE) {
  cor.mat <- matrix(rep(0, length(set.ls)^2), ncol=length(set.ls))
  for (i in seq_along(set.ls)) {
    cor.mat[i,i] <- 1
    if (i == length(set.ls)) {
      break
    }
    for (j in (i+1):length(set.ls)) {
      j.idx <- JacIdx(rownames(set.ls[[i]]), rownames(set.ls[[j]]))
      cor.mat[i,j] <- j.idx
      cor.mat[j,i] <- j.idx
    }
  }
  
  if (pair) {
    pairwise.similarity <- cor.mat
  } else {
    cor.vec <- unique(as.vector(cor.mat))
    cor.vec <- sort(cor.vec, decreasing = TRUE)
    overall.similarity <- mean(cor.vec[-1])  
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
#' @param object A count matrix
#' @param labels character vector of sample names or labels. Defaults to colnames(object).
#' @param vst.norm if TRUE perform vst transformation.
#'
#' @return ggplot2 object
#' @export 
#' 
#' @importFrom DESeq2 vst
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @importFrom stats prcomp
#' @importFrom stringr str_remove
ggPCA <- function(object, labels = colnames(object), vst.norm=TRUE) {
  if (vst.norm) {
    counts_norm <- DESeq2::vst(as.matrix(object))
  } else {
    counts_norm <- object
  }
  
  pca <- prcomp(t(counts_norm))
  rownames(pca$x) <- labels
  
  pc.var <- round(summary(pca)$importance[2,], 3)
  
  as.data.frame(pca$x) %>%
    mutate(group = str_remove(labels, '\\.\\d$')) %>%
    ggplot(aes(x=PC1, y=PC2)) +
    geom_point(aes(color=group), size=3) +
    ggrepel::geom_text_repel(label=labels, max.overlaps = 20) +
    geom_vline(xintercept=0, color='grey80', lty=2) +
    geom_hline(yintercept=0, color='grey80', lty=2) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          legend.position = 'top', 
          axis.text = element_text(color='black')) +
    labs(x=paste0('PC1: ', pc.var[1]*100, '%'),
         y=paste0('PC2: ', pc.var[2]*100, '%'))
}

#' Combine list of DE results
#'
#' @param res.ls Named list of DE results. Each elements in the list correspond to a table of DE results between groups.
#' @param fc.col Name of the fold change column.
#' @param levels Factor levels of the DE group.
#'
#' @return data.frame
#' @export
#'
#' @import dplyr
reduceRes <- function(res.ls, fc.col, levels=names(res.ls)) {
  df <- data.frame()
  for (id in names(res.ls)) {
    curr <- res.ls[[grep(id, names(res.ls), value=T)]] 
    curr$GeneID <- rownames(curr)
    df1 <- curr %>% 
      dplyr::mutate(Group = factor(rep(id, nrow(curr)), levels = levels)) %>% 
      dplyr::select(GeneID, !!sym(fc.col), Group)
    df <- rbind(df, df1)
  }
  return(df)
}

#' Box-violin plot comparing values between groups
#'
#' @param data A dataframe (or a tibble).
#' @param x The grouping variable from the \code{data}.
#' @param y The value variable from the \code{data}.
#' @param color The color variable from the \code{data}.  
#' @param palette The color palette for different groups.
#'
#' @return ggplot2 object
#' @export
#' 
#' @import dplyr
#' @import ggplot2 
#' @import rstatix 
#' @import paintingr 
#' @import ggpubr
#' @importFrom stats as.formula
BetweenStatPlot <- function(data, x, y, color, palette = NULL) {
  stat.formula <- as.formula(paste(y, "~", x))
  stat_dat <- data %>% 
    wilcox_test(stat.formula) %>% 
    adjust_pvalue() %>%
    p_format(p.adj, digits = 2, leading.zero = FALSE, add.p = T, accuracy = 2e-16) %>% 
    add_xy_position(x = x, dodge=0.8, step.increase=0.5) 
  
  x.labs <- paste0(unique(data[,x]), "\n(n=", tabulate(data[,x]),")")
  
  if (is.null(palette)) palette <- paint_palette("Spring")
  
  data %>% 
    ggplot(aes_string(x, y, color = color)) +
    geom_violin(width = 0.8) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.text = element_text(color='black')) +
    scale_color_manual(values = palette) +
    scale_x_discrete(labels = x.labs) +
    stat_pvalue_manual(data = stat_dat, label = "p.adj", tip.length = 0.01, size = 3)+
    labs(x='')  
}

