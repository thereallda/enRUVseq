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
#' @importFrom rlang .data
BetweenStatPlot <- function(data, x, y, color, palette = NULL) {
  stat.formula <- as.formula(paste(y, "~", x))
  stat_dat <- data %>% 
    wilcox_test(stat.formula) %>% 
    adjust_pvalue() %>%
    p_format(.data$p.adj, digits = 2, leading.zero = FALSE, 
             trailing.zero = TRUE, add.p = TRUE, accuracy = 2e-16) %>% 
    add_xy_position(x = x, dodge=0.8, step.increase=0.5) 
  
  x.labs <- paste0(unique(data[,x]), "\n(n=", tabulate(data[,x]),")")
  x.num <- length(unique(data[,x])) # number of x types
  if (is.null(palette)) palette <- paint_palette("Spring", x.num, 'continuous')
  
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


#' Assessment of normalization
#'
#' @param data.raw A un-normalized count data matrix of shape n x p, where n is 
#' the number of samples and p is the number of features. 
#' @param data.normalized A normalized count data matrix (n x p). 
#' @param pos.controls Vector of positive enriched genes. 
#' @param enrich.idx Matrix with two rows indicating the column index of 
#' enrichment and input samples in the raw/normalized count data matrix. 
#' The first row is the column index of input and the second row is the 
#' column index of enrichment samples. 
#'
#' @return List of metrics that assessing normalization effect. 
#' @export
#'
#' @importFrom stats cor
AssessNormalization <- function(data.raw, 
                                data.normalized,
                                pos.controls,
                                enrich.idx
                                # de.raw.res, 
                                # de.norm.res
                                ) {
  ## Constants
  # Minimum number of controls (hard threshold)
  min.num.pos.controls <- 100
  
  # control must be present in raw data
  if (!all(pos.controls %in% rownames(data.raw))) {
    stop("Could not find all positive control features in raw data.")
  }
  # control must be present in normalized data
  if (!all(pos.controls %in% rownames(data.normalized))) {
    stop("Could not find all positive control features in normalized data.")
  }
  # number of control must greater than min.num.pos.controls
  if (!all(length(pos.controls) >= min.num.pos.controls)) {
    stop("Number of controls must have at least 100. ")
  }
  
  # perform log2 transformation with prior.count offset
  # prior.count is the average count, proportional to the library size, to be 
  # added to count data to avoid taking the log of zero
  prior.count <- mean(1e6*2/colSums(data.raw))
  data.raw.log <- log2(data.raw + 1)
  data.normalized.log <- log2(data.normalized + 1)
  
  # calculate log-Fold-Change
  lfc.raw <- data.raw.log[, enrich.idx[2,] ] - data.raw.log[, enrich.idx[1,] ]
  lfc.norm <- data.normalized.log[, enrich.idx[2,] ] - data.normalized.log[, enrich.idx[1,] ]
  
  # calculate improvement of correlation between samples among positive controls 
  # before and after normalization
  raw.cor <- stats::cor(lfc.raw[pos.controls,], method = 'pearson')
  norm.cor <- stats::cor(lfc.norm[pos.controls,], method = 'pearson')
  diff.cor <- calcCorDiff(raw.cor, norm.cor)
  
  # # based on enriched counts
  # raw.cor.en <- stats::cor(data.raw.log[pos.controls, enrich.idx[2,]])
  # norm.cor.en <- stats::cor(data.normalized.log[pos.controls, enrich.idx[2,]])
  # diff.cor.en <- calcCorDiff(raw.cor.en, norm.cor.en)
  # 
  # # based on all counts
  # raw.cor.all <- stats::cor(data.raw.log[pos.controls, ])
  # norm.cor.all <- stats::cor(data.normalized.log[pos.controls, ])
  # diff.cor.all <- calcCorDiff(raw.cor.all, norm.cor.all)
  
  # # assess differential analysis performance
  # raw.sim <- setSimilarity(de.raw.res)
  # norm.sim <- setSimilarity(de.norm.res)
  # diff.sim <- calcCorDiff(raw.sim, norm.sim)
  # 
  metrics <- list(
    DC = diff.cor
    # DC.en = diff.cor.en,
    # DC.all = diff.cor.all
    # DiffSim = diff.sim
  )
  
  return(metrics)
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
  
  if ("DESeq" %in% normalization.methods) {
    design.formula <- as.formula("~condition")
    de.ls[["DESeq"]] <- DESeq2DE(counts = data.raw, 
                                 group = group,
                                 design.formula = design.formula, 
                                 contrast.df = contrast_df)
    
    de.tmp <- DESeq2DE(counts = data.raw,
                      group = group_ei,
                      design.formula = design.formula,
                      contrast.df = contrast_ei)
    de.ls[["DESeq"]]$res.sig.ls <- c(de.tmp$res.sig.ls, de.ls[["DESeq"]]$res.sig.ls)
  }
  
  if ("RLE" %in% normalization.methods) {
    design.formula <- as.formula("~0+condition")
    de.ls[["RLE"]] <- edgeRDE(counts = data.raw,
                              group = group,
                              design.formula = design.formula,
                              contrast.df = contrast_df,
                              norm.factors = data.norm.ls[["RLE"]]$normFactor)
    
    de.tmp <- edgeRDE(counts = data.raw,
                      group = group_ei,
                      design.formula = design.formula,
                      contrast.df = contrast_ei,
                      norm.factors = data.norm.ls[["RLE"]]$normFactor)
    de.ls[["RLE"]]$res.sig.ls <- c(de.tmp$res.sig.ls, de.ls[["RLE"]]$res.sig.ls)
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
                    contrast.df) {
  
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
  contrast.vec <- apply(contrast.df, 1, function(x) { paste(paste0('condition',x), collapse='-') })
  contrast.mat <- limma::makeContrasts(contrasts = contrast.vec, levels = design.mat)
  lrt.ls <- apply(contrast.mat, 2, function(x) { edgeR::glmLRT(fit.glm, contrast = x) })
  res.ls <- lapply(lrt.ls, function(x) {
    res1 <- edgeR::topTags(x, n = Inf, adjust.method = "BH")
    res.tab <- res1$table %>% tibble::rownames_to_column() %>% dplyr::rename(GeneID = rowname)
    return(res.tab)
    })
  names(res.ls) <- gsub('condition', '', contrast.vec)
  # cutoff for significant DEGs
  res.sig.ls <- lapply(res.ls, function(x) { x[x$logFC >= 1 & x$FDR < 0.05,] })
  
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
#' @param data A dataframe (or a tibble).
#' @param x The grouping variable from the \code{data}.
#' @param y The value variable from the \code{data}.
#' @param fill The fill variable from the \code{data}.  
#' @param palette The fill palette for different groups.
#'
#' @return ggplot2 object
#' @export
#'
ggDotPlot <- function(data, x, y, fill, palette = NULL) {
  
  if (is.null(palette)) {
    palette <- paintingr::paint_palette("Splash",length(unique(data[,fill])),'continuous')
  }
  
  data %>% 
    ggplot(aes_string(x, y, fill = fill)) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', color = NA, 
                 dotsize = 0.8, position = 'dodge') +
    stat_summary(fun.data = mean_sd, size = 0.5, shape = 19, 
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
#' @export
#'
#' @importFrom stats sd
mean_sd <- function(x) {
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
## edgeR
utils::globalVariables(c('rowname'))