#' Normalization and normalization assessment in one function
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is the number of samples and p is the number of features. 
#' @param group Vector of samples group, e.g., c("Young.Input","Young.Enrich")
#' @param spike.in.prefix A character specify the prefix of spike-in id. 
#' @param n.neg.control Number of negative control genes for RUV normalization. 
#' @param n.pos.eval Number of positive evaluation genes for wanted variation assessment.
#' @param n.neg.eval Number of negative evaluation genes for unwanted variation assessment.
#' @param scaling.method Vector of normalization methods that are applied to the data.
#'   Available methods are: \code{c("TC", "UQ", "TMM", "DESeq")}. 
#'   Select one or multiple methods. By default all normalization methods will be applied.
#' @param ruv.norm Whether to perform RUV normalization. 
#' @param ruv.k The number of factors of unwanted variation to be estimated from the data.
#' @param ruv.drop The number of singular values to drop in the estimation of 
#' unwanted variation, default drop the first singular value that represent the 
#' difference between enrichment and input. 
#' @param pam_krange Integer or vector of integers indicates the number of 
#' clusters for PAM clustering, default: 2:6. 
#' @param pc_k Integer indicates the metrics will be calculated in the first kth PCs, default: 3.
#'
#' @return list
#' @export
#'
#' @importFrom utils head
#' @importFrom stringr str_extract
#' @importFrom stats as.formula model.matrix
enONE <- function(data, group, spike.in.prefix, 
                  n.neg.control = 1000, n.pos.eval = 1000, n.neg.eval = 1000,
                  scaling.method = c("TC", "UQ", "TMM", "DESeq"),
                  ruv.norm = TRUE, ruv.k = 1, ruv.drop = 0,
                  pam_krange = 2:6, pc_k = 3) {
  
  sc_idx <-  t(sapply(unique(group), function(i) grep(i, group)))
  enrich_idx <- matrix(c(grep('Input', group, ignore.case = TRUE), 
                         grep('Enrich', group, ignore.case = TRUE)), 
                       nrow = 2, byrow = TRUE)

  counts_nsp <- data[grep(spike.in.prefix, rownames(data), invert = TRUE),]
  counts_sp <- data[grep(spike.in.prefix, rownames(data)),]
  ## gene selection 
  ### 1. negative control genes for RUV
  cat(paste("The number of negative control genes for RUV:",n.neg.control,"\n"))
  enrich_group <- str_extract(group, "(Input)|(Enrich)")
  designMat <- model.matrix(~0+enrich_group)
  deg.en <- edgeRDE(counts_sp,
                    group = enrich_group,
                    design.formula = as.formula("~0+condition"),
                    contrast.df = data.frame(Group1='Enrich',Group2='Input')
  )
  # non-sig de top 1000
  res_tab <- deg.en$res.ls$Enrich_Input
  neg.control <- head(res_tab[order(res_tab$FDR, decreasing = TRUE),]$GeneID, n=n.neg.control)
  
  ### 2. positive evaluation genes
  de.all <- edgeRDE(counts_nsp[!rownames(counts_nsp) %in% c('Syn1','Syn2'),],
                    group = group,
                    design.formula = as.formula("~condition"),
                    coef = 2:length(unique(group)))
  res_tab <- de.all$res.ls[[1]]
  
  cat(paste("The number of positive evaluation genes:",n.pos.eval,"\n"))
  pos.eval.set <- head(res_tab[order(res_tab$FDR),]$GeneID, n=n.pos.eval)
  ### 3. negative evaluation genes
  cat(paste("The number of negative evaluation genes:",n.neg.eval,"\n"))
  neg.eval.set <- head(res_tab[order(res_tab$FDR, decreasing = TRUE),]$GeneID, n=n.neg.eval)
  
  ## apply normalization 
  cat("Apply normalization...\n")
  norm.nsp.ls <- ApplyNormalization(data, 
                                    scaling.method = scaling.method,  
                                    ruv.norm = ruv.norm, ruv.k = ruv.k, ruv.drop = ruv.drop,
                                    spike.in.prefix = spike.in.prefix,
                                    # below parameters are created inside function
                                    control.idx = neg.control, 
                                    sc.idx = sc_idx, 
                                    enrich.idx = enrich_idx)
  ## assessment 
  bio_group_index <- as.numeric(factor(group, levels=unique(group)))
  assay_group_index <- as.numeric(factor(enrich_group, levels=unique(enrich_group)))
  cat("Perform assessment...\n")
  norm.nsp.eval <- AssessNormalization(norm.nsp.ls, 
                                       pam_krange = pam_krange,
                                       pc_k = pc_k,
                                       # below parameters are created inside function
                                       bio_group = bio_group_index, 
                                       assay_group = assay_group_index, 
                                       pos.eval.set = pos.eval.set,
                                       neg.eval.set = neg.eval.set)
  return(list(
    gene.set = list('NC'=neg.control, 'PE'=pos.eval.set, 'NE'=neg.eval.set),
    norm.data.ls = norm.nsp.ls,
    norm.assessment = norm.nsp.eval
  ))
  
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
#' @param labels character vector of sample names or labels. Defaults to colnames(object).
#' @param vst.norm if TRUE perform vst transformation.
#' @param palette The color palette for different groups.
#'
#' @return ggplot2 object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom DESeq2 vst
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats prcomp
#' @importFrom stringr str_remove
#' @importFrom paintingr paint_palette
ggPCA <- function(object, labels = colnames(object), vst.norm=FALSE, palette = NULL) {
  if (vst.norm) {
    counts_norm <- DESeq2::vst(as.matrix(object))
  } else {
    counts_norm <- object
  }

  pca <- prcomp(t(counts_norm))
  rownames(pca$x) <- labels

  pc.var <- round(summary(pca)$importance[2,], 3)

  pca_dat <- as.data.frame(pca$x) %>%
    mutate(group = str_remove(labels, '\\.\\d$'))

  if (is.null(palette)) {
    palette <- paint_palette("Spring", length(unique(pca_dat$group)), 'continuous')
  }

  pca_dat %>%
    ggplot(aes(x=PC1, y=PC2)) +
    geom_point(aes(color=group), size=3) +
    geom_text_repel(label=labels, max.overlaps = 20) +
    geom_vline(xintercept=0, color='grey80', lty=2) +
    geom_hline(yintercept=0, color='grey80', lty=2) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = 'top',
          axis.text = element_text(color='black')) +
    scale_color_manual(values = palette) +
    labs(x=paste0('PC1: ', pc.var[1]*100, '%'),
         y=paste0('PC2: ', pc.var[2]*100, '%'))
}

# wrapper of fviz_pca_biplot

#' Biplot of individuals and variables
#'
#' @param object PCA object returned by \code{prcomp}. 
#' @param performance_score Vector of performance scores from ass. 
#' @param ... Additional parameters can be passed to \code{factoextra::fviz_pca_biplot()}.
#'
#' @return ggplot2 object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom factoextra fviz_pca_biplot
#' @importFrom tibble rownames_to_column
#' @importFrom paintingr paint_palette
ggPCA_Biplot <- function(object, performance_score, ...) {
  pc.score <- as.data.frame(object$x) %>% rownames_to_column() %>% mutate(Performance=performance_score)
  
  factoextra::fviz_pca_biplot(object, geom.ind = 'text', repel = TRUE, ...) +
    geom_point(data=pc.score, aes(PC1, PC2, color=Performance), size=3) +
    scale_color_gradientn(colors = paint_palette('Vesuvius', 100, 'continuous'))
  
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
#' @param data A data frame (or a tibble).
#' @param x The grouping variable from the \code{data}.
#' @param y The value variable from the \code{data}.
#' @param color The color variable from the \code{data}.  
#' @param palette The color palette for different groups.
#' @param test Perform "wilcox.test" or "t.test" or not test.  
#'
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
BetweenStatPlot <- function(data, x, y, color, palette = NULL, test = c('wilcox.test', 't.test', 'none')) {
  stat.formula <- as.formula(paste(y, "~", x))
  
  test <- match.arg(test, choices = c('wilcox.test', 't.test', 'none'))

  if (test != 'none') {
    if (test == 'wilcox.test') {
      stat_dat <- data %>% 
        wilcox_test(stat.formula) 
    } 
    if (test == 't.test') {
      stat_dat <- data %>% 
        t_test(stat.formula) 
    }
    stat_dat <- stat_dat %>% 
    adjust_pvalue() %>%
    p_format(.data$p.adj, digits = 2, leading.zero = FALSE, 
             trailing.zero = TRUE, add.p = TRUE, accuracy = 2e-16) %>% 
    add_xy_position(x = x, dodge=0.8, step.increase=0.5) 
  }
  
  x.labs <- paste0(unique(data[,x]), "\n(n=", tabulate(data[,x]),")")
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
#' @param coef integer or character vector indicating which coefficients of the 
#' linear model are to be tested equal to zero. Values must be columns or column names of design. 
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
                    coef = NULL) {
  
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
## ggPCA_Biplot
utils::globalVariables(c("Performance"))
## edgeR
utils::globalVariables(c('rowname'))