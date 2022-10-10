#' Accessor for the 'counts' slot of Enone object. 
#'
#' @param object Enone. 
#' @param slot Which slot to get, one of \code{sample} or \code{spike_in}.  
#' @param method Which counts matrix to get, must be one of the raw or normalized counts matrix presented in the selected slot. 
#' @param value Raw or normalized counts matrix. 
#' @name Counts
#' @aliases Counts Counts,Enone,character,character-method 
#' Counts<-,Enone,character,character,matrix-method
#' 
#' @return matrix. 
#' @export
#'
setMethod("Counts", signature = signature(object="Enone", slot="character", method="character"), 
          function(object, slot=c("sample","spike_in"), method) {
            
            slot <- match.arg(slot, choices = c("sample","spike_in"))
            
            if (is.null(names(object@counts[[slot]]))) {
              stop("Normalizations for ", slot, " not found. At least one normalization should be performed.")
            }
            method <- match.arg(method, choices = names(object@counts[[slot]]))
            
            object@counts[[slot]][[method]]
            })

#' @rdname Counts
#' @name Counts
#' @export "Counts<-"
setReplaceMethod("Counts", signature = signature(object="Enone", slot="character", method="character", value="matrix"),
                 function(object, slot=c("sample","spike_in"), method, value) {
                   slot <- match.arg(slot, choices = c("sample","spike_in"))
                   object@counts[[slot]][[method]] <- value
                   methods::validObject(object)
                   return(object)
                 })

#' Accessor of enONE normalization factors
#'
#' @param object Enone. 
#' @param slot Which slot to get, one of \code{sample} or \code{spike_in}.  
#' @param method Which normalization methods to get, must be one of the methods presented in the selected slot. 
#' @name getFactor
#' @aliases getFactor getFactor,Enone,character,character-method
#'
#' @return vector or list of factors. 
#' @export
#'
setMethod("getFactor", signature = signature(object="Enone", slot="character", method="character"),
          function(object, slot=c("sample","spike_in"), method) {
            slot <- match.arg(slot, choices = c("sample","spike_in"))
            object@enone_factor[[slot]][[method]]
          })

#' Accessor of enONE metrics
#'
#' @param object Enone. 
#' @name getMetrics
#' @aliases getMetrics getMetrics,Enone-method
#'
#' @return data.frame
#' @export
#'
setMethod("getMetrics", signature = signature(object="Enone"),
          function(object) {
            object@enone_metrics
          })

#' Accessor of enONE score
#'
#' @param object Enone. 
#' @name getScore
#' @aliases getScore getScore,Enone-method
#' 
#' @return data.frame
#' @export
#'
setMethod("getScore", signature = signature(object="Enone"),
          function(object) {
            object@enone_score
          })

#' Find enriched genes between enrich and input samples
#'
#' @param object Enone.
#' @param slot Which slot to get, one of \code{sample} or \code{spike_in}.  
#' @param method Which normalization methods to get, must be one of the methods presented in the selected slot.  
#' @param fc.cutoff Integer indicate the cutoff for log2-Fold-Change, only 
#' consider genes with log2-Fold-Change greater than this cutoff, default: 1. 
#' @param p.cutoff Numeric indicate the cutoff for adjusted p-value (e.g., FDR), 
#' genes with adjusted p-value smaller than cutoff were consider, default: 0.05. 
#' 
#' @name FindEnrichment 
#' @aliases FindEnrichment  FindEnrichment,Enone,character,character,numeric,numeric-method
#' 
#' @return data.frame
#' @export
#'
setMethod("FindEnrichment", 
          signature = signature(object="Enone", slot="character", method="character", fc.cutoff="numeric", p.cutoff="numeric"),
          function(object, slot=c("sample","spike_in"), method, fc.cutoff=1, p.cutoff=0.05) {
            
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
            } else {
              counts_df <- SummarizedExperiment::assay(object)[!SummarizedExperiment::rowData(object)$SpikeIn & !SummarizedExperiment::rowData(object)$Synthetic,]
            }
            
            # get list of factors 
            factor.ls <- getFactor(object, slot=slot, method=method)
            if ("normFactor" %in% names(factor.ls)) {
              if ("adjustFactor" %in% names(factor.ls)) {
                # if norm factors and adjust fators were both provided
                de <- edgeRDE(counts = counts_df,
                              group = object$condition,
                              contrast.df = contrast_df,
                              norm.factors = factor.ls$normFactor,
                              adjust.factors = factor.ls$adjustFactor,
                              fc.cutoff=fc.cutoff, p.cutoff=p.cutoff
                              )
              } else {
                # if only norm factors were provided
                de <- edgeRDE(counts = counts_df,
                              group = object$condition,
                              contrast.df = contrast_df,
                              norm.factors = factor.ls$normFactor,
                              design.formula = "~0+condition",
                              fc.cutoff=fc.cutoff, p.cutoff=p.cutoff)
              }
            } else {
              stop("One or both of 'normFactor' and 'adjustFacotr' should be provided.")
            }
            # save enrichment in Enone object
            object@enrichment[[slot]] <- de$res.ls
            object@enrichment_filtered[[slot]] <- de$res.sig.ls
            
            return(object)
            })
