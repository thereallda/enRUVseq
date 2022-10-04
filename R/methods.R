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
          function(object, slot, method) {
            object@counts[[slot]][[method]]
            })

#' @rdname Counts
#' @name Counts
#' @export "Counts<-"
setReplaceMethod("Counts", signature = signature(object="Enone", slot="character", method="character", value="matrix"),
                 function(object, slot, method, value) {
                   object@counts[[slot]][[method]] <- value
                   methods::validObject(object)
                   return(object)
                 })

#' Accessor of enONE normalization factors
#'
#' @param object Enone. 
#' @param slot Which slot to get, one of \code{sample} or \code{spike_in}.  
#' @param method Which counts matrix to get, must be one of the raw or normalized counts matrix presented in the selected slot. 
#' @name getFactor
#' @aliases getFactor getFactor,Enone,character,character-method
#'
#' @return vector or list of factors. 
#' @export
#'
setMethod("getFactor", signature = signature(object="Enone", slot="character", method="character"),
          function(object, slot, method) {
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
