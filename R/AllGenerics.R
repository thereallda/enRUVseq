#' @rdname Counts
#' @export
setGeneric("Counts", function(object, slot, method) standardGeneric("Counts"))

#' @rdname Counts
#' @export
setGeneric("Counts<-", function(object, slot, method, value) standardGeneric("Counts<-"))

#' @rdname getFactor
#' @export
setGeneric("getFactor", function(object, slot, method) standardGeneric("getFactor"))

#' @rdname getMetrics
#' @export
setGeneric("getMetrics", function(object) standardGeneric("getMetrics"))

#' @rdname getScore
#' @export
setGeneric("getScore", function(object) standardGeneric("getScore"))

#' @rdname FindEnrichment
#' @export
setGeneric("FindEnrichment", function(object, slot, method, fc.cutoff, p.cutoff) standardGeneric("FindEnrichment"))
