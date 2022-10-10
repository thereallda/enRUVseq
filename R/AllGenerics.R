#' @rdname Counts
#' @export
setGeneric("Counts", function(object, slot=c("sample","spike_in"), method) standardGeneric("Counts"))

#' @rdname Counts
#' @export
setGeneric("Counts<-", function(object, slot=c("sample","spike_in"), method, value) standardGeneric("Counts<-"))

#' @rdname getFactor
#' @export
setGeneric("getFactor", function(object, slot=c("sample","spike_in"), method) standardGeneric("getFactor"))

#' @rdname getMetrics
#' @export
setGeneric("getMetrics", function(object) standardGeneric("getMetrics"))

#' @rdname getScore
#' @export
setGeneric("getScore", function(object) standardGeneric("getScore"))

#' @rdname FindEnrichment
#' @export
setGeneric("FindEnrichment", function(object, slot=c("sample","spike_in"), method, fc.cutoff=1, p.cutoff=0.05) standardGeneric("FindEnrichment"))
