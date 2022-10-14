#' @rdname createEnone
#' @export
#'
setClass(
  Class = "Enone",
  contains = "SummarizedExperiment",
  slots = list(
    counts = "list",
    enone_factor = "list",
    enone_metrics = "data.frame",
    enone_score = "data.frame",
    enrichment = "list",
    enrichment_filtered = "list",
    parameter = "list"
  )
)

# setClass(
#   Class = "Assay",
#   slots = list(
#     sample = "ANY",
#     spike_in = "ANY"
#   )
# )



#' Enone object and constructor
#' 
#' @description \code{Enone} object extends the \code{SummarizedExperiment} class. 
#' The \code{createEnone} is a easy constructor of \code{Enone} object
#' 
#' @param data A un-normalized count data matrix of shape n x p, where n is the 
#' number of samples and p is the number of features. 
#' @param bio.group Vector of samples group, e.g., c("Young.Input","Young.Enrich","Old.Input","Old.Enrich").
#' @param enrich.group Vector of enrichment group, e.g., c("Input","Enrich","Input","Enrich").
#' @param batch.group Vector of samples batch, e.g., c("A","A","B","B"), default=NULL. 
#' @param spike.in.prefix A character specify the prefix of spike-in id, e.g., "FB" stands for fly spike-in id. 
#' @param input.id Input library id, must be consistent with \code{enrich.group}, e.g., "Input". 
#' @param enrich.id Enrich library id, must be consistent with \code{enrich.group}, e.g., "Enrich".  
#' @param synthetic.id Vector of synthetic RNA id, e.g. c("Syn1","Syn2"), default=NULL. 
#'
#' Description of each slot:
#' \code{assay} \code{SummarizedExperiment::Assays} object, contains all counts. 
#' \code{counts} list. 
#' \code{enone_factor} list. 
#' \code{enone_metrics} data.frame. 
#' \code{enone_score} data.frame. 
#' \code{enrichment} list. 
#' \code{enrichment_filtered} list. 
#' \code{parameter} list.
#' 
#' @return A Enone S4 object
#' @export
#'
createEnone <- function(data, bio.group, enrich.group, batch.group=NULL,
                        spike.in.prefix, input.id="Input", enrich.id="Enrich",
                        synthetic.id=NULL
                        ) {
  
  # parameters 
  params <- list(
    spike.in.prefix = spike.in.prefix,
    synthetic.id = synthetic.id,
    input.id = input.id,
    enrich.id = enrich.id
  )
  
  # create Assay object
  # ## assay contains Assay with raw counts of sample (and spike_in)
  # counts_sample <- data[grep(spike.in.prefix, rownames(data), invert = TRUE),]
  # counts_spike_in <- data[grep(spike.in.prefix, rownames(data)),]
  # assay <- createAssay(sample = list("Raw"=counts_sample),
  #                      spike_in = list("Raw"=counts_spike_in))
  # assay <- list(sample = list("Raw"=counts_sample),
  #               spike_in = list("Raw"=counts_spike_in))
  assay <- list(sample = list(), spike_in = list())
  ## enrichment and enrichment filtered contains empty Assay
  # enrichment_assay <- enrichment_filtered_assay <- createAssay()
  enrichment_assay <- enrichment_filtered_assay <- list(sample=list(), spike_in=list())
  
  # factor slot
  enone_factor <- list(sample=list(), spike_in=list())
  
  # SummarizedExperiment object
  ## rowData for mapping gene id 
  rowDf <- S4Vectors::DataFrame(GeneID=rownames(data),
                                SpikeIn=(rownames(data) %in% grep(spike.in.prefix,rownames(data),value=TRUE)) 
  )
  # if synthetic RNA id provided
  if (!is.null(synthetic.id)) rowDf$Synthetic <- rowDf$GeneID %in% grep(paste(synthetic.id,collapse = "|"),rowDf$GeneID,value=TRUE)
  ## colData for mapping samples
  colDf <- S4Vectors::DataFrame(
    id = colnames(data),
    condition = bio.group,
    enrich = enrich.group,
    replicate = countReplicate(bio.group),
    batch = NA_character_
  )
  
  if (!is.null(batch.group)) colDf$batch <- batch.group
  # create SummarizedExperiment object
  se <- SummarizedExperiment::SummarizedExperiment(S4Vectors::SimpleList(as.matrix(data)),
                                                   colData = colDf,
                                                   rowData = rowDf)
  # create Enone object
  Enone <- methods::new("Enone", 
                       se,
                       counts = assay,
                       enone_factor = enone_factor,
                       enone_metrics = data.frame(),
                       enone_score = data.frame(),
                       enrichment = enrichment_assay,
                       enrichment_filtered = enrichment_filtered_assay,
                       parameter = params
                       )
  # validObject(Enone)
  return(Enone)
}

# createAssay <- function(sample=NULL, spike_in=NULL) {
#   
#   if (is.null(sample)) sample <- list()
#   if (is.null(spike_in)) spike_in <- list()
#   
#   new("Assay", sample=sample, spike_in=spike_in)
# }

