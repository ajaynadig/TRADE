#' TRADE Function
#'
#' This function performs univariate or bivariate TRADE analysis on differential expression (DE) summary statistics.
#'
#' @param mode The mode of analysis, either "univariate" or "bivariate".
#' @param results1 The first set of DE summary statistics.
#' @param results2 The second set of DE summary statistics (required only for bivariate mode).
#' @param annot_table Table of gene-wise annotations for enrichment analysis.
#' @param log2FoldChange The name of the column with log2FoldChanges in the DE summary statistics.
#' @param lfcSE The name of the column with log2FoldChange standard errors in the DE summary statistics.
#' @param pvalue The name of the column with p values in the DE summary statistics.
#' @param model_significant Should TRADE run the model separately on significant genes and calculate e.g. fraction significant? (default is TRUE)
#' @param genes_exclude Which genes should be excluded (e.g. perturbed genes).
#' @param estimate_sampling_covariance In cases of shared samples in bivariate analysis, use ash utility to estimate sampling covariance matrices (default is FALSE).
#' @param covariance_matrix_set The basis set of covariance matrices for bivariate analyses, either "mash_default", "adaptive_grid", or "combined" (default is "combined").
#' @param component_varexplained_threshold For adaptive grid, the threshold variance explained for TRADE to retain in analysis (default is 0).
#' @param weight_nocorr The prior on bivariate component with 0 correlation. 1 corresponds to no penalty (default is 1).
#' @param n_sample How many samples to draw from distribution and include in output.
#'
#' @return The output of the TRADE analysis.
#'
#' @examples
#' # Example usage of TRADE function
#' results1 <- read.csv("results1.csv")
#' results2 <- read.csv("results2.csv")
#' annot_table <- read.csv("annot_table.csv")
#' trade_output <- TRADE(mode = "bivariate", results1 = results1, results2 = results2, annot_table = annot_table)
#'
#' @export
TRADE <- function(mode = NULL, 
                  results1, 
                  results2 = NULL,
                  annot_table = NULL,
                  log2FoldChange = "log2FoldChange",
                  lfcSE = "lfcSE", 
                  pvalue = "pvalue",
                  model_significant = TRUE,
                  genes_exclude = NULL, 
                  estimate_sampling_covariance = FALSE,
                  covariance_matrix_set = "combined", 
                  component_varexplained_threshold = 0,
                  weight_nocorr = 1,
                  n_sample = NULL) { 

  if (!(mode %in% c("univariate","bivariate"))) {
    stop("mode improperly specified; please use 'univariate' or 'bivariate")
  } else if (mode == "bivariate" & is.null(results2)) {
    stop("for bivariate mode, please provide second set of DE results")
  }

  # Wrangle first set of input DE summary statistics
  names(results1)[names(results1) == log2FoldChange] <- "log2FoldChange"
  names(results1)[names(results1) == lfcSE] <- "lfcSE"
  if (!is.null(pvalue)) {
    names(results1)[names(results1) == pvalue] <- "pvalue"
  }

  # Wrangle second set of input DE summary statistics if in bivariate mode
  if (!is.null(results2)) {
    names(results2)[names(results2) == log2FoldChange] <- "log2FoldChange"
    names(results2)[names(results2) == lfcSE] <- "lfcSE"
    if (!is.null(pvalue)) {
      names(results2)[names(results2) == pvalue] <- "pvalue"
    }
  }

  if(!is.null(annot_table)){
    ncols_annot <- ncol(annot_table)
    annot_table <- annot_table[,colSums(annot_table) > 0]
    if(ncols_annot != ncol(annot_table)){
      message(paste(ncols_annot - ncol(annot_table) ,"columns in annotation table have no annotations and were removed for having no genes"))
    }
  }

  if (mode == "univariate") {
    output = TRADE_univariate(results = results1,
                              annot_table = annot_table,
                              model_significant = model_significant,
                              n_sample = n_sample,
                              genes_exclude = genes_exclude)

  } else if (mode == "bivariate") {
    output = TRADE_bivariate(results1 = results1,
                             results2 = results2,
                             genes_exclude = genes_exclude,
                             estimate_sampling_covariance = estimate_sampling_covariance,
                             covariance_matrix_set = covariance_matrix_set,
                             component_varexplained_threshold = component_varexplained_threshold,
                             weight_nocorr = weight_nocorr,
                             n_sample = n_sample)
  }

  return(output)



}
