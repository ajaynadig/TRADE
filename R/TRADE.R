TRADE <- function(mode = NULL, #univariate or bivariate
                  results1, #first set of DE summary statistics
                  results2 = NULL, #second set of DE summary statistics
                  annot_table = NULL, #table of gene-wise annotations for enrichment analysis
                  log2FoldChange = "log2FoldChange", #name of column with log2FoldChanges
                  lfcSE = "lfcSE", #name of column with log2FoldChange standard errors
                  pvalue = "pvalue", #name of column with p values
                  model_significant = TRUE, #should trade run model separately on significant genes and calculate e.g. fraction significant?
                  genes_exclude = NULL, #which genes should be excluded (e.g. perturbed genes?)
                  estimate_sampling_covariance = FALSE, #in cases of shared samples in bivariate analysis, use ash utility to estimate sampling covariance matrices
                  covariance_matrix_set = "adaptive_grid", #what basis set of covariance matrices for bivariate analyses, "mash_default" or "adaptive_grid"
                  component_varexplained_threshold = 0, # for adaptive grid, what threshold variance explained for TRADE to retain in analysis
                  weight_nocorr = 1, # prior on bivariate component with 0 correlation. 1 corresponds to no penalty
                  n_samples = NULL) { #how many samples to draw from distribution and include in output

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

  if (mode == "univariate") {
    output = TRADE_univariate(results = results1,
                              annot_table = annot_table,
                              model_significant = model_significant,
                              n_sample = n_samples,
                              genes_exclude = genes_exclude)

  } else if (mode == "bivariate") {
    output = TRADE_bivariate(results1 = results1,
                             results2 = results2,
                             genes_exclude = NULL,
                             estimate_sampling_covariance = FALSE,
                             covariance_matrix_set = covariance_matrix_set,
                             component_varexplained_threshold = 0,
                             weight_nocorr = 1,
                             n_sample = n_samples)
  }

  return(output)



}
