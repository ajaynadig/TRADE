require(ggplot2)
require(doBy)
require(ashr)

TRADE_univariate <- function(results = NULL,
                           annot_table = NULL,
                           model_significant = TRUE,
                           n_sample = 10000,
                           genes_exclude = NULL,
                           verbose = FALSE) {
  
  if (model_significant) {
    na.filter = !is.finite(results$log2FoldChange) | !is.finite(results$lfcSE) | !is.finite(results$pvalue)
  } else {
    na.filter = !is.finite(results$log2FoldChange) | !is.finite(results$lfcSE)
    
  }
  num_na = sum(na.filter)

  extreme_log2FC = abs(results$log2FoldChange) > 10
  num_extreme = sum(extreme_log2FC)

  num_exclude = sum(na.filter | abs(results$log2FoldChange) > 10)

  results = results[!na.filter & !extreme_log2FC,]


  message(paste0(num_na, " genes excluded with NA entries"))
  message(paste0(num_extreme, " genes excluded with extreme log2FC (abs value > 10)"))
  message(paste0(nrow(results), " genes retained for analysis"))

  #Gene Exclusion
  if (!is.null(genes_exclude) & any(genes_exclude %in% rownames(results))){

    genes_exclude = genes_exclude[genes_exclude %in% rownames(results)]
    results = results[!(rownames(results) %in% genes_exclude),]
    annot_table = annot_table[!(rownames(annot_table) %in% genes_exclude),]
  } else if (!is.null(genes_exclude) & !(any(genes_exclude %in% rownames(results)))) {
    message("No Exclusion Genes in Results")
  }

  #Set number of samples to be drawn to estimate distribution features
  if (is.null(n_sample)) {n_sample = nrow(results)}
  message("Fitting mixture model")
  main_output = get_distribution_output(results,n_sample, annot_table = annot_table, model_significant = model_significant, verbose = verbose)

  main_output$qc = list(num_na = num_na,
                        num_extreme = num_extreme,
                        num_exclude = num_exclude)
  return(main_output)
}

fit_ash <- function(l2fc,l2fc_se,min_l2fc,max_l2fc,n_sample) {
  ash_model = ash(betahat = l2fc,
                  sebetahat = l2fc_se,
                  mixcompdist = "halfuniform",
                  outputlevel = 3,
                  grange = c(min_l2fc,max_l2fc),
                  prior = "uniform")

  weights =  ash_model$fitted_g$pi
  uniform_a = ash_model$fitted_g$a
  uniform_b = ash_model$fitted_g$b


  component = sample(c(1:length(weights)),size = n_sample, prob = weights, replace = TRUE)
  samples = runif(n_sample, min = uniform_a[component], max = uniform_b[component])


  return(list(fit = ash_model,
              samples = samples))
}

get_distribution_output <- function(results,n_sample, annot_table = NULL, model_significant = TRUE, verbose = FALSE) {
  min_effsize = min(results$log2FoldChange, na.rm = TRUE)
  max_effsize = max(results$log2FoldChange, na.rm = TRUE)

  total_output = fit_ash(results$log2FoldChange,
                         results$lfcSE,
                         min_effsize,
                         max_effsize,
                         n_sample)
  
  if (verbose) {message("Computing Distribution Features")}
  

  #calculate var and mean
  means = (total_output$fit$fitted_g$a + total_output$fit$fitted_g$b)/2
  vars = (1/12) * (total_output$fit$fitted_g$b - total_output$fit$fitted_g$a)^2

  mixture_mean = sum(total_output$fit$fitted_g$pi * means)

  variance_expectation = sum(total_output$fit$fitted_g$pi * (means - mixture_mean)^2)
  expectation_variance =  sum(total_output$fit$fitted_g$pi * vars)

  mixture_variance = variance_expectation + expectation_variance

  #calculate Me
  fourthmoment_numerator = (((total_output$fit$fitted_g$b - mixture_mean)^5) - ((total_output$fit$fitted_g$a - mixture_mean)^5))
  fourthmoment_denominator = (total_output$fit$fitted_g$b - mixture_mean - total_output$fit$fitted_g$a + mixture_mean)
  fourthmoments = (1/5) * (fourthmoment_numerator/fourthmoment_denominator)

  #Set fourth moment 0 for point mass
  fourthmoments[1] <- 0

  kappa_numerator = sum(total_output$fit$fitted_g$pi * fourthmoments)
  kappa_denominator = mixture_variance^2

  kappa = kappa_numerator/kappa_denominator

  Me = 3 * nrow(results)/kappa

  df_viz = data.frame(log2FoldChange = c(results$log2FoldChange,total_output$samples),
                      distribution = c(rep("Empirical",nrow(results)),
                                       rep("Inferred",length(total_output$samples))))
  plot <- ggplot(data = df_viz,
                 mapping = aes(x = log2FoldChange,
                               fill = distribution,
                               y = after_stat(scaled)))+
    geom_density(alpha = 0.2)+
    theme_bw()+
    labs(x = "Log2(FC)")

  if (model_significant) {
    if (verbose) {message("Computing Significant Gene Features")}
    
    #significant genes: Bonferroni
    sig_filter_Bonferroni = results$pvalue < 0.05/sum(!is.na(results$pvalue))
    sig_filter_Bonferroni[is.na(sig_filter_Bonferroni)] <- FALSE
    if (any(sig_filter_Bonferroni)){

      sig_tim_Bonferroni <- frac_subset(results,sig_filter_Bonferroni,!sig_filter_Bonferroni,min_effsize,max_effsize,n_sample)
      var_sig_Bonferroni = sig_tim_Bonferroni$var_subset
      var_nonsig_Bonferroni = sig_tim_Bonferroni$var_complement
      frac_sig_Bonferroni = sig_tim_Bonferroni$frac_subset
      sig_results_Bonferroni = results[sig_filter_Bonferroni,]
      num_sig_Bonferroni = sig_tim_Bonferroni$n_subset
      num_nonsig_Bonferroni = sig_tim_Bonferroni$n_complement

    } else {
      var_sig_Bonferroni = 0
      frac_sig_Bonferroni = 0
      sig_results_Bonferroni = NA
      num_sig_Bonferroni = 0
      var_nonsig_Bonferroni = NA
      num_nonsig_Bonferroni = nrow(results)
    }
    
    #significant genes: FDR
    sig_filter_FDR = p.adjust(results$pvalue, method = "fdr") < 0.05
    sig_filter_FDR[is.na(sig_filter_FDR)] <- FALSE
    if (any(sig_filter_FDR)){

      sig_tim_FDR <- frac_subset(results,sig_filter_FDR,!sig_filter_FDR,min_effsize,max_effsize,n_sample)
      var_sig_FDR = sig_tim_FDR$var_subset
      var_nonsig_FDR = sig_tim_FDR$var_complement
      frac_sig_FDR = sig_tim_FDR$frac_subset
      sig_results_FDR = results[sig_filter_FDR,]
      num_sig_FDR = sig_tim_FDR$n_subset
      num_nonsig_FDR = sig_tim_FDR$n_complement

    } else {
      var_sig_FDR = 0
      frac_sig_FDR = 0
      sig_results_FDR = NA
      num_sig_FDR = 0
      var_nonsig_FDR = NA
      num_nonsig_FDR = nrow(results)
    }
  } else {
    var_sig_Bonferroni = NA
    frac_sig_Bonferroni = NA
    sig_results_Bonferroni = NA
    num_sig_Bonferroni = NA
    var_nonsig_Bonferroni = NA
    num_nonsig_Bonferroni = NA
    var_sig_FDR = NA
    frac_sig_FDR = NA
    sig_results_FDR = NA
    num_sig_FDR = NA
    var_nonsig_FDR = NA
    num_nonsig_FDR = NA
  }


#  gene-set enrichments
  if (!is.null(annot_table)) {
    if (verbose) {message("Computing Enrichments")}
    
    frac_var_annot = sapply(1:ncol(annot_table),
                        function(annot) {
                          output_annot =  frac_subset(results,
                                                      annot_table[,annot] == 1,
                                                      annot_table[,annot] == 0,
                                             min_effsize,
                                             max_effsize,
                                             n_sample)
                          return(list(var_sig = output_annot$var_subset,
                                   frac_sig = output_annot$frac_subset))
                        })

    frac_genes_annot = sapply(1:ncol(annot_table),
                              function(annot) {
                                frac = sum(annot_table[,annot],na.rm = TRUE)/sum(!is.na(annot_table[,annot]))
                              })
    annot_output_table = data.frame(annot = colnames(annot_table),
                                    var = unlist(frac_var_annot[1,]),
                                    frac_var = unlist(frac_var_annot[2,]),
                                    frac_genes = frac_genes_annot,
                                    enrichment = unlist(frac_var_annot[2,])/frac_genes_annot)

  } else {
    annot_output_table = NA
  }
  
  if (verbose) {message(paste("TI =",round(mixture_variance,4)))}
  

  return(list(plot = plot,
              fit = list(distribution = total_output$fit$fitted_g,
                         #posteriormeans = total_output$fit$result$PosteriorMean,
                         loglik = total_output$fit$loglik),
              distribution_summary = list(mean = mixture_mean,
                                  transcriptome_wide_impact = mixture_variance,
                                  Me = Me),
              significant_genes_Bonferroni = list(significant_gene_results_Bonferroni = sig_results_Bonferroni,
                                       var_sig_Bonferroni = var_sig_Bonferroni,
                                       var_nonsig_Bonferroni = var_nonsig_Bonferroni,
                                       frac_sig_Bonferroni = frac_sig_Bonferroni,
                                       num_sig_Bonferroni = num_sig_Bonferroni,
                                       num_nonsig_Bonferroni = num_nonsig_Bonferroni),
              significant_genes_FDR = list(significant_gene_results_FDR = sig_results_FDR,
                                                  var_sig_FDR = var_sig_FDR,
                                           var_nonsig_FDR = var_nonsig_FDR,
                                                  frac_sig_FDR = frac_sig_FDR,
                                                  num_sig_FDR = num_sig_FDR,
                                           num_nonsig_FDR = num_nonsig_FDR),
              samples = total_output$samples,
              annot_output = annot_output_table))
}


frac_subset <- function(results,filter,filter_complement,min_effsize,max_effsize,n_sample) {
  subset_output = fit_ash(results$log2FoldChange[filter],
                       results$lfcSE[filter],
                       min_effsize,
                       max_effsize,
                       n_sample)

  complement_output = fit_ash(results$log2FoldChange[filter_complement],
                          results$lfcSE[filter_complement],
                          min_effsize,
                          max_effsize,
                          n_sample)

  means_subset = (subset_output$fit$fitted_g$a + subset_output$fit$fitted_g$b)/2
  vars_subset = (1/12) * (subset_output$fit$fitted_g$b - subset_output$fit$fitted_g$a)^2

  mixture_mean_subset = sum(subset_output$fit$fitted_g$pi * means_subset)

  variance_expectation_subset = sum(subset_output$fit$fitted_g$pi * (means_subset - mixture_mean_subset)^2)
  expectation_variance_subset = sum(subset_output$fit$fitted_g$pi * vars_subset)
  mixture_variance_subset = variance_expectation_subset + expectation_variance_subset

  means_complement = (complement_output$fit$fitted_g$a + complement_output$fit$fitted_g$b)/2
  vars_complement = (1/12) * (complement_output$fit$fitted_g$b - complement_output$fit$fitted_g$a)^2

  mixture_mean_complement = sum(complement_output$fit$fitted_g$pi * means_complement)

  variance_expectation_complement = sum(complement_output$fit$fitted_g$pi * (means_complement - mixture_mean_complement)^2)
  expectation_variance_complement = sum(complement_output$fit$fitted_g$pi * vars_complement)
  mixture_variance_complement = variance_expectation_complement + expectation_variance_complement

  var_subset = sum(filter,na.rm = TRUE) * mixture_variance_subset
  var_complement = sum(filter_complement, na.rm = TRUE) * mixture_variance_complement

  frac_subset = var_subset/(var_subset + var_complement)

  return(list(var_subset = mixture_variance_subset,
              var_complement = mixture_variance_complement,
         frac_subset = frac_subset,
         n_subset = sum(filter),
         n_complement = sum(filter_complement)))
}




