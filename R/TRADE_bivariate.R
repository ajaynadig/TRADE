require(mashr)
require(ashr)

TRADE_bivariate <- function(results1 = NULL,
                          results2 = NULL,
                          genes_exclude = NULL,
                          estimate_sampling_covariance = FALSE,
                          covariance_matrix_set,
                          component_varexplained_threshold = 0,
                          weight_nocorr = 1,
                          n_sample = NULL) {
  start_time = Sys.time()

  #Gene Exclusion
  if (!is.null(genes_exclude) & any(genes_exclude %in% c(rownames(results1),rownames(results2)))){
    genes_exclude = genes_exclude[genes_exclude %in% c(rownames(results1),rownames(results2))]
    results1 = results1[!(rownames(results1) %in% genes_exclude),]
    results2 = results2[!(rownames(results2) %in% genes_exclude),]
  } else if (!is.null(genes_exclude) & !(any(genes_exclude %in% c(rownames(results1),rownames(results2))))) {
    message("No Exclusion Genes in Results")
  }

  message("Fitting mixture model")


  intersect_genes = Reduce(intersect,
                           list(rownames(results1),
                                rownames(results2)))

  betahat_df <- cbind(results1$log2FoldChange[match(intersect_genes,rownames(results1))],
                      results2$log2FoldChange[match(intersect_genes,rownames(results2))])

  betahat_df_complete = betahat_df[complete.cases(betahat_df),]
  cor_raw = cor(betahat_df_complete[,1],betahat_df_complete[,2])

  se_df <- cbind(results1$lfcSE[match(intersect_genes,rownames(results1))],
                 results2$lfcSE[match(intersect_genes,rownames(results2))])

  data = mash_set_data(betahat_df, se_df)

  #Generate grid of covariance matrices
  #First, run univariate TIM to help identify mixture components

  if (covariance_matrix_set == "grid") {
    slopes = logseq(0.1,10,n=5)
    corrs = seq(-1,1,length.out = 9)
    U_grid = list()
    i = 1
    for (slope in slopes) {
      for (corr in corrs) {
        U_i = diag(2)
        U_i[1,1] = 1/(slope + 1)
        U_i[2,2] = slope/(slope + 1)
        U_i[1,2] = corr * sqrt(U_i[1,1] * U_i[2,2])
        U_i[2,1] = corr * sqrt(U_i[1,1] * U_i[2,2])

        U_grid[[i]] = U_i
        i = i + 1
      }
    }

    U = U_grid

    if (estimate_sampling_covariance) {
      message("Estimate sampling correlation")
      V.em = mash_estimate_corr_em(data, U)
      data = mash_update_data(data, V=V.em$V)
    }

    m   = mash(data, U)
  } else if (covariance_matrix_set == "mash_default") {

    message("Get data-driven covariance matrices")
    m.1by1 = mash_1by1(data)
    strong = unique(which.minn(m.1by1$result$lfsr[,1],round(0.05 * nrow(m.1by1$result$lfsr))),
                    which.minn(m.1by1$result$lfsr[,2],round(0.05 * nrow((m.1by1$result$lfsr)))))
    U.pca = cov_pca(data,2,subset=strong)
    U.ed = cov_ed(data, U.pca, subset=strong)

    message("Get canonical covariance matrices")
    U.c = cov_canonical(data)


    U = c(U.c,U.ed)

    if (estimate_sampling_covariance) {
      message("Estimate sampling correlation")
      V.em = mash_estimate_corr_em(data, U)
      data = mash_update_data(data, V=V.em$V)
    }

    m   = mash(data, U)
  } else if (covariance_matrix_set == "adaptive_grid") {

    results1_ash <- ash(betahat_df[,1],
                        se_df[,1],
                        mixcompdist = "halfnormal",
                        prior = "uniform")

    results1_sds = unique(results1_ash$fitted_g$sd[-1])
    results1_sds_weights <- sapply(results1_sds, function(x) sum(results1_ash$fitted_g$pi[results1_ash$fitted_g$sd == x]))
    #results1_sds_varexplained = (results1_sds_weights * sum(complete.cases(betahat_df))) *(results1_sds^2)
    #results1_sds_proportionexplained = results1_sds_varexplained/sum(results1_sds_varexplained)
    results1_vars = (results1_sds^2)[results1_sds_weights > component_varexplained_threshold]

    results2_ash <- ash(betahat_df[,2],
                        se_df[,2],
                        mixcompdist = "halfnormal",
                        prior = "uniform")

    results2_sds = unique(results2_ash$fitted_g$sd[-1])
    results2_sds_weights <- sapply(results2_sds, function(x) sum(results2_ash$fitted_g$pi[results2_ash$fitted_g$sd == x]))
    #results2_sds_varexplained = (results2_sds_weights * sum(complete.cases(betahat_df))) *(results2_sds^2)
    #results2_sds_proportionexplained = results2_sds_varexplained/sum(results2_sds_varexplained)
    results2_vars = (results2_sds^2)[results2_sds_weights > component_varexplained_threshold]

    U_adaptive_grid = list()
    i = 1
    corrs = seq(-1,1,length.out = 9)
    for (results1_var in c(0,results1_vars)) {
      for (results2_var in c(0,results2_vars)) {
        for (corr in corrs) {
          U_i = diag(2)
          U_i[1,1] = results1_var
          U_i[2,2] = results2_var
          U_i[1,2] = corr * sqrt(U_i[1,1] * U_i[2,2])
          U_i[2,1] = corr * sqrt(U_i[1,1] * U_i[2,2])
          U_adaptive_grid[[i]] = U_i
          i = i + 1
        }
      }
    }

    U = unique(U_adaptive_grid)
    prior = c(weight_nocorr,
              sapply(1:length(U),
                   function(x) {
                     covarmat = U[[x]]
                     if (covarmat[1,2] !=0) {return(1)} else {return(weight_nocorr)}
                   }))
    if (estimate_sampling_covariance) {
      message("Estimate sampling correlation")
      V.em = mash_estimate_corr_em(data, U, grid = 1,normalizeU = FALSE)
      data = mash_update_data(data, V=V.em$V)
    }

    m   = mash(data, U, grid = 1,normalizeU = FALSE,prior=prior)
  }

  message("Run mash")

  message("Compute mixture covariance matrix")
  covariance_matrices = array(data = NA, dim = c(2,2,length(m$fitted_g$pi)))
  covariance_matrices[,,1] <- 0

  for (weightnum in 1:length(m$fitted_g$grid)){
    for (ulistnum in 1:length(U)) {
      covariance_matrices[,,((weightnum - 1)*length(U)) + ulistnum + 1] <- m$fitted_g$Ulist[ulistnum][[1]] * m$fitted_g$grid[weightnum]
    }
  }


  weighted_covariance_matrices = array(data = NA, dim = c(2,2,length(m$fitted_g$pi)))

  for (matrixnum in 1:length(m$fitted_g$pi)) {
    weighted_covariance_matrices[,,matrixnum] <- covariance_matrices[,,matrixnum] * m$fitted_g$pi[matrixnum]
  }


  mixture_covariance_matrix <- rowSums(weighted_covariance_matrices,dims = 2)

  if (any(diag(mixture_covariance_matrix) == 0)) {
    mixture_correlation_matrix = matrix(NA, nrow = 2, ncol = 2)
    message("One Perturbation has 0 Effect Size Variance; Correlation not defined")
  } else {
    mixture_correlation_matrix = solve(diag(sqrt(diag(mixture_covariance_matrix)))) %*% mixture_covariance_matrix %*% solve(diag(sqrt(diag(mixture_covariance_matrix))))
  }

  #generate samples

  if (!is.null(n_sample)) {
    component = sample(1:length(m$fitted_g$pi), size = n_sample, prob = m$fitted_g$pi, replace = TRUE)

    component_table = table(component)

    posterior_samples = Reduce(rbind,
                               lapply(1:length(component_table),
                               function(i) {
                                 componentnum = as.numeric(names(component_table))[i]
                                 mvtnorm::rmvnorm(component_table[i],sigma = covariance_matrices[,,componentnum])
                               }))
  } else {
    posterior_samples = NA
  }

  #Return output
  end_time = Sys.time()
  timediff = difftime(end_time,start_time, units = "mins")
  return(list(TI_correlation = mixture_correlation_matrix[1,2],
              correlation_matrix = mixture_correlation_matrix,
              covariance_matrix = mixture_covariance_matrix,
              cor_raw = cor_raw,
              fitted = m$fitted_g,
              posterior_samples = posterior_samples,
              V = data$V,
              loglik = m$loglik,
              runtime = timediff))

}

