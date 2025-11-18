#' Simulation code: MOORT imputation risk computation
#' 
#' Computed the mask-one-out rank transformation (MOORT) for a collection of 
#' imputation models.
#'
#' @param dat
#'   A complete data matrix or data.frame. This is the fully observed dataset
#'   before missingness is generated.
#'
#' @param rho
#'   MCAR missingness rate. Each entry is independently set to \code{NA} 
#'   with probability \code{rho}.
#'
#' @param K
#'   Number of folds for cross-fitting. 
#'
#' @param imputers
#'   A list of imputation methods.
#'   A typical list looks like:
#'   \preformatted{
#'   imputers = list(
#'     mean            = imputer_mean,
#'     em_mvn          = imputer_em_gaussian,
#'     nn_hotdeck      = imputer_nn_hotdeck,
#'     g_mmg           = imputer_gmmg,
#'     CCMV            = imputer_ccmv_gaussian,
#'     mice            = imputer_mice
#'   )
#'   }
#' @param graph
#'   A graph structure required by graph-based imputers (e.g., Markov
#'   missing graph). 
#'
#' @param M
#'   Number of multiple imputation.
#'
#' @return
#'   MOORT imputation risk for each imputer in the imputers list.
#'
#'

MOORT_crossfit = function(dat, rho, K, M, imputers, graph) {
  
  n = nrow(dat)
  dat_NA = mmg::genMCAR(dat, rho = rho)
  folds = sample(rep(1:K, length.out = n))
  results = vector("list", length(imputers))
  names(results) = names(imputers)
  
  for (name in names(imputers)) {
    results[[name]] = list(S = numeric(0))
  }
  
  # Cross-fitting
  for (k in 1:K) {
    test_idx  = which(folds == k)
    train_idx = which(folds != k)
    
    dat_train = dat_NA[train_idx, , drop=F]
    dat_test  = dat_NA[test_idx, , drop=F]
    
    # train
    models = list()
    for (name in names(imputers)) {
      imp = imputers[[name]]
      if (name == "g_mmg") {
        models[[name]] = imp$fit(dat_train, graph = graph)
      } else {
        models[[name]] = imp$fit(dat_train)
      }
    }
    
    for (i in 1:nrow(dat_test)) {
      obs_j = as.numeric(which(!is.na(dat_test[i, ])))
      if (!length(obs_j)) next
      if (length(obs_j) == 1) {
        j = obs_j
      } else {
        j = sample(obs_j, 1)
      }
      
      x_true = dat_test[i, j]
      
      # mask
      dat_mask = dat_test
      dat_mask[i, j] = NA
      
      for (name in names(imputers)) {
        imp = imputers[[name]]
        imp_vals = imp$predict(models[[name]], dat_mask, i, j, m = M)
        G_hat = ecdf(imp_vals)
        S_val = G_hat(x_true)
        results[[name]]$S = c(results[[name]]$S, S_val)
      }
    }
  }
  
  
  out_tbl = data.frame(method = character(), ks_stat = numeric())
  
  for (name in names(results)) {
    S_vals = results[[name]]$S
    ks_test = ks.test(S_vals, "punif", 0, 1)
    out_tbl = rbind(out_tbl, data.frame(method = name, ks_stat = ks_test$statistic))
  }
  
  return(out_tbl)
}
