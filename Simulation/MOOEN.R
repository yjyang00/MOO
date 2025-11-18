#' Simulation code: MOOEN imputation risk computation
#' 
#' Computed the mask-one-out with energy distance (MOOEN) for a collection of 
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
#'   MOOEN imputation risk for each imputer in the imputers list.
#'
#'


MOOEN_crossfit = function(dat, rho, K, M, imputers = list(), graph = NULL) {
  
  n = nrow(dat)
  dat_NA = mmg::genMCAR(dat, rho = rho)
  folds = sample(rep(1:K, length.out = n))
  results = lapply(names(imputers), function(name) list(E = numeric(0)))
  names(results) = names(imputers)
  
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
      obs_idx = which(!is.na(dat_test[i, ]))
      if (length(obs_idx) == 0) next
      
      for(j in obs_idx) {
        x_true = dat_test[i, j]
        
        # mask
        dat_mask = dat_test
        dat_mask[i, j] = NA
        
        for (name in names(imputers)) {
          imp = imputers[[name]]
          imp1 = imp$predict(models[[name]], dat_mask, i, j, m = M)
          imp2 = imp$predict(models[[name]], dat_mask, i, j, m = M)
          
          # energy distance
          term1 = mean(abs(x_true - imp1))
          D = abs(outer(imp1, imp2, "-"))
          term2 = sum(D) / (2 * M * (M - 1))
          E_val = term1 - term2
          results[[name]]$E = c(results[[name]]$E, E_val)
        }
      }
    }
  }
  
  out_tbl = data.frame(method = names(results), MOOEN = sapply(results, function(x) sum(x$E)/n))
  
  return(out_tbl)
}
