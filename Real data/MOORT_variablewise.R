# perform variable-wise MOORT

MOORT_varwise = function(dat, folds, crossfit = 3, M = 20, imputers = list(), graph = NULL) {
  
  n = nrow(dat)
  d = ncol(dat)
  folds = folds
  results = lapply(names(imputers), function(name) {
    list(S_by_var = vector("list", d))
  })
  names(results) = names(imputers)
  
  cat("Starting cross-fitting with", crossfit, "folds and", length(imputers), "imputers.\n")
  cat("----------------------------------------------------------\n")
  
  for (fold_k in 1:crossfit) {
    test_idx  = which(folds == fold_k)
    train_idx = which(folds != fold_k)
    
    dat_train = dat[train_idx, , drop=FALSE]
    dat_test  = dat[test_idx, , drop=FALSE]
    
    # train
    models = list()
    for (name in names(imputers)) {
      imp = imputers[[name]]
      if (name == "mmg") {
        models[[name]] = imp$fit(dat_train, graph = graph)
      } else {
        models[[name]] = imp$fit(dat_train)
      }
    }
    
    cat(sprintf("Finished training all imputers for fold %d.\n", fold_k))
    cat(" Starting variable-wise loop.\n")
    
    # --- variable-wise loop ---
    for (j in 1:d) {
      Dj = which(!is.na(dat_test[, j]))   # indices where X_j is observed
      for (i in Dj) {
        x_true = dat_test[i, j]
        dat_mask = dat_test
        dat_mask[i, j] = NA
        
        for (name in names(imputers)) {
          imp = imputers[[name]]
          imp_vals = imp$predict(models[[name]], dat_mask, i, j, m = M)
          
          G_hat = ecdf(imp_vals)
          S_val = G_hat(x_true)
          
          if (is.null(results[[name]]$S_by_var[[j]])) {
            results[[name]]$S_by_var[[j]] = numeric(0)
          }
          results[[name]]$S_by_var[[j]] = c(results[[name]]$S_by_var[[j]], S_val)
        }
      }
    }
    
    cat(sprintf("\n Completed fold %d / %d.\n", fold_k, crossfit))
  }
  
  out_tbl = data.frame(method = character(), var = integer(), ks_stat = numeric())
  for (name in names(results)) {
    for (j in 1:d) {
      S_vals = results[[name]]$S_by_var[[j]]
      ks_stat = ks.test(S_vals, "punif", 0, 1)$statistic
      out_tbl = rbind(out_tbl, data.frame(method = name, var = j, ks_stat = ks_stat))
    }
  }
  
  return(list(results = results, per_variable = out_tbl))
}