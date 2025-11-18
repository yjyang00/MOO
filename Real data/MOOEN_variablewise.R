# variable-wise MOOEN
MOOEN_varwise = function(dat, folds, crossfit, M, imputers, graph) {
  
  n = nrow(dat)
  d = ncol(dat)
  folds = folds
  
  # Store results per method
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
        x_true = dat_test[[i, j]]
        dat_mask = dat_test
        dat_mask[i, j] = NA
        
        for (name in names(imputers)) {
          imp = imputers[[name]]
          imp1 = imp$predict(models[[name]], dat_mask, i, j, m = M) 
          imp2 = imp$predict(models[[name]], dat_mask, i, j, m = M)
          
          term1 = mean(abs(x_true - imp1))
          D = abs(outer(imp1, imp2, "-"))
          term2 = sum(D) / (2 * M * (M - 1))
          E_val = term1 - term2
          results[[name]]$S_by_var[[j]] = c(results[[name]]$S_by_var[[j]], E_val)
        }
      }
    }
    
    cat(sprintf("\n Completed fold %d / %d.\n", fold_k, crossfit))
  }
  
  out_tbl = data.frame(method = character(), var = integer(), E_score = numeric())
  for (name in names(results)) {
    for (j in 1:d) {
      S_vals = results[[name]]$S_by_var[[j]]
      out_tbl = rbind(out_tbl, data.frame(method = name, var = j, E_score = sum(S_vals)/n))
    }
  }
  # aggregate
  summary_tbl = aggregate(E_score ~ method, data = out_tbl, sum)
  colnames(summary_tbl)[2] = "total_risk"
  
  return(list(results = results, per_variable = out_tbl, summary = summary_tbl))
}

