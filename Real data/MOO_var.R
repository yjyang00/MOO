# Perform variable-wise MOO
MOO_varwise = function(dat, folds, crossfit, M, imputers = list(), graph = NULL) {
 
  n = nrow(dat)
  d = ncol(dat)
  folds = folds
  
  # Store results per method
  total_loss = setNames(vector("list", length(imputers)), names(imputers))
  for (name in names(imputers)) {
    total_loss[[name]] = numeric(d) 
  }
  
  cat("Starting cross-fitting with", crossfit, "folds and", length(imputers), "imputers.\n")
  cat("----------------------------------------------------------\n")
  
  for (fold_k in 1:crossfit) {
    test_idx  = which(folds == fold_k)
    train_idx = which(folds != fold_k)
    
    dat_train = dat[train_idx, , drop=FALSE]
    dat_test  = dat[test_idx, , drop=FALSE]
    
    # train
    models = setNames(vector("list", length(imputers)), names(imputers))
    for (name in names(imputers)) {
      imp = imputers[[name]]
      if (name == "mmg") {
        models[[name]] = imp$fit(dat_train, graph = graph)
      } else {
        models[[name]] = imp$fit(dat_train)
      }
    }
    
    cat(sprintf("Finished training all imputers for fold %d.\n", fold_k))
    cat(" Starting variable-wise imputation and loss computation...\n")
    
    for (j in 1:d) {                           # fix variable j
      Dj = which(!is.na(dat_test[, j]))        # rows where X_ij is observed
      
      for (i in Dj) {                          # loop over D_j
        x_true = dat_test[[i, j]]
        
        dat_mask = dat_test
        dat_mask[i, j] = NA
        
        for (name in names(imputers)) {
          x_hat =  imputers[[name]]$predict(models[[name]], dat_mask, i, j, m = M)
          total_loss[[name]][j] = total_loss[[name]][j] + mean( (x_true - x_hat)^2 )

        }
      }
    }
    
    cat(sprintf("\n Completed fold %d / %d.\n", fold_k, crossfit))
  }
  
  avg_loss = lapply(total_loss, function(vec) vec / n)
  
  # tidy result
  results = do.call(rbind, lapply(names(avg_loss), function(name) {
    data.frame(
      method   = name,
      variable = seq_len(d),
      loss     = avg_loss[[name]]
    )
  }))
  
  return(results)
}
