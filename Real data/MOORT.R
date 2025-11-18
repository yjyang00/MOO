# Implement MOORT on real data

MOOQT_crossfit = function(dat, folds, crossfit = 5, M = 20, imputers = list(), graph = NULL) {
  
  n = nrow(dat)
  dat_NA = dat
  folds = folds 
  
  # Store
  results = vector("list", length(imputers))
  names(results) = names(imputers)
  for (name in names(imputers)) {
    results[[name]] = list(S = numeric(0))
  }
  
  cat("Starting cross-fitting with", crossfit, "folds and", length(imputers), "imputers.\n")
  cat("----------------------------------------------------------\n")
  
  for (fold_k in 1:crossfit) {
    test_idx  = which(folds == fold_k)
    train_idx = which(folds != fold_k)
    
    dat_train = dat_NA[train_idx, , drop=F]
    dat_test  = dat_NA[test_idx, , drop=F]
    
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
    
    for (i in 1:nrow(dat_test)) {
      obs_j = as.numeric(which(!is.na(dat_test[i, ])))
      if (!length(obs_j)) next
      if (length(obs_j) == 1) {
        j = obs_j
      } else {
        j = sample(obs_j, 1)
      }
      
      x_true = dat_test[i, j]
      
      # mask entry
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
    cat(sprintf("\n Completed fold %d / %d.\n", fold_k, crossfit))
  }
  
  methods_vec = character()
  ks_stats_vec = numeric()
  out_tbl = data.frame(method = character(), ks_stat = numeric(), 
                       fail_rate = numeric(), mark = character())
  
  for (name in names(results)) {
    S_vals = results[[name]]$S
    ks_test = ks.test(S_vals, "punif", 0, 1)
    out_tbl = rbind(out_tbl, data.frame(method = name, ks_stat = ks_test$statistic))
  }
  
  return(list(results, out_tbl))
}