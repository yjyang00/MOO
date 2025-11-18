# combined

imp_risk = function(dat, rho, K, M, imputers, graph=NULL){
  
  dat_NA = mmg::genMCAR(dat, rho = rho)
  n = nrow(dat)
  folds = sample(rep(1:K, length.out = n))
  
  # Initialize 
  total_loss_MOO = setNames(as.list(rep(0, length(imputers))), names(imputers))
  results_MOORT = lapply(names(imputers), function(name) list(S = numeric(0)))
  names(results_MOORT) = names(imputers)
  results_MOOEN = lapply(names(imputers), function(name) list(E = numeric(0)))
  names(results_MOOEN) = names(imputers)
  
  for (k in 1:K) {
    test_idx  = which(folds == k)
    train_idx = which(folds != k)
    
    dat_train = dat_NA[train_idx, , drop=F]
    dat_test  = dat_NA[test_idx, , drop=F]
    
    # train
    models = setNames(vector("list", length(imputers)), names(imputers))
    for (name in names(imputers)) {
      imp = imputers[[name]]
      if (name == "g_mmg"){
        models[[name]] = imp$fit(dat_train, graph=graph)
      }else { 
        models[[name]] = imp$fit(dat_train)
      }
    }
    
    # ----------- compute MOO, MOORT, MOOEN -------------
    ## MOO
    for (i in 1:nrow(dat_test)) {
      row_loss = setNames(as.list(rep(0, length(imputers))), names(imputers))
      obs_idx = which(!is.na(dat_test[i, ]))
      if (!length(obs_idx)) next
      
      for (j in obs_idx) {
        x_true = dat_test[i, j]
        
        # mask
        dat_mask = dat_test
        dat_mask[i, j] = NA
        
        for (name in names(imputers)) {
          x_hat = imputers[[name]]$predict(models[[name]], dat_mask, i, j, m = M)
          row_loss[[name]] = row_loss[[name]] + mean( (x_true - x_hat)^2 )
        }
      }
      
      for (name in names(imputers)) {
        total_loss_MOO[[name]] = total_loss_MOO[[name]] + row_loss[[name]]
      }
    }
    
    ## MOORT
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
        results_MOORT[[name]]$S = c(results_MOORT[[name]]$S, S_val)
      }
    }
    
    ## MOOEN
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
          results_MOOEN[[name]]$E = c(results_MOOEN[[name]]$E, E_val)
        }
      }
    }
    
  }
  
  avg_MOO = setNames(as.numeric(unlist(total_loss_MOO)) / n, names(imputers))
  avg_MOOEN = sapply(results_MOOEN, function(x) sum(x$E) / n)
  avg_MOORT = sapply(names(results_MOORT), function(name) {
    S_vals = results_MOORT[[name]]$S
    ks.test(S_vals, "punif", 0, 1)$statistic
  })
  
  # Combine into a single output table
  out_tbl = data.frame(
    Method = names(imputers),
    MOO = avg_MOO,
    MOORT = avg_MOORT,
    MOOEN = avg_MOOEN
  )
  
  return(out_tbl)
}


#### example ####
library(tidyverse)
dat = iris[,1:4]
iris_model = qgraph::EBICglasso(cov(dat), n = nrow(dat), threshold = T)
adj_mat = abs(iris_model) != 0
diag(adj_mat) = 0
colnames(adj_mat) = rownames(adj_mat) = colnames(dat)
qgraph::qgraph(adj_mat, labels = colnames(dat))
g_iris = igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
dat_scale = scale(dat)
imputers = list(
  mean       = imputer_mean,
  em_mvn     = imputer_em_gaussian,
  nn_hotdeck = imputer_nn_hotdeck,
  g_mmg      = imputer_gmmg,
  CCMV       = imputer_ccmv_gaussian
)
## note on efficiency: the implementation in simulation can be time-expensive depending on the 
## dimension of data and number of iterations. We thus recommend parallel computing.
res = imp_risk(dat = dat_scale, rho = 0.3, K = 5, M = 20, imputers = imputers, graph = g_iris) # one-run, use for loop for multiple runs.

# visualize PI diagram
ggplot(res) +
  geom_point(aes(x = MOO, y = MOOEN, color = Method), size = 3) +
  labs(title = "PI Diagram", x = "MOO", y = "MOOEN") +
  theme_minimal()

ggplot(res) +
  geom_point(aes(x = MOO, y = MOORT, color = Method), size = 3) +
  labs(title = "PI Diagram", x = "MOO", y = "MOORT") +
  theme_minimal()
