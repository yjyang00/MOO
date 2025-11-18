#' Simulation code: MOO imputation risk computation
#'
#' Computes the mask-one-out (MOO) imputation risk for a collection of
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
#' @param m
#'   Number of multiple imputation.
#'
#' @return
#'   MOO imputation risk for each imputer in the imputers list.
#'
#'


MOO_crossfit = function(dat, rho, K, imputers = list(), graph,  m) {
  
  # MCAR
  dat_NA = mmg::genMCAR(dat, rho = rho)
  n = nrow(dat)
  
  folds = sample(rep(1:K, length.out = n))
  total_loss = setNames(as.list(rep(0, length(imputers))), names(imputers))
  
  for (k in 1:K) {
    test_idx  = which(folds == k)
    train_idx = which(folds != k)
    
    dat_train = dat_NA[train_idx, ]
    dat_test  = dat_NA[test_idx, ]
    
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
          x_hat = imputers[[name]]$predict(models[[name]], dat_mask, i, j, m = m)
          row_loss[[name]] = row_loss[[name]] + mean( (x_true - x_hat)^2 )
        }
      }
      
      for (name in names(imputers)) {
        total_loss[[name]] = total_loss[[name]] + row_loss[[name]]
      }
    }
  }
  
  # average
  avg_loss = setNames(as.numeric(unlist(total_loss)) / n, names(imputers))
  return(avg_loss)
}


#### Example
source("gmmg.R")
source("EM.R")
source("imputer.R")

dat = read.table("yacht.data")  # (available from UCI ML repo)
dat = dat[, -1] # remove the first column
colnames(dat) = paste0("x", 1:6)
yachto_model = EBICglasso(cov(dat), n = nrow(dat), threshold = T)
adj_mat = abs(yachto_model) != 0
diag(adj_mat) = 0
colnames(adj_mat) = rownames(adj_mat) = colnames(dat)
qgraph(adj_mat, labels = colnames(dat))
g_yacht = igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
dat_scale = scale(dat)
colnames(dat_scale) = colnames(dat)

# specify imputers
imputers = list(
  mean       = imputer_mean,
  em_mvn     = imputer_em_gaussian,
  nn_hotdeck = imputer_nn_hotdeck,
  g_mmg      = imputer_gmmg,
  CCMV       = imputer_ccmv_gaussian,
  mice       = imputer_mice
)

# one run, use for loop for multiple runs.
result = MOO_crossfit(dat_scale, rho = 0.3, K = 5, graph = g_yacht, 
                      imputers = imputers, m = 20) 


