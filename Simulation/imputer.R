# mean imputer
imputer_mean = list(
  fit = function(dat_train) {
    colMeans(dat_train, na.rm = TRUE) # compute and store column means from training set
  },
  predict = function(model, dat_mask, i, j, m) {
    rep(model[j], m) # impute test entry using stored column mean
  }
)


# nearest neighbor hot deck
imputer_nn_hotdeck = list(
  fit = function(dat_train, k_nn = 10) {
    # k_nn: number of nearest neighbors (default to 10)
    list(train = dat_train, k_nn = k_nn)
  },
  predict = function(model, dat_mask, i, j, m) {
    Xtr = model$train
    
    # observed features in test row i (exclude j)
    feats  = which(!is.na(dat_mask[i, ]) & seq_len(ncol(dat_mask)) != j)
    
    # candidate donors
    donors = if (length(feats)) {
      which(!is.na(Xtr[, j]) &
              rowSums(is.na(Xtr[, feats, drop = FALSE])) == 0)
    } else {
      which(!is.na(Xtr[, j]))
    }
    if (!length(donors)) return(NA)
    if (!length(feats)) return(Xtr[sample(donors, m), j])
    
    # Euclidean distance
    xi = as.numeric(dat_mask[i, feats, drop = FALSE])
    Xd = as.matrix(Xtr[donors, feats, drop = FALSE])
    d2 = rowSums((Xd - matrix(xi, nrow = length(donors),
                              ncol = length(feats), byrow = TRUE))^2)
    
    # select nearest neighbors
    ord   = order(d2)
    k_eff = min(model$k_nn, length(donors)) 
    nn    = donors[ord[seq_len(k_eff)]]
    Xtr[sample(nn, m, replace=T), j]
  }
)

# random hot deck
imputer_random_hotdeck = list(
  fit = function(dat_train) {
    p = ncol(dat_train)
    donors_by_col = lapply(1:p, function(j) which(!is.na(dat_train[, j]))) # X_j observed 
    list(
      donors_by_col = donors_by_col,
      train_values  = dat_train
    )
  },
  predict = function(model, dat_mask, i, j, m=1) {
    donors = model$donors_by_col[[j]]
    if (!length(donors)) return(rep(NA, m))
    d = sample(donors, m, replace=T)
    model$train_values[d, j]
  }
)

# g-mmg
imputer_gmmg = list(
  fit = function(dat_train, graph) {
    out = g_mmg(data = dat_train, graph = graph, m = 1)
    out$graph = graph
    out
  },
  predict = function(model, dat_mask, i, j, m=1) {
    est.MLE = model$est.MLE
    orig_R = paste(as.integer(!is.na(dat_mask[i, ])), collapse = "")
    model_ids = model$pat %>% filter(orig_pattern == orig_R) %>% dplyr::select(model_pattern)
    if (nrow(model_ids) == 0) {
      d          = ncol(dat_mask)
      graph      = model$graph
      var_names  = colnames(dat_mask)
      missingpattern = data.frame(pattern = orig_R, n = 1)
      new_pat_list = model_patterns_list(
        missing_patterns = missingpattern,
        graph_structure  = coerce_graph(graph, dat_mask),
        d                = d,
        var_names        = var_names
      )
      model_id = apply(new_pat_list[[1]]$model_patterns, 1, paste0, collapse = "")
      model_id = model_id[substr(model_id, j, j) == "1"]
    } else {
      model_id = model_ids %>%
        filter(substr(model_pattern, j, j) == "1") %>%
        pull(model_pattern)
    }

    unique_id = model$unique_model_pat %>% rowwise() %>%
      mutate(pattern = paste(c_across(-unique_id), collapse = "")) %>% ungroup() %>%
      filter(pattern %in% model_id) %>% pull(unique_id)
    
    if(length(unique_id) == 0){
      pat.df = model$unique_model_pat %>% rowwise() %>%
        mutate(pattern = paste(c_across(-unique_id), collapse = "")) %>% ungroup()
      
      # when length(unique_id) == 0: only consider patterns where target j is observed
      candidates = pat.df[substr(pat.df$pattern, j, j) == "1", ]
      if (nrow(candidates) > 0) {
        overlap = sapply(candidates$pattern, function(p)
          sum(strsplit(p,"")[[1]] == strsplit(orig_R,"")[[1]] & strsplit(orig_R,"")[[1]] == "1"))
        unique_id = candidates$unique_id[which.max(overlap)]
      } else {
        # if none exist, fall back to complete case
        unique_id = which(pat.df$pattern == "11111111")
      }
    }
    
    mu_star = est.MLE[[unique_id]]$mu.hat  # MLE for mu
    Sigma_star = est.MLE[[unique_id]]$Sigma.hat  # MLE for Sigma
    
    observed_indices = which(as.numeric(strsplit(orig_R, "")[[1]]) == 1 & 
                               as.numeric(strsplit(model_id, "")[[1]]) == 1)  # Variables observed in both
    imp_indices = j
    
    observed_local = match(observed_indices, which(as.numeric(strsplit(model_id, "")[[1]]) == 1))
    imp_local = match(imp_indices, which(as.numeric(strsplit(model_id, "")[[1]]) == 1))
    
    # Impute
    x_obs = as.numeric(dat_mask[i, observed_indices])
    mu_obs = mu_star[observed_local]
    mu_imp = mu_star[imp_local]
    
    Sigma_oo = Sigma_star[observed_local, observed_local, drop = FALSE]
    Sigma_mo = Sigma_star[imp_local, observed_local, drop = FALSE]
    Sigma_mm = Sigma_star[imp_local, imp_local, drop = FALSE]
    
    if (length(observed_local) == 0) {
      # no observed vars: marginal
      cond_mean = mu_imp
      cond_var = Sigma_mm
    } else {
      Sigma_oo_inv = MASS::ginv(Sigma_oo)
      cond_mean = mu_imp + Sigma_mo %*% Sigma_oo_inv %*% (x_obs - mu_obs)
      cond_var = Sigma_mm - Sigma_mo %*% Sigma_oo_inv %*% t(Sigma_mo)
    }
    cond_var = (cond_var + t(cond_var)) / 2
    eigs = eigen(cond_var, symmetric = TRUE)$values
    if (any(eigs < 1e-8)) {
      cond_var = cond_var + diag(abs(min(eigs)) + 1e-6, nrow(cond_var))
    }
    
    # sample
    rnorm(n=m, mean=cond_mean, sd=sqrt(cond_var))
  }
)

# CCMV
imputer_ccmv_gaussian = list(
  fit = function(dat_train){
    X = as.data.frame(dat_train)
    num_cols = vapply(X, is.numeric, logical(1))
    Xnum = X[, num_cols, drop=F]
    CC = Xnum[complete.cases(Xnum), , drop=F]
    list(
      num_cols = num_cols,
      mu       = colMeans(CC),
      Sigma    = cov(CC)
    )
  },
  predict = function(model, dat_mask, i, j, m=1){
    num_cols = model$num_cols
    j_num = match(j, which(num_cols))
    mu    = model$mu
    Sigma = model$Sigma
    
    # observed numeric features in row i
    Xnum = dat_mask[, num_cols]
    obs_idx_num = which(!is.na(Xnum[i, ]) & seq_len(ncol(Xnum)) != j_num)
    mu_j = mu[j_num]
    
    if (length(obs_idx_num) == 0) {
      # fall back to marginal mean
      var_j = Sigma[j_num, j_num]
      if (is.na(var_j) || !is.finite(var_j) || var_j < .Machine$double.eps) return(mu_j)
      return(rnorm(m, mean = mu_j, sd = sqrt(var_j)))
    }
    
    mu_obs   = mu[obs_idx_num]
    S_j_obs  = Sigma[j_num, obs_idx_num, drop=F]
    S_obsobs = Sigma[obs_idx_num, obs_idx_num, drop=F]
    if (any(!is.finite(S_obsobs)) || 
        (is.matrix(S_obsobs) && qr(S_obsobs)$rank < ncol(S_obsobs)) ||
        (is.numeric(S_obsobs) && length(S_obsobs) == 1 && abs(S_obsobs) < .Machine$double.eps)) {
      return(rep(mu_j, m)) # fall back to the marginal mean if matrix is singular
    }
    x_obs = as.numeric(Xnum[i, obs_idx_num])
    S_obsobs_inv = solve(S_obsobs)
    cond_mean = as.numeric(mu_j + S_j_obs %*% S_obsobs_inv %*% (x_obs - mu_obs))
    cond_var = as.numeric(Sigma[j_num, j_num] - S_j_obs %*% S_obsobs_inv %*% t(S_j_obs))
    cond_var = max(cond_var, 1e-6) # guard tiny negatives from numerics
    if (is.na(cond_var) || !is.finite(cond_var) || cond_var < .Machine$double.eps) {
      return(rep(cond_mean, m))
    }
    return(rnorm(m, mean = cond_mean, sd = sqrt(cond_var)))
  }
)

# EM Gaussian
imputer_em_gaussian = list(
  fit = function(dat_train, tol = 1e-4, max_iter = 100, reg = 1e-6) {
    X = as.matrix(dat_train)
    fit = em_mvn(X, tol = tol, max_iter = max_iter, reg = reg)
    list(mu = fit$mu, Sigma = fit$Sigma)
  },
  predict = function(model, dat_mask, i, j, m = 1) {
    mu = model$mu
    S  = model$Sigma
    X  = as.matrix(dat_mask)
    D  = ncol(X)
    
    # observed indices in row i
    obs = which(!is.na(X[i, ]) & seq_len(D) != j)
    
    if (length(obs) == 0) {
      return(rnorm(m, mean = mu[j], sd = sqrt(S[j, j])))
    }
    
    # conditional Gaussian
    S_oo = S[obs, obs, drop = F]
    S_oj = S[obs, j, drop = F]
    beta = tryCatch(solve(S_oo, S_oj), error = function(e) rep(0, length(obs)))
    
    cond_mean = mu[j] + t(beta) %*% (X[i, obs] - mu[obs])
    cond_var  = S[j, j] - t(S_oj) %*% beta
    cond_var  = max(cond_var, 0)  # numerical safeguard
    
    rnorm(m, mean = cond_mean, sd = sqrt(cond_var))
  }
)

# MICE
imputer_mice = list(
  fit = function(dat_train) {
    list(train = dat_train)
  },
  predict = function(model, dat_mask, i, j, m=1){
    combined = rbind(model$train, dat_mask[i,])
    imp = mice::mice(combined, m=m, printFlag = FALSE)
    draws = lapply(1:m, function(k) {
      comp = mice::complete(imp, k)
      comp[nrow(combined), j]
    })
    return (unlist(draws))
  }
)

