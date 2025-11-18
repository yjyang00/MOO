# mean imputer
imputer_mean = list(
  fit = function(dat_train) {
    colMeans(dat_train, na.rm = TRUE) # compute and store column means from training set
  },
  predict = function(model, dat_mask, i, j, m=1) {
    rep(model[j], m) # impute test entry using stored column mean
  }
)

# NN hot deck
imputer_nn_hotdeck = list(
  fit = function(dat_train, k_nn = 10) {
    list(train = dat_train, k_nn = k_nn)
  },
  predict = function(model, dat_mask, i, j, m) {
    Xtr   = model$train
    k_nn  = model$k_nn
    
    # observed features in test row i (excluding j)
    feats = which(!is.na(dat_mask[i, ]) & seq_len(ncol(dat_mask)) != j)
    
    if (length(feats) > 0) {
      donors = which(!is.na(Xtr[, j]) & rowSums(is.na(Xtr[, feats, drop = FALSE])) == 0)
    } else {
      donors = which(!is.na(Xtr[, j]))
    }
    
    # if no donors found, fall back to marginal distribution of observed Xj
    if (length(donors) == 0) {
      observed_vals = Xtr[[j]][!is.na(Xtr[[j]])]
      if (length(observed_vals) == 0) return(NA)
      return(sample(observed_vals, m, replace = TRUE))
    }
    
    # if no features, sample from marginal donors
    if (length(feats) == 0) {
      return(sample(Xtr[donors, j, drop = TRUE], m, replace = TRUE))
    }
    
    # Euclidean distance
    xi = as.numeric(dat_mask[i, feats, drop = FALSE])
    Xd = as.matrix(Xtr[donors, feats, drop = FALSE])
    d2 = rowSums((Xd - matrix(xi, nrow = length(donors),
                              ncol = length(feats), byrow = TRUE))^2)
    
    # select nearest neighbors
    ord   = order(d2)
    k_eff = min(k_nn, length(donors))
    nn    = donors[ord[seq_len(k_eff)]]
    
    # sample from nearest neighbors
    val = Xtr[nn, j, drop = TRUE]
    return(sample(val, m, replace = TRUE))
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
    return(as.data.frame(model$train_values)[d, j, drop=T])
  }
)

# mmg
imputer_mmg = list(
  fit = function(dat_train, graph) {
    X = data.frame(intercept = rep(1, nrow(dat_train)))
    Y = as.data.frame(lapply(dat_train, as.integer))
    out = mbpe_mmg(X=X, Y=Y, K=5, em_num_init=10, em_max_iter = 300, em_tol = 1e-04, 
                   graph = graph, m = 1, s_max = rep(max(dat_train[1], na.rm=T), 5))
    out$graph = graph
    out
  },
  predict = function(model, dat_mask, i, j, m=1) {
    X = data.frame(intercept = rep(1, nrow(dat_mask)))
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
      
      candidates = pat.df[substr(pat.df$pattern, j, j) == "1", ]
      if (nrow(candidates) > 0) {
        overlap = sapply(candidates$pattern, function(p)
          sum(strsplit(p,"")[[1]] == strsplit(orig_R,"")[[1]] & strsplit(orig_R,"")[[1]] == "1"))
        unique_id = candidates$unique_id[which.max(overlap)]
      }
    }
    
    beta_star = est.MLE[[unique_id]]$beta_star  
    theta_star = est.MLE[[unique_id]]$theta_star 
    
    observed_indices = which(as.numeric(strsplit(orig_R, "")[[1]]) == 1 & 
                               as.numeric(strsplit(model_id, "")[[1]]) == 1)  # Variables observed in both
    
    imp_indices = j
    
    observed_local = match(observed_indices, which(as.numeric(strsplit(model_id, "")[[1]]) == 1))
    imp_local = match(imp_indices, which(as.numeric(strsplit(model_id, "")[[1]]) == 1))
    
    # Impute
    x_obs = as.numeric(dat_mask[i, observed_indices])
    logW_im = predict_probs(beta_star, X[i, , drop=FALSE], log=TRUE)
    
    if (length(observed_indices) > 0) {
      model_vec  = as.numeric(strsplit(model_id, "")[[1]])
      local_vars = which(model_vec == 1)
      for (k_iter in 1:nrow(theta_star)) {
        for (t in observed_indices) {
          theta_col = match(t, local_vars)
          logW_im[, k_iter] = logW_im[, k_iter] +
            dbinom(as.numeric(dat_mask[i, t]),
                   size = rep(max(dat_train[1], na.rm=T), 5)[t],
                   prob = theta_star[k_iter, theta_col],
                   log = TRUE)
        }
      }
    }
    
    W_im = exp(logW_im - apply(logW_im, 1, logsumexp))
    Z_star = sample(1:nrow(theta_star), size = 1, prob = W_im)
    
    model_vec = as.numeric(strsplit(model_id, "")[[1]])
    theta_col = match(j, which(model_vec == 1))
    Z_star = sample(1:nrow(theta_star), size = m, replace = TRUE, prob = W_im)
    probs  = theta_star[Z_star, theta_col]
    draws  = rbinom(m, size = rep(max(dat_train[1], na.rm=T), 5)[j], prob = probs)
    return(draws)
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
    
    # observed numeric features in row i (excluding j)
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

# EM-MBPE
imputer_em_mbpe = list(
  fit = function(dat_train, tol = 1e-4, max_iter = 300, reg = 1e-4) {
    X = data.frame(intercept = rep(1, nrow(dat_train)))
    Y = as.data.frame(lapply(dat_train, as.integer))
    # fit_mar_randinit is from from the mixturebpe package (more details: https://github.com/danielsuen/mixturebpe)
    fit = fit_mar_randinit(X, Y, K=5, M_mar=max_iter, numRandInit=3, s_max=rep(max(dat_train[1], na.rm=T), 5), 
                           stop_eps_mar=1e-4, stop_eps_latent=1e-4, num_imps=1) 
    list(beta = fit$beta_star, theta = fit$theta_star)  # store parameters
  },
  predict = function(model, dat_mask, i, j, m = 1) {
    beta  = as.matrix(model$beta)
    theta = as.matrix(model$theta) 
    Y = as.matrix(dat_mask)
    s_max = rep(max(dat_train[1], na.rm=T), 5)
    D  = ncol(Y)
    
    # observed indices in row i (excluding j)
    obs = which(!is.na(Y[i, ]) & seq_len(D) != j)
    X_i = matrix(1, 1, 1)
    logW_im = predict_probs(beta, X_i, log = TRUE) 
    
    if (length(obs) > 0) {
      for (t in obs) {
        for (k_iter in 1:5) {
          logW_im[1, k_iter] = logW_im[1, k_iter] +
            dbinom(as.numeric(Y[i, t]),
                   size = s_max[t],
                   prob = theta[k_iter, t],
                   log = TRUE)
        }
      }
    }
    
    logW_im = logW_im - max(logW_im)
    W_im = exp(logW_im)
    W_im = W_im / sum(W_im)
    
    Z_star = sample(1:nrow(theta), size = m, replace = TRUE, prob = W_im)
    probs  = theta[Z_star, j]
    draws  = rbinom(m, size = s_max[j], prob = probs)
    return(draws)
    
  }
)

imputer_mice = list(
  fit = function(dat_train) {
    list(train = dat_train)
  },
  predict = function(model, dat_mask, i, j, m){
    combined = rbind(model$train, dat_mask[i,])
    imp = mice::mice(combined, m=m, printFlag = FALSE)
    draws = lapply(1:m, function(k) {
      comp = mice::complete(imp, k)
      comp[nrow(combined), j]
    })
    return (unlist(draws))
  }
)


