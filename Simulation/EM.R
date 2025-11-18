########### EM ############
library(mvtnorm)

em_mvn = function(X, tol = 1e-4, max_iter = 200, reg = 1e-6) {
  X = as.matrix(X)
  N = nrow(X)
  D = ncol(X)
  
  # initialize
  mu = colMeans(X, na.rm = TRUE)
  complete_cases = na.omit(X)
  
  if (nrow(complete_cases) >= D) {
    Sigma = cov(complete_cases)
  } else {
    v = apply(X, 2, function(col) {
      v = var(col, na.rm = TRUE)
      ifelse(is.finite(v), v, 1)
    })
    Sigma = diag(v)
  }
  
  diag(Sigma) = diag(Sigma) + reg
  ll_old = -Inf
  
  # EM loop
  for (t in 1:max_iter) {
    # E-Step
    E_X_sum = numeric(D)
    E_XXT_sum = matrix(0, D, D)
    
    for (n in 1:N) {
      obs_idx = which(!is.na(X[n, ]))
      mis_idx = which(is.na(X[n, ]))
      
      if (length(mis_idx) == 0) {
        x_n = X[n, ]
        E_X_sum = E_X_sum + x_n
        E_XXT_sum = E_XXT_sum + tcrossprod(x_n)
        next
      }
      
      if (length(obs_idx) == 0) {
        E_X_sum = E_X_sum + mu
        E_XXT_sum = E_XXT_sum + tcrossprod(mu) + Sigma
        next
      }
      
      mu_o = mu[obs_idx]
      mu_m = mu[mis_idx]
      Sigma_oo = Sigma[obs_idx, obs_idx, drop = FALSE]
      Sigma_om = Sigma[obs_idx, mis_idx, drop = FALSE]
      Sigma_mm = Sigma[mis_idx, mis_idx, drop = FALSE]
      
      Sigma_oo_reg = Sigma_oo + diag(reg, length(obs_idx))
      A = solve(Sigma_oo_reg, Sigma_om)
      
      x_obs = X[n, obs_idx]
      cond_mu = mu_m + t(A) %*% (x_obs - mu_o)
      cond_Sigma = Sigma_mm - t(A) %*% Sigma_om
      
      E_X_n = numeric(D)
      E_X_n[obs_idx] = x_obs
      E_X_n[mis_idx] = cond_mu
      E_X_sum = E_X_sum + E_X_n
      
      E_XXT_n = tcrossprod(E_X_n)
      E_XXT_n[mis_idx, mis_idx] = E_XXT_n[mis_idx, mis_idx] + cond_Sigma
      E_XXT_sum = E_XXT_sum + E_XXT_n
    }
    
    # M-Step
    mu_new = E_X_sum / N
    Sigma_new = E_XXT_sum / N - tcrossprod(mu_new)
    
    diag(Sigma_new) = diag(Sigma_new) + reg
    Sigma_new = 0.5 * (Sigma_new + t(Sigma_new))
    
    # Calculate log-likelihood for convergence check
    ll_new = 0
    for (n in 1:N) {
      obs_idx = which(!is.na(X[n, ]))
      if (length(obs_idx) > 0) {
        ll_new = ll_new + mvtnorm::dmvnorm(x = X[n, obs_idx], mean = mu_new[obs_idx], 
                                           sigma = Sigma_new[obs_idx, obs_idx, drop = FALSE], 
                                           log = TRUE)
      }
    }
    
    # Convergence check
    if ((ll_new - ll_old) < tol) {
      break
    }
    
    # Update parameters for next iteration
    mu = mu_new
    Sigma = Sigma_new
    ll_old = ll_new
  }
  
  return(list(mu = mu, Sigma = Sigma))
}

