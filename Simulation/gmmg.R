########### GMMG ############
coerce_graph = function(graph, data, tol = 1e-12) {
  stopifnot(is.data.frame(data) || is.matrix(data))
  vars = colnames(data)
  if (is.null(vars)) vars = paste0("x", seq_len(ncol(data)))
  
  # helper functions
  as_adjlist = function(lst) {
    if (is.null(names(lst))) {
      if (length(lst) != length(vars)) stop("Graph list must match number of variables.")
      names(lst) = vars
    }
    if (!all(names(lst) %in% vars)) stop("Graph list has names not in data.")
    lst = lapply(lst, \(v) unique(as.character(v)))
    for (v in setdiff(vars, names(lst))) lst[[v]] = character()
    for (u in vars) {
      for (v in lst[[u]]) {
        if (!(u %in% lst[[v]])) lst[[v]] = c(lst[[v]], u)
      }
    }
    lapply(lst, sort)[vars]
  }
  
  
  if (is.list(graph) && !inherits(graph, "igraph") && !is.data.frame(graph)) {
    return(as_adjlist(graph))
  }
  
  # matrix
  if (is.matrix(graph)) {
    if (nrow(graph) != length(vars) || ncol(graph) != length(vars))
      stop("Matrix graph must have dims equal to ncol(data).")
    if (!is.null(rownames(graph)) && !all(rownames(graph) == vars)) {
      graph = graph[vars, vars, drop = FALSE]
    }
    is_sym = isTRUE(all(abs(graph - t(graph)) < tol))
    if (is_sym) {
      adj = (abs(graph) > tol) * 1
      diag(adj) = 0
      lst = lapply(seq_len(ncol(adj)), function(j) vars[which(adj[j, ] != 0)])
      names(lst) = vars
      return(as_adjlist(lst))
    }
  }
  
  # edge list
  if (is.data.frame(graph) || (is.matrix(graph) && ncol(graph) == 2)) {
    ed = as.data.frame(graph, stringsAsFactors = FALSE)
    if (ncol(ed) != 2) stop("Edge list must have exactly 2 columns.")
    colnames(ed) = c("from", "to")
    ed$from = as.character(ed$from); ed$to = as.character(ed$to)
    if (!all(ed$from %in% vars) || !all(ed$to %in% vars))
      stop("Edge list contains names not found in data.")
    lst = setNames(vector("list", length(vars)), vars)
    for (v in vars) lst[[v]] = character()
    for (k in seq_len(nrow(ed))) {
      u = ed$from[k]; v = ed$to[k]
      if (u == v) next
      lst[[u]] = c(lst[[u]], v)
      lst[[v]] = c(lst[[v]], u)
    }
    return(as_adjlist(lst))
  }
  
  # igraph
  if (inherits(graph, "igraph")) {
    if (!requireNamespace("igraph", quietly = TRUE))
      stop("Please install 'igraph' to pass an igraph object.")
    g = graph
    g_names = igraph::V(g)$name
    if (is.null(g_names)) stop("igraph object must have vertex names.")
    missing = setdiff(vars, g_names)
    if (length(missing)) {
      g = igraph::add_vertices(g, nv = length(missing), name = missing)
    }
    g = igraph::induced_subgraph(g, vids = vars)
    lst = lapply(vars, function(v) sort(igraph::neighbors(g, v, mode = "all")$name))
    names(lst) = vars
    return(as_adjlist(lst))
  }
  
  stop("Unsupported 'graph' input.")
}

model_patterns_list = function(missing_patterns, graph_structure, d, var_names) {
  if (is.null(var_names)) {
    var_names = names(graph_structure)
  }
  pattern_matrix = as.data.frame(do.call(rbind, strsplit(as.character(missing_patterns[[1]]), "")))
  pattern_matrix = as.data.frame(lapply(pattern_matrix, as.numeric))
  colnames(pattern_matrix) = var_names
  
  model_patterns = vector("list", nrow(missing_patterns))
  
  model_patterns = purrr::map(1:nrow(missing_patterns), function(i){
    r = pattern_matrix[i, ]  # current missing pattern
    missing_vars = names(r)[r == 0]  # identify missing variables
    observed_vars = names(r)[r == 1]  # identify observed variables
    
    # identify missing variable blocks
    block_patterns = list()
    while (length(missing_vars) > 0) {
      queue = missing_vars[1]
      block = character(0)
      
      while (length(queue) > 0) {
        current = queue[1]
        queue = queue[-1]
        
        if (!(current %in% block)) {
          block = c(block, current)
          
          # find missing neighbors and add them to queue
          missing_neighbors = intersect(graph_structure[[current]], missing_vars)
          queue = c(queue, missing_neighbors)
        }
      }
      
      block_patterns[[length(block_patterns) + 1]] = block
      missing_vars = setdiff(missing_vars, block)  # remove processed variables
    }
    
    body_patterns = purrr::map(block_patterns, function(block) {
      observed_neighbors = unique(unlist(purrr::map(block, ~ intersect(graph_structure[[.x]], observed_vars))))
      
      body_pattern = rep(0, length(r))
      names(body_pattern) = names(r)
      
      if (length(observed_neighbors) > 0) {
        body_pattern[c(block, observed_neighbors)] = 1
      } else {
        body_pattern[block] = 1  # no neighbors observed: marginal
      }
      
      return(body_pattern)
    })
    
    # store model patterns
    model_patterns = if (length(body_patterns) == 0) {
      matrix(rep(1, length(r)), nrow = 1)  # If no decomposition, full pattern
    } else {
      unique(do.call(rbind, body_patterns))  # store unique body patterns
    }
    
    # store results in a structured format
    list(
      pattern = missing_patterns$pattern[i],
      count = missing_patterns$n[i],
      model_patterns = model_patterns
    )
  })
  return(model_patterns)
}


# process model patterns and returns a data frame
model_patterns_df = function(model_patterns, d, var_names) {
  model_patterns_df_list = list()
  
  for (i in 1:length(model_patterns)) {
    orig_pattern = model_patterns[[i]]$pattern
    model_df = as.data.frame(as.matrix(model_patterns[[i]]$model_patterns))
    colnames(model_df) = var_names
    model_df = cbind(orig_pattern = orig_pattern, model_df)
    model_patterns_df_list[[i]] = model_df
  }
  
  final_model_patterns_df = do.call(rbind, model_patterns_df_list)
  return(final_model_patterns_df)
}

# extract fully observed subsets of data for each model pattern
extract_observed_units = function(dat, d, unique_model_patterns, var_names) {
  extracted_Y = list()
  for (i in 1:nrow(unique_model_patterns)) {
    observed_mask = as.logical( unique_model_patterns[i, var_names] ) # Y_observed
    complete_rows = complete.cases(dat[, observed_mask, drop = F]) # Select rows whose columns in obs_mask are observed
    extracted_Y[[i]] = dat[complete_rows, observed_mask, drop = F]
  }
  return(extracted_Y)
}


g_imp = function(dat, d, est.MLE, orig_pattern, patterns_join, unique_model_patterns){
  X_imputed = dat
  
  for(i in 1:nrow(dat)){
    orig_R = orig_pattern[i]
    
    if(orig_R == paste(rep(1, d), collapse = "")){
      X_imputed[i, ] = dat[i, ]
      next # complete case: jump to next iteration
    }
    
    model_ids = patterns_join %>% filter(orig_pattern == orig_R) %>% pull(unique_id) # model pattern IDs for the current row
    
    if (length(model_ids) == 0) {
      warning(paste("No model patterns found at row", i))
      next
    }
    
    for(model_id in model_ids){
      model_R = unique_model_patterns[model_id, 1:d]  # Get the model pattern
      mu_star = est.MLE[[model_id]]$mu.hat  # MLE for mu
      Sigma_star = est.MLE[[model_id]]$Sigma.hat  # MLE for Sigma
      
      observed_indices = which(as.numeric(strsplit(orig_R, "")[[1]]) == 1 & model_R == 1)  # variables observed in both
      imp_indices = which(model_R == 1 & as.numeric(strsplit(orig_R, "")[[1]]) == 0)  # variables to impute
      
      observed_local = match(observed_indices, which(model_R == 1))
      imp_local = match(imp_indices, which(model_R == 1))
      
      # Impute
      x_obs = as.numeric(dat[i, observed_indices])
      mu_obs = mu_star[observed_local]
      mu_imp = mu_star[imp_local]
      
      Sigma_oo = Sigma_star[observed_local, observed_local, drop = FALSE]
      Sigma_mo = Sigma_star[imp_local, observed_local, drop = FALSE]
      Sigma_mm = Sigma_star[imp_local, imp_local, drop = FALSE]
      
      if (length(observed_local) == 0) {
        # no observed vars in this model pattern: marginal
        cond_mean = mu_imp
        cond_var = Sigma_mm
      } else {
        cond_mean = mu_imp + Sigma_mo %*% solve(Sigma_oo) %*% (x_obs - mu_obs)
        cond_var = Sigma_mm - Sigma_mo %*% solve(Sigma_oo) %*% t(Sigma_mo)
      }
      
      imputed_value = mvrnorm(1, mu = cond_mean, Sigma = cond_var)
      X_imputed[i, imp_indices] = imputed_value
    }
    
  }
  
  return(X_imputed)
}


g_mmg = function(data, graph, m = 1, parallel = FALSE) {
  # ---- Check inputs ----
  if (!(is.matrix(data) || is.data.frame(data))) {
    stop("Data should be a matrix or data frame", call. = FALSE)
  }
  if (ncol(data) < 2) stop("Data should contain at least two columns", call. = FALSE)
  if (!is.numeric(m)) stop("Argument 'm' must be numeric", call. = FALSE)
  m = floor(m)
  if (m < 1L) stop("Number of imputations (m) must be at least 1.", call. = FALSE)
  
  data = as.data.frame(data)
  d = ncol(data)
  graph = coerce_graph(graph, data)
  var_names = colnames(data)
  
  # ---- Derive patterns ----
  pattern_strings = apply(data, 1, function(x) paste0(ifelse(is.na(x), "0", "1"), collapse = ""))
  
  missingpattern = data %>%
    mutate(pattern = pattern_strings) %>%
    count(pattern) %>%
    arrange(desc(n))
  
  model_pat_list = model_patterns_list(missing_patterns = missingpattern,
                                       graph_structure = graph, d = d, var_names = var_names)
  model_pat_df = model_patterns_df(model_pat_list, d = d, var_names = var_names)
  
  unique_model_patterns = model_pat_df %>%
    dplyr::select(-orig_pattern) %>%
    distinct() %>%
    mutate(unique_id = 1:nrow(.))
  
  patterns_join = model_pat_df %>%
    left_join(unique_model_patterns, by = var_names)
  
  observed_units = extract_observed_units(data, d, unique_model_patterns, var_names = var_names)
  
  patterns_cleaned = patterns_join %>%
    dplyr::select(-unique_id) %>%
    mutate(model_pattern = apply(.[, -1], 1, paste0, collapse = "")) %>%
    dplyr::select(orig_pattern, model_pattern)
  
  est.MLE = lapply(seq_along(observed_units), function(k) {
    X = observed_units[[k]]
    list(mu.hat = colMeans(X, na.rm = TRUE), Sigma.hat = cov(X))
  })
  
  
  for (k in seq_along(est.MLE)) {
    if (anyNA(est.MLE[[k]]$mu.hat)) {
      cat("NaN in mu.hat for pattern", k, "\n")
    }
    if (anyNA(est.MLE[[k]]$Sigma.hat)) {
      cat("NaN in Sigma.hat for pattern", k, "\n")
    }
  }
  
  return(list(
    pat = patterns_cleaned,
    unique_model_pat = unique_model_patterns,
    est.MLE = est.MLE
  ))
}

