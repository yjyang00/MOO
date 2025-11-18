mbpe_mmg = function(X, Y, K, graph, em_num_init = 5, em_max_iter = 100, em_tol = 1e-5, m, s_max = NULL){
  
  d = ncol(Y)
  p = ncol(X)
  graph = coerce_graph(graph, Y)
  var_names = colnames(Y)
  
  # ---- Derive patterns ----
  res_vec = matrix(as.integer(!is.na(Y) & Y >= 0 & Y <= s_max[col(Y)]), nrow(Y), d)
  missingpattern = data.frame(pattern = apply(res_vec, 1, \(x) paste0(ifelse(is.na(x), 0, x), collapse = ""))) %>%
    count(pattern) %>%
    arrange(desc(n))
  
  model_pat_list = model_patterns_list(missingpattern, graph, d, var_names)
  model_pat_df = model_patterns_df(model_pat_list, d, var_names = var_names)
  
  unique_model_patterns = model_pat_df %>%
    dplyr::select(-orig_pattern) %>%
    distinct() %>%
    mutate(unique_id = dplyr::row_number())
  
  patterns_join = model_pat_df %>% left_join(unique_model_patterns, by = var_names)
  
  patterns_cleaned = patterns_join %>%
    dplyr::select(-unique_id) %>%
    mutate(model_pattern = apply(.[, -1], 1, paste0, collapse = "")) %>%
    dplyr::select(orig_pattern, model_pattern)
  
  data2 = cbind(X, Y, orig_pattern=apply(res_vec, 1, \(x) paste0(x, collapse="")))
  
  observed_units = lapply(1:nrow(unique_model_patterns), function(i) {
    mask = as.logical(unique_model_patterns[i, var_names])
    rows = complete.cases(Y[, mask, drop = FALSE])
    Xi = X[rows, , drop = FALSE]
    Yi = Y[rows, mask, drop = FALSE]
    cbind(Xi, Yi)
  })
  
  
  # ---- EM ----
  res_EM = list()
  for (i in 1:nrow(unique_model_patterns)) {
    Ui = observed_units[[i]]
    Xi = Ui[, seq_len(p), drop = FALSE]
    Yi = Ui[, -seq_len(p), drop = FALSE]
    s_max_i = s_max[which(unique_model_patterns[i, -ncol(unique_model_patterns)] == 1)]
    
    # EM_fit_serial is from the mixturebpe package (more details: https://github.com/danielsuen/mixturebpe)
    res_EM[[i]] = EM_fit_serial(
      X = Xi,
      Y = Yi,
      num_initializations = em_num_init,
      K = K,
      s_max = s_max_i,
      M_em = em_max_iter,
      tol = em_tol,
      printFlag = FALSE
    )
  }
  
  return(list(
    pat = patterns_cleaned,
    unique_model_pat = unique_model_patterns,
    est.MLE = res_EM
  ))
}

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
