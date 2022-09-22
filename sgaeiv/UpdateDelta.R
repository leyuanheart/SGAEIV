# L2 norm
L2Norm = function(x) {
  sqrt(sum((x)^2))
}

# groupwise thresholding operator
S = function(z, t) {
  max(0, 1 - t / L2Norm(z)) * z  #apply(z, 2, function(z) {max(0, 1 - t / L2Norm(z)) * z})
}

# iteration for variable delta under lasso penalty
delta_lasso = function(zeta_mat, lambda, vartheta = 1) {
  t = lambda / vartheta
  delta_lasso_mat = apply(zeta_mat, 2, S, t = t)
  c(delta_lasso_mat)
}

# iteration for variable delta under MCP penalty
delta_MCP = function(zeta_mat, gamma = 3, lambda, vartheta = 1) {
  t = lambda / vartheta
  l2_norm = apply(zeta_mat, 2, L2Norm)
  threshold = gamma * lambda
  delta_MCP_mat = zeta_mat
  delta_MCP_mat[, l2_norm <= threshold] =  apply(matrix(delta_MCP_mat[, l2_norm <= threshold], nrow = p), 2, S, t = t) / (1 - (1 / (gamma * vartheta)))
  c(delta_MCP_mat)
}

# iteration for variable delta under SACD penalty
delta_SCAD = function(zeta_mat, gamma = 3, lambda, vartheta = 1) {
  t1 = lambda / vartheta
  t2 = (gamma * lambda) / ((gamma - 1) * vartheta)
  threshold1 = lambda + lambda / vartheta
  threshold2 = gamma * lambda
  l2_norm = apply(zeta_mat, 2, L2Norm)
  delta_SCAD_mat = zeta_mat
  delta_SCAD_mat[, l2_norm <= threshold1] =  apply(matrix(delta_SCAD_mat[, l2_norm <= threshold1], nrow = p), 2, S, t = t1)
  delta_SCAD_mat[, (l2_norm > threshold1 & l2_norm <= threshold2)] =  apply(matrix(delta_SCAD_mat[, (l2_norm > threshold1 & l2_norm <= threshold2)], nrow = p), 
                                                                            2, S, t = t2) / (1 - 1 / ((gamma - 1) * vartheta))
  c(delta_SCAD_mat)
}

# iteration for variable delta
UpdateDelta = function(zeta_mat, gamma = 3, lambda, vartheta = 1, penalty) {
  switch(penalty, delta_lasso(zeta_mat = zeta_mat, lambda = lambda, vartheta = vartheta), 
         delta_MCP(zeta_mat = zeta_mat, gamma = gamma, lambda = lambda, vartheta = vartheta), 
         delta_SCAD(zeta_mat = zeta_mat, gamma = gamma, lambda = lambda, vartheta = vartheta))
  # penalty = 1 represents lasso, penalty = 2 represents MCP, penalty = 3 represents SCAD
}