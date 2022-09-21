# iteration procedures 
Iteration_sgaeiv = function(lambda, W1, W2, Z1, Z2, y, beta_initial, eta_initial, gamma = 3, vartheta = 1, penalty, Cn = log(n * p + q)) {
  # Algorithm of ADMM for concave fusion
  max_tol = 10^(-3)
  iter = 0
  diff = 1
  beta = beta_initial
  eta = eta_initial
  upsilon = matrix(0, nrow = choose(n, 2) * p, ncol = 1)
  while(diff > max_tol & iter < 1000) {
    iter = iter + 1
    zeta = A %*% beta + 1 / vartheta * upsilon
    zeta_mat = matrix(zeta, nrow = p, ncol = choose(n, 2))
    delta = UpdateDelta(zeta_mat, gamma = gamma, lambda = lambda, vartheta = vartheta, penalty = penalty)
    beta = solve(V_ZW + 2 * vartheta * t(A) %*% A, U_ZW %*% y + 2 * vartheta * t(A) %*% (delta - 1 / vartheta * upsilon))
    eta = solve(P_Z, t(2 * Z_a) %*% y - t(WZ) %*% beta)
    upsilon = upsilon + vartheta * (A %*% beta - delta)
    diff = L2Norm(A %*% beta - delta)
    # print(paste('iter = ', iter, 'diff = ', diff, sep = ''))
  }
  # construct adjacent matrix
  delta_mat = matrix(delta, nrow = p, ncol = choose(n, 2))
  k = 0
  adjacent_matrix = matrix(0, nrow = n, ncol = n)
  for (j in 2:n) {
    for (i in 1:(j - 1)) {
      k = k + 1
      if (L2Norm(delta_mat[, k]) == 0) {adjacent_matrix[i, j] = 1}
    }
  }
  # compute the number of groups
  K = 0
  res_index = 1:n
  group_index = list()
  while(length(res_index) != 0) {
    K = K + 1
    group_index[[K]] = c(res_index[1], which(adjacent_matrix[res_index[1], ] == 1))
    res_index = setdiff(res_index, group_index[[K]])
  }
  # compute the grouped beta ----------
  # alpha = matrix(nrow = p, ncol = K)
  # MeanBeta = function(beta_mat, index) {
  #   ifelse(length(index) == 1, beta_mat[, index], colMeans(matrix(beta_mat[, index], nrow = p)))
  # }
  # alpha = sapply(group_index, MeanBeta, beta_mat = temp)
  
  # compute the modified Bayes Information Criterion
  BIC = log(L2Norm(y - Z_a %*% eta - W_a %*% beta)^2 / n) + Cn * log(n) / n * (K * p + q)

  return(list(beta = beta, eta = eta, delta = delta, K = K, BIC = BIC))
  #return(list(beta = beta, eta = eta, delta = delta, K = K))
}
