source("../UpdateDelta.R")
source("../iteration_ma.R")
source("../iteration_procedures.R")

library(MASS)
library(snow)
library(fossil)
library(rlecuyer)

K_est_scad = vector()
K_est_mcp = vector()
K_est_ma = vector()

RI_scad = vector()
RI_mcp = vector()
RI_ma = vector()

alpha1_est_scad = vector()
alpha2_est_scad = vector()
alpha3_est_scad = vector()
alpha1_est_mcp = vector()
alpha2_est_mcp = vector()
alpha3_est_mcp = vector()
alpha1_est_ma = vector()
alpha2_est_ma = vector()
alpha3_est_ma = vector()

eta_est_scad = vector()
eta_est_mcp = vector()
eta_est_ma = vector()

alpha1_oracle = vector()
alpha2_oracle = vector()
alpha3_oracle = vector()
eta_or = vector()

a1 = c(1, 0, 0)
a2 = c(0, 1, 0)
a3 = c(0, 0, 1)

sd_alpha1_oracle = vector()
sd_alpha2_oracle = vector()
sd_alpha3_oracle = vector()


NN = 200

timestart = Sys.time()
## ---- circling -----------------
n = 90
q = 3
p = 1
rho_z = 0.3
sigma = 1
sigma_me_sqr = 0.25
corr_Z = matrix(rho_z, nrow = q, ncol = q)
diag(corr_Z) = 1
eta_true = matrix(c(1, 1, 1), nrow = q)
K_true = 3
alpha1_true = 3
alpha2_true = 0
alpha3_true = -3
# generate matrix D = {e_i - e_j, i<j}^T
I_n = diag(n)
D = matrix(nrow = choose(n, 2), ncol = n)  # choose(n, 2) = (n * (n - 1)) / 2
k = 0
for (j in 2:n) {
  for (i in 1:(j - 1)) {
    k = k + 1
    D[k, ] = I_n[i, ] - I_n[j, ]
  }
}

# generate matrix A = D \otimes I_p
A = kronecker(D, diag(p))   # (n * (n - 1)) / 2 * p by n * p

# tuning parameter for BIC
Cn = 5 * log(log(n * p + q))

cl <- makeCluster(6)
clusterExport(cl, c('n', 'p', 'q', 'A'))
clusterExport(cl,c('L2Norm','S','delta_SCAD','delta_MCP','UpdateDelta'))
set.seed(20)
clusterSetupRNG(cl, type = "RNGstream", seed=c(1, 22, 333, 444, 55, 6))

for (e in 1:NN) {
  print(e)
  # -----------------------------------------------------------------------------------------------------
  Z = mvrnorm(n, rep(1, q), corr_Z)
  
  w = matrix(rep(1, n))
  # w = cbind(rep(1, n), rnorm(n, 0, sigma))
  # u = rnorm(n)
  
  
  xi_z1 = rnorm(n, mean = 0, sd = sqrt(sigma_me_sqr))
  xi_z2 = rnorm(n, mean = 0, sd = sqrt(sigma_me_sqr))
  xi_z3 = rnorm(n, mean = 0, sd = sqrt(sigma_me_sqr))
  xi_z4 = rnorm(n, mean = 0, sd = sqrt(sigma_me_sqr))
  xi_z5 = rnorm(n, mean = 0, sd = sqrt(sigma_me_sqr))
  xi_z6 = rnorm(n, mean = 0, sd = sqrt(sigma_me_sqr))
  
  
  
  Z1 = cbind(Z[, 1]+ xi_z5, Z[, 2] + xi_z1, Z[, 3] + xi_z2)
  Z2 = cbind(Z[, 1]+ xi_z6, Z[, 2] + xi_z3, Z[, 3] + xi_z4)
  
  w1 = w
  w2 = w
  
  ## W = diag(w_1^T, w_2^T, ... , w_n^T)  n by n * p
  W = matrix(0, nrow = n, ncol = n * p)
  for (i in 1:n) {
    W[i, ((i - 1) * p + 1):(i * p)] = w[i, ]
  }
  
  ## W1, W2
  W1 = matrix(0, nrow = n, ncol = n * p)
  for (i in 1:n) {
    W1[i, ((i - 1) * p + 1):(i * p)] = w1[i, ]
  }
  
  W2 = matrix(0, nrow = n, ncol = n * p)
  for (i in 1:n) {
    W2[i, ((i - 1) * p + 1):(i * p)] = w2[i, ]
  }
  
  a = c(rep(0,n/3), rep(1, n/3), rep(2, n/3))
  b = sample(1:n, n, replace = FALSE)
  group_index_true = a[b]
  # group_index_true = matrix((Z[, 1]^2 + u < 1) * 1, nrow = n)
  real_tag = rep(group_index_true, each = p)
  beta_true = matrix(alpha1_true * (real_tag == 0) + alpha2_true * (real_tag == 1) + alpha3_true * (real_tag == 2), nrow = n)
  y = Z %*% eta_true + W %*% beta_true + rnorm(n, 0, 0.5)
  
  # generate matrix W_a
  Z_a = 0.5 * (Z1 + Z2)
  
  # generate matrix W_a
  W_a = 0.5 * (W1 + W2)
  
  # generate matrix P_Z = t(Z1)*Z2 + t(Z2)*Z1
  P_Z = t(Z1) %*% Z2 + t(Z2) %*% Z1
  
  # generate matrix W_qta = t(W1)*W2 + t(W2)*W1
  W_qta = t(W1) %*% W2 + t(W2) %*% W1
  
  # generate matrix WZ = t(W1)*Z2 + t(W2)*Z1
  WZ = t(W1) %*% Z2 + t(W2) %*% Z1
  
  # generate matrix V_ZW = W_qta - WZ*P_Z^{-1}*t(WZ)
  V_ZW = W_qta - WZ %*% solve(P_Z) %*% t(WZ)
  
  # generate matrix U_ZW = t(2*W_a) - WZ*P_Z^{-1}*t(2*Z_a)
  U_ZW = t(2 * W_a) - WZ %*% solve(P_Z) %*% t(2 * Z_a)
  
  # --------------------------------------------------------
  
  ## initial value  (It is very important to choose a good initial value!!!!!!!!!)
  lambda_star = 0.002
  Q_Z = (diag(n) - Z_a %*% solve(t(Z_a) %*% Z_a) %*% t(Z_a))
  beta_initial = solve(t(W_a) %*% Q_Z %*% W_a + lambda_star * t(A) %*% A, t(W_a) %*% Q_Z %*% y)
  eta_initial  = solve(t(Z_a) %*% Z_a, t(Z_a) %*% (y - W_a %*% beta_initial))         # (q + 1) by 1
  
  # -------------------------------------------------------------------------
  
  # calculate the oracle estimator
  ## the group information matrix
  H_qta = model.matrix(~factor(group_index_true)-1)
  H = kronecker(H_qta, diag(p))
  # -------------------------------------------------------
  M = solve(P_Z - t(WZ) %*% H %*% solve(t(H) %*% W_qta %*% H) %*% t(H) %*% WZ) %*% (t(2 * Z_a) - t(WZ) %*% H %*% solve(t(H) %*% W_qta %*% H) %*% t(2 * W_a %*% H))
  N = solve(t(H) %*% W_qta %*% H - t(H) %*% WZ %*% solve(P_Z) %*% t(WZ) %*% H) %*% (t(2 * W_a %*% H) - t(H) %*% WZ %*% solve(P_Z) %*% t(2 * Z_a))
  eta_oracle = M %*% y
  alpha_oracle = N %*% y
  alpha1_oracle = append(alpha1_oracle, alpha_oracle[1])
  alpha2_oracle = append(alpha2_oracle, alpha_oracle[2])
  alpha3_oracle = append(alpha3_oracle, alpha_oracle[3])
  eta_or = rbind(eta_or, t(eta_oracle))
  ## compute the sigma_hat
  sigma_sqr = 1 / (n - q - K_true * p) * sum((y - Z_a %*% eta_oracle - W_a %*% H %*% alpha_oracle)^2)
  ## compute the asymptotic standard error of alpha_oracle
  sd1 = sqrt(sigma_sqr  * t(a1) %*% N %*% t(N) %*% a1)
  sd2 = sqrt(sigma_sqr  * t(a2) %*% N %*% t(N) %*% a2)
  sd3 = sqrt(sigma_sqr  * t(a3) %*% N %*% t(N) %*% a3)
  sd_alpha1_oracle = append(sd_alpha1_oracle, sd1)
  sd_alpha2_oracle = append(sd_alpha2_oracle, sd2)
  sd_alpha3_oracle = append(sd_alpha3_oracle, sd3)
  
  ## Solution Path
  lambda_min = 1/n + 0.02
  # lambda_max = max(abs(t(cbind(Z_a, W_a)) %*% y) / n)
  lambda_max = 1.5
  k = 20
  lambda_list = exp(seq(log(lambda_min), log(lambda_max), length.out = k))
  # lambda_list = c(0.45, 0.5, 0.55, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.75, 0.8, 1)

  
  
  clusterExport(cl, c('Z_a', "W_a", "P_Z", "W_qta", "WZ", "V_ZW", "U_ZW", "Q_Z"))
  
  ## --------------------SCAD --------------------------------------------------
  beta_hat = matrix(nrow = n * p, ncol = k)
  eta_hat = matrix(nrow = q, ncol = k)
  delta_hat = matrix(nrow = choose(n, 2) * p, ncol = k)
  K_hat = matrix(nrow = 1, ncol = k)
  BIC = matrix(nrow = 1, ncol = k)
  
  results <- parLapply(cl, lambda_list, Iteration_sgaeiv, W1, W2, Z1, Z2, y, beta_initial, eta_initial, 3, 1, 3, Cn)
  for (i in 1:k) {
    beta_hat[, i] = results[[i]]$beta
    eta_hat[, i] = results[[i]]$eta
    delta_hat[, i] = results[[i]]$delta
    K_hat[, i] = results[[i]]$K
    BIC[, i] = results[[i]]$BIC
  }
  
  matplot(lambda_list, t(beta_hat), type = 'l', main = 'SCAD',
          xlab = expression(symbol(l)), ylab = expression(symbol(m)))
  
  plot(lambda_list, BIC, type = 'l')
  # how to choose the tuning parameter lambda
  index = which.min(BIC)
  # index = which.max(BIC[2: k] - BIC[1: (k - 1)])
  K_star = K_hat[, index]
  K_est_scad = append(K_est_scad, K_star)
  eta_star = eta_hat[, index]
  beta_star = beta_hat[,index]
  delta_star = delta_hat[,index]
  # ------------------------------------------------------------
  
  # construct adjacent matrix
  delta_mat = matrix(delta_star, nrow = p, ncol = choose(n, 2))
  l = 0
  adjacent_matrix = matrix(0, nrow = n, ncol = n)
  for (j in 2:n) {
    for (i in 1:(j - 1)) {
      l = l + 1
      if (L2Norm(delta_mat[, l]) == 0) {adjacent_matrix[i, j] = 1}
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
  group_index_est <- vector(length = n)
  for(i in 1:length(group_index)){
    group_index_est[1:n %in% group_index[[i]]] <- i-1
  }
  RI = rand.index(group_index_est,group_index_true)
  RI_scad = append(RI_scad, RI)
  ## judge whether K_star = 2
  if(K_star == 3) {
    # compute the grouped beta ----------
    alpha = matrix(nrow = p, ncol = K)
    beta_mat = matrix(beta_star, nrow = p)
    
    MeanBeta = function(beta_mat, index) {
      if(length(index) == 1){
        beta_mat[, index]
      } else {
        rowMeans(matrix(beta_mat[, index], nrow = p))
      }
    }
    
    alpha = sapply(group_index, MeanBeta, beta_mat = beta_mat)
    b = vector()
    if(alpha[1] > 1 && abs(alpha[2]) < 1) {
      b = c(alpha[1], alpha[2], alpha[3])
    } else if(alpha[1] > 1 && abs(alpha[3]) < 1){
      b = c(alpha[1], alpha[3], alpha[2])
    } else if(alpha[2] > 1 && abs(alpha[3]) < 1) {
      b = c(alpha[2], alpha[3], alpha[1])
    } else if(alpha[2] > 1 && abs(alpha[1]) < 1) {
      b = c(alpha[2], alpha[1], alpha[3])
    } else if(alpha[3] > 1 && abs(alpha[2]) < 1) {
      b = c(alpha[3], alpha[2], alpha[1])
    } else if(alpha[3] > 1 && abs(alpha[1]) < 1) {
      b = c(alpha[3], alpha[1], alpha[2])
    }
    alpha1_est_scad = append(alpha1_est_scad, b[1])
    alpha2_est_scad = append(alpha2_est_scad, b[2])
    alpha3_est_scad = append(alpha3_est_scad, b[3])
    
    eta_est_scad = rbind(eta_est_scad, eta_star)
  }
  
  ## --------------------MCP --------------------------------------------------
  beta_hat = matrix(nrow = n * p, ncol = k)
  eta_hat = matrix(nrow = q, ncol = k)
  delta_hat = matrix(nrow = choose(n, 2) * p, ncol = k)
  K_hat = matrix(nrow = 1, ncol = k)
  BIC = matrix(nrow = 1, ncol = k)
  
  results <- parLapply(cl, lambda_list, Iteration_sgaeiv, W1, W2, Z1, Z2, y, beta_initial, eta_initial, 3, 1, 2, Cn)
  for (i in 1:k) {
    beta_hat[, i] = results[[i]]$beta
    eta_hat[, i] = results[[i]]$eta
    delta_hat[, i] = results[[i]]$delta
    K_hat[, i] = results[[i]]$K
    BIC[, i] = results[[i]]$BIC
  }
  
  matplot(lambda_list, t(beta_hat), type = 'l', main = 'MCP',
          xlab = expression(symbol(l)), ylab = expression(symbol(m)))
  
  # how to choose the tuning parameter lambda
  index = which.min(BIC)
  # index = which.max(BIC[2: k] - BIC[1: (k - 1)])
  K_star = K_hat[, index]
  K_est_mcp = append(K_est_mcp, K_star)
  eta_star = eta_hat[, index]
  beta_star = beta_hat[,index]
  delta_star = delta_hat[,index]
  # ------------------------------------------------------------

  # construct adjacent matrix
  delta_mat = matrix(delta_star, nrow = p, ncol = choose(n, 2))
  l = 0
  adjacent_matrix = matrix(0, nrow = n, ncol = n)
  for (j in 2:n) {
    for (i in 1:(j - 1)) {
      l = l + 1
      if (L2Norm(delta_mat[, l]) == 0) {adjacent_matrix[i, j] = 1}
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
  group_index_est <- vector(length = n)
  for(i in 1:length(group_index)){
    group_index_est[1:n %in% group_index[[i]]] <- i-1
  }
  RI = rand.index(group_index_est,group_index_true)
  RI_mcp = append(RI_mcp, RI)
  ## judge whether K_star = 2
  if(K_star == 3) {
    # compute the grouped beta ----------
    alpha = matrix(nrow = p, ncol = K)
    beta_mat = matrix(beta_star, nrow = p)
    
    MeanBeta = function(beta_mat, index) {
      if(length(index) == 1){
        beta_mat[, index]
      } else {
        rowMeans(matrix(beta_mat[, index], nrow = p))
      }
    }
    
    alpha = sapply(group_index, MeanBeta, beta_mat = beta_mat)
    b = vector()
    if(alpha[1] > 1 && abs(alpha[2]) < 1) {
      b = c(alpha[1], alpha[2], alpha[3])
    } else if(alpha[1] > 1 && abs(alpha[3]) < 1){
      b = c(alpha[1], alpha[3], alpha[2])
    } else if(alpha[2] > 1 && abs(alpha[3]) < 1) {
      b = c(alpha[2], alpha[3], alpha[1])
    } else if(alpha[2] > 1 && abs(alpha[1]) < 1) {
      b = c(alpha[2], alpha[1], alpha[3])
    } else if(alpha[3] > 1 && abs(alpha[2]) < 1) {
      b = c(alpha[3], alpha[2], alpha[1])
    } else if(alpha[3] > 1 && abs(alpha[1]) < 1) {
      b = c(alpha[3], alpha[1], alpha[2])
    }
    alpha1_est_mcp = append(alpha1_est_mcp, b[1])
    alpha2_est_mcp = append(alpha2_est_mcp, b[2])
    alpha3_est_mcp = append(alpha3_est_mcp, b[3])
    
    eta_est_mcp = rbind(eta_est_mcp, eta_star)
  }
  
  # -----------------Ma and Huang method------------------------------------
  beta_hat = matrix(nrow = n * p, ncol = k)
  eta_hat = matrix(nrow = q, ncol = k)
  delta_hat = matrix(nrow = choose(n, 2) * p, ncol = k)
  K_hat = matrix(nrow = 1, ncol = k)
  BIC = matrix(nrow = 1, ncol = k)

  results <- parLapply(cl, lambda_list, Iteration_Ma, W_a, Z_a, y, beta_initial, eta_initial, 3, 1, 2, Cn)
  for (i in 1:k) {
    beta_hat[, i] = results[[i]]$beta
    eta_hat[, i] = results[[i]]$eta
    delta_hat[, i] = results[[i]]$delta
    K_hat[, i] = results[[i]]$K
    BIC[, i] = results[[i]]$BIC
  }

  # how to choose the tuning parameter lambda
  index = which.min(BIC)
  # index = which.max(BIC[2: k] - BIC[1: (k - 1)])
  K_star = K_hat[, index]
  K_est_ma = append(K_est_ma, K_star)
  eta_star = eta_hat[, index]
  beta_star = beta_hat[,index]
  delta_star = delta_hat[,index]
  # ------------------------------------------------------------
  ## judge whether K_star = 2
  # construct adjacent matrix
  delta_mat = matrix(delta_star, nrow = p, ncol = choose(n, 2))
  l = 0
  adjacent_matrix = matrix(0, nrow = n, ncol = n)
  for (j in 2:n) {
    for (i in 1:(j - 1)) {
      l = l + 1
      if (L2Norm(delta_mat[, l]) == 0) {adjacent_matrix[i, j] = 1}
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
  group_index_est <- vector(length = n)
  for(i in 1:length(group_index)){
    group_index_est[1:n %in% group_index[[i]]] <- i-1
  }
  RI = rand.index(group_index_est,group_index_true)
  RI_ma = append(RI_ma, RI)
  ## judge whether K_star = 2
  if(K_star == 2) {
    # compute the grouped beta ----------
    alpha = matrix(nrow = p, ncol = K)
    beta_mat = matrix(beta_star, nrow = p)
    
    MeanBeta = function(beta_mat, index) {
      if(length(index) == 1){
        beta_mat[, index]
      } else {
        rowMeans(matrix(beta_mat[, index], nrow = p))
      }
    }
    
    alpha = sapply(group_index, MeanBeta, beta_mat = beta_mat)
    b = vector()
    if(alpha[1] > 1 && abs(alpha[2]) < 1) {
      b = c(alpha[1], alpha[2], alpha[3])
    } else if(alpha[1] > 1 && abs(alpha[3]) < 1){
      b = c(alpha[1], alpha[3], alpha[2])
    } else if(alpha[2] > 1 && abs(alpha[3]) < 1) {
      b = c(alpha[2], alpha[3], alpha[1])
    } else if(alpha[2] > 1 && abs(alpha[1]) < 1) {
      b = c(alpha[2], alpha[1], alpha[3])
    } else if(alpha[3] > 1 && abs(alpha[2]) < 1) {
      b = c(alpha[3], alpha[2], alpha[1])
    } else if(alpha[3] > 1 && abs(alpha[1]) < 1) {
      b = c(alpha[3], alpha[1], alpha[2])
    }
    alpha1_est_ma = append(alpha1_est_ma, b[1])
    alpha2_est_ma = append(alpha2_est_ma, b[2])
    alpha3_est_ma = append(alpha3_est_ma, b[3])
    
    eta_est_ma = rbind(eta_est_ma, eta_star)
  }

}
stopCluster(cl)
## circling -------------------------------------------------------
timeend = Sys.time()
timeend-timestart
#
save.image(file = "example2.RData")
