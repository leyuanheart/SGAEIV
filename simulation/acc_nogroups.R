source("../sgaeiv/UpdateDelta.R")
source("../sgaeiv/iteration_procedures.R")
library(MASS)
# library(dummies)
library(fossil)
library(rlecuyer)
library(snow)
# ----------------------------------------------------------------
NN = 200
K_est_scad = vector(length = NN)
K_est_mcp = vector(length = NN)
# K_est_ma = vector(length = NN)
# fcr_scad = vector()
# fcr_mcp = vector()
# fcr_ma = vector()
beta_est_scad = vector()
beta_est_mcp = vector()
# alpha1_est_ma = vector()
# alpha2_est_ma = vector()
eta_est_scad = vector()
eta_est_mcp = vector()
# eta_est_ma = vector()
eta_oracle = vector()
beta_oracle = vector()
sd_beta_scad = vector()
sd_eta_scad = vector()
sd_beta_mcp = vector()
sd_eta_mcp = vector()
sd_beta_oracle = vector()
sd_eta_oracle = vector()

n = 150
q = 3
p = 1
rho_z = 0.3
sigma_z = 1
sigma_w = 1
sigma = 0.25
corr_Z = matrix(rho_z, nrow = q-1, ncol = q-1)
diag(corr_Z) = 1
eta_true = matrix(c(1, 1.5, 2), nrow = q)
beta_true = 2

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
Cn = 2 * log(n * p + q)
# -------------------------------------------------------------------------------------------------------
timestart = Sys.time()
cl <- makeCluster(6)
clusterExport(cl, c('n', 'p', 'q', 'A'))
clusterExport(cl,c('L2Norm','S','delta_SCAD','delta_MCP','UpdateDelta'))
set.seed(30)
clusterSetupRNG(cl, type = "RNGstream", seed=c(1, 22, 333, 444, 55, 6))
for (e in 1:NN) {
  print(e)
  # -----------------------------------------------------------------------------------------------------
  z = mvrnorm(n, rep(1, q-1), sigma_z * corr_Z)   # generate multivariate normal random numbers
  Z = cbind(rep(1, n), z)
  w = matrix(rnorm(n, mean = 1, sd = sigma_w), nrow = n)
  
  xi_z11 = rnorm(n, mean = 0, sd = sqrt(sigma))
  xi_z12 = rnorm(n, mean = 0, sd = sqrt(sigma))
  xi_z21 = rnorm(n, mean = 0, sd = sqrt(sigma))
  xi_z22 = rnorm(n, mean = 0, sd = sqrt(sigma))
  xi_w1 = rnorm(n, mean = 0, sd = sqrt(sigma))
  xi_w2 = rnorm(n, mean = 0, sd = sqrt(sigma))
  
  Z1 = cbind(Z[, 1], Z[, 2] + xi_z11, Z[, 3] + xi_z12)
  Z2 = cbind(Z[, 1], Z[, 2] + xi_z21, Z[, 3] + xi_z22)
  
  w1 = w + xi_w1
  w2 = w + xi_w2
  
  ## 构造W = diag(w_1^T, w_2^T, ... , w_n^T)  n by n * p
  W = matrix(0, nrow = n, ncol = n * p)
  for (i in 1:n) {
    W[i, ((i - 1) * p + 1):(i * p)] = w[i, ]
  }
  
  ## 构造W1, W2
  W1 = matrix(0, nrow = n, ncol = n * p)
  for (i in 1:n) {
    W1[i, ((i - 1) * p + 1):(i * p)] = w1[i, ]
  }
  
  W2 = matrix(0, nrow = n, ncol = n * p)
  for (i in 1:n) {
    W2[i, ((i - 1) * p + 1):(i * p)] = w2[i, ]
  }
  
  
  y = Z %*% eta_true + w %*% beta_true + rnorm(n, 0, 0.5)
  # -----------------------------------------------------------------------------
  ## related matrix
  
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
  
  # --------------------------------------------------------------------
  ## initial value  (It is very important to choose a good initial value!!!!!!!!!)
  Q_Z = (diag(n) - Z_a %*% solve(t(Z_a) %*% Z_a) %*% t(Z_a))
  lambda_star = 0.002
  beta_initial = solve(t(W_a) %*% Q_Z %*% W_a + lambda_star * t(A) %*% A, t(W_a) %*% Q_Z %*% y) 
  eta_initial  = solve(t(Z_a) %*% Z_a, t(Z_a) %*% (y - W_a %*% beta_initial))         # (q + 1) by 1
  # -----------------------------------------------------------------------
  # calculate the oracle estimator
  
  w_a = 0.5 * (w1 + w2)
  
  w_qta = t(w1) %*% w2 + t(w2) %*% w1
  
  wZ = t(w1) %*% Z2 + t(w2) %*% Z1
  
  # ------------------------------------------------------------------------
  M = 2*solve(P_Z-t(wZ)%*%solve(w_qta)%*%wZ)%*%(t(Z_a)-t(wZ)%*%solve(w_qta)%*%t(w_a))
  N = 2*solve(w_qta-wZ%*%solve(P_Z)%*%t(wZ))%*%(t(w_a)-wZ%*%solve(P_Z)%*%t(Z_a))
  eta_or = M %*% y
  beta_or = N %*% y
  
  eta_oracle = rbind(eta_oracle, c(eta_or))
  beta_oracle = append(beta_oracle, beta_or)
  
  ## compute the sigma_hat
  sigma_hat = 1 / (n - q - p) * sum((y - Z_a %*% eta_or - w_a %*% beta_or)^2)
  ## compute the asymptotic standard error of beta_oracle, eta_oracle
  sd_eta_oracle = rbind(sd_eta_oracle, sqrt(diag(sigma_hat * M %*% t(M))))
  sd_beta_oracle = append(sd_beta_oracle, sqrt(sigma_hat * N %*% t(N)))
  # --------------------------------------------
  ## Solution Path
  lambda_min = 1/n + 0.02
  # lambda_max = max(abs(t(cbind(Z_a, W_a)) %*% y) / n)
  lambda_max = 1.5
  k = 15
  lambda_list = exp(seq(log(lambda_min), log(lambda_max), length.out = k))
  
  
  # ----------parallel
  clusterExport(cl, c('Z_a', "W_a", "P_Z", "W_qta", "WZ", "V_ZW", "U_ZW", 'Q_Z'))
  # -----------------------------MCP----------------------------------------------------------
  beta_hat = matrix(0, nrow = n * p, ncol = k)
  eta_hat = matrix(0, nrow = q, ncol = k)
  delta_hat = matrix(0, nrow = choose(n, 2) * p, ncol = k)
  K_hat = matrix(0, nrow = 1, ncol = k)
  BIC = matrix(0, nrow = 1, ncol = k)
  
  results <- parLapply(cl, lambda_list, Iteration_sgaeiv, W1, W2, Z1, Z2, y, beta_initial, eta_initial, 3, 1, 2, Cn)
  for (i in 1:k) {
    beta_hat[, i] = results[[i]]$beta
    eta_hat[, i] = results[[i]]$eta
    delta_hat[, i] = results[[i]]$delta
    K_hat[, i] = results[[i]]$K
    BIC[, i] = results[[i]]$BIC
  }
  matplot(lambda_list, t(beta_hat), type = 'l', main = 'MCP',
          xlab = expression(symbol(l)), ylab = expression(symbol(b)))
  # plot(lambda_list, BIC, type = 'l', main = 'BIC')
  # -------------------------------------------------------------------------------------------
  # how to choose the tuning parameter lambda
  index = which.min(BIC)
  # index = which.max(BIC[2: k] - BIC[1: (k - 1)])
  K_star = K_hat[, index]
  K_est_mcp[e] = K_star 
  eta_star = eta_hat[, index]
  beta_star = beta_hat[,index]
  # ------------------------------------------------------------
  eta_est_mcp = rbind(eta_est_mcp, eta_star)
  ## judge whether K_star = 1
  if(K_star == 1) {
    beta_est_mcp = append(beta_est_mcp, beta_star)
    
    
    ## compute the sigma_hat
    sigma_hat = 1 / (n - q - p) * sum((y - Z_a %*% eta_star - w_a * beta_star)^2)
    ## compute the asymptotic standard error of beta_oracle, eta_oracle
    sd_eta_mcp = rbind(sd_eta_mcp, sqrt(diag(sigma_hat * M %*% t(M))))
    sd_beta_mcp = append(sd_beta_mcp, sqrt(sigma_hat * N %*% t(N)))
    
  }
  # -------------------------------------------------------------------------------------
  # -----------------------------SCAD-------------------------------------------------------
  beta_hat = matrix(0, nrow = n * p, ncol = k)
  eta_hat = matrix(0, nrow = q, ncol = k)
  delta_hat = matrix(0, nrow = choose(n, 2) * p, ncol = k)
  K_hat = matrix(0, nrow = 1, ncol = k)
  BIC = matrix(0, nrow = 1, ncol = k)
  
  results <- parLapply(cl, lambda_list, Iteration_sgaeiv, W1, W2, Z1, Z2, y, beta_initial, eta_initial, 3, 1, 3, Cn)
  for (i in 1:k) {
    beta_hat[, i] = results[[i]]$beta
    eta_hat[, i] = results[[i]]$eta
    delta_hat[, i] = results[[i]]$delta
    K_hat[, i] = results[[i]]$K
    BIC[, i] = results[[i]]$BIC
  }
  matplot(lambda_list, t(beta_hat), type = 'l', main = 'SCAD',
          xlab = expression(symbol(l)), ylab = expression(symbol(b)))
  # plot(lambda_list, BIC, type = 'l', main = 'BIC')
  # -------------------------------------------------------------------------------------------
  # how to choose the tuning parameter lambda
  index = which.min(BIC)
  # index = which.max(BIC[2: k] - BIC[1: (k - 1)])
  K_star = K_hat[, index]
  K_est_scad[e] = K_star 
  eta_star = eta_hat[, index]
  beta_star = beta_hat[,index]
  delta_star = delta_hat[,index]
  # ------------------------------------------------------------
  eta_est_scad = rbind(eta_est_scad, eta_star)
  ## judge whether K_star = 1
  if(K_star == 1) {
    beta_est_scad = append(beta_est_scad, beta_star)
    
    ## compute the sigma_hat
    sigma_hat = 1 / (n - q - p) * sum((y - Z_a %*% eta_star - w_a * beta_star)^2)
    ## compute the asymptotic standard error of beta_oracle, eta_oracle
    sd_eta_scad = rbind(sd_eta_scad, sqrt(diag(sigma_hat * M %*% t(M))))
    sd_beta_scad = append(sd_beta_scad, sqrt(sigma_hat * N %*% t(N)))
  }
}
stopCluster(cl)
timeend = Sys.time()


save.image("example3.RData")