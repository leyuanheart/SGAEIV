source("../UpdateDelta.R")
source("../iteration_ma.R")
source("../iteration_procedures.R")


library(parallel)
# library(dummies) # for R version 3.6.0
library(dummy)


load('nine_dat.RData')

# ======================== Training ===================================


n = length(nine_dat$ID)
q = 6
p = 1

Z1 = cbind(rep(1, n),
           nine_dat$Age_cal,
           nine_dat$race, 
           nine_dat$college, 
           nine_dat$SBP1,
           nine_dat$DBP1)
# scale(nine_dat$SBP1, center = F, scale = T),
# scale(nine_dat$DBP1, center = F, scale = T))
#nine_dat$SBP1_change, 
#nine_dat$DBP1_change)

Z2 = cbind(rep(1, n),
           nine_dat$Age_cal,
           nine_dat$race, 
           nine_dat$college, 
           nine_dat$SBP3,
           nine_dat$DBP3)
# scale(nine_dat$SBP3, center = F, scale = T),
# scale(nine_dat$DBP3, center = F, scale = T))
#nine_dat$SBP3_change, 
#nine_dat$DBP3_change)


group <- nine_dat$Group
# g <- dummy(group)
# g <- g[, -1]
g <- dummy(data.frame(factor(group)))
g <- as.numeric(g[, -1]) - 1

#w1 = cbind(rep(1, n), g)
#w2 = cbind(rep(1, n), g)

w1 = matrix(g, nrow = n)
w2 = matrix(g, nrow = n)


##------------------------


W1 = matrix(0, nrow = n, ncol = n * p)
for (i in 1:n) {
  W1[i, ((i - 1) * p + 1):(i * p)] = w1[i, ]
}

W2 = matrix(0, nrow = n, ncol = n * p)
for (i in 1:n) {
  W2[i, ((i - 1) * p + 1):(i * p)] = w2[i, ]
}

y = as.vector(nine_dat$BMI_change)

##----------------------------------------------------------------------------
## related matrix
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
Q_Z = (diag(n) - Z_a %*% solve(t(Z_a) %*% Z_a) %*% t(Z_a))
lambda_star = 0.002
beta_initial = solve(t(W_a) %*% Q_Z %*% W_a + lambda_star * t(A) %*% A, t(W_a) %*% Q_Z %*% y) 
eta_initial  = solve(t(Z_a) %*% Z_a, t(Z_a) %*% (y - W_a %*% beta_initial))         # (q + 1) by 1
# ---------------------------------------------------------------------------


## Solution Path for MCP
lambda_min = 0.01
# lambda_max = max(abs(t(cbind(Z_a, W_a)) %*% y) / n)
lambda_max = 0.09
k = 15
lambda_list = exp(seq(log(lambda_min), log(lambda_max), length.out = k))
beta_hat = matrix(nrow = n * p, ncol = k)
eta_hat = matrix(nrow = q, ncol = k)
delta_hat = matrix(nrow = choose(n, 2) * p, ncol = k)
K_hat = matrix(nrow = 1, ncol = k)
BIC = matrix(nrow = 1, ncol = k)


cl <- makeCluster(detectCores() - 2)
clusterExport(cl, c('n', 'p', 'q', 'A'))
clusterExport(cl,c('L2Norm','S','delta_lasso', 'delta_SCAD','delta_MCP','UpdateDelta'))
clusterExport(cl, c('Z_a', "W_a", "P_Z", "W_qta", "WZ", "V_ZW", "U_ZW", "Q_Z"))

Cn = log(n * p + q)

results_ma <- parLapply(cl, lambda_list, Iteration_Ma, W_a, Z_a, y, beta_initial, eta_initial, 3, 1, 2, Cn)

results_mcp <- parLapply(cl, lambda_list, Iteration_sgaeiv, W1, W2, Z1, Z2, y, beta_initial, eta_initial, 3, 1, 2, Cn)
results_scad <- parLapply(cl, lambda_list, Iteration_sgaeiv, W1, W2, Z1, Z2, y, beta_initial, eta_initial, 3, 1, 3, Cn)

stopCluster(cl)

# if you want to get the results of the proposed method using scad 
# or the method in \cite{Ma et al., 2020}
# uncomment the correspond lines.

results = results_mcp
# results = results_scad
# results = results_ma



for (i in 1:k) {
  beta_hat[, i] = results[[i]]$beta
  eta_hat[, i] = results[[i]]$eta
  delta_hat[, i] = results[[i]]$delta
  K_hat[, i] = results[[i]]$K
  BIC[, i] = results[[i]]$BIC
}

# png(file="path.png")
matplot(lambda_list, t(beta_hat), type = 'l', las=1, col='black', #main = 'MCP',
        xlab = expression(symbol(l)), ylab = expression(symbol(b)))
# dev.off()
# dev.new()
# png(file="lambda.png")
plot(lambda_list, BIC, type = 'l', main = 'BIC')
# dev.off()

# ------------------------------------------------------------------------------------

# how to choose the tuning parameter lambda
index = which.min(BIC)
# index = which.max(BIC[2: k] - BIC[1: (k - 1)])
K_star = K_hat[, index]
eta_star = eta_hat[, index]
beta_star = beta_hat[,index]
delta_star = delta_hat[,index]

# ------------------------------------------------------------
# construct adjacent matrix
delta_mat = matrix(delta_star, nrow = p, ncol = choose(n, 2))
l = 0
adjacent_matrix = matrix(0, nrow = n, ncol = n)
for (j in 2:n) {
  for (h in 1:(j - 1)) {
    l = l + 1
    if (L2Norm(delta_mat[, l]) == 0) {adjacent_matrix[h, j] = 1}
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

# ======================================================================


save.image(file = "real_data.RData")
