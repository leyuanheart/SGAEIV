# ====================== Table 1 and Table 2 Start ========================================

## ============== Example 1 ==================
load('./simulation/results/example1.RData')

scad = c(mean(K_est_scad),
        median(K_est_scad),
        sd(K_est_scad),
        sum(K_est_scad == K_true)/NN,
        mean(RI_scad),
        mean(alpha1_est_scad - alpha1_true),
        sd(alpha1_est_scad),
        mean(alpha2_est_scad - alpha2_true),
        sd(alpha2_est_scad))


mcp = c(mean(K_est_mcp),
        median(K_est_mcp),
        sd(K_est_mcp),
        sum(K_est_mcp == K_true)/NN,
        mean(RI_mcp),
        mean(alpha1_est_mcp - alpha1_true),
        sd(alpha1_est_mcp),
        mean(alpha2_est_mcp - alpha2_true),
        sd(alpha2_est_mcp))
        
        
ma = c(mean(K_est_ma),
      median(K_est_ma),
      sd(K_est_ma),
      sum(K_est_ma == K_true)/NN,
      mean(RI_ma),
      mean(alpha1_est_ma - alpha1_true),
      sd(alpha1_est_ma),
      mean(alpha2_est_ma  - alpha2_true),
      sd(alpha2_est_ma))


df1 = cbind(scad, mcp, ma)
rownames(df1) = c('K_mean', 'K_median', 'K_std', 'prop', 'RI',
                  'alpha1_bias', 'alpha1_std', 'alpha2_bias', 'alpha2_std')

## ============== Example 2 ==================
load('./simulation/results/example2.RData')


scad = c(mean(K_est_scad),
         median(K_est_scad),
         sd(K_est_scad),
         sum(K_est_scad == K_true)/NN,
         mean(RI_scad),
         mean(alpha1_est_scad - alpha1_true),
         sd(alpha1_est_scad),
         mean(alpha2_est_scad - alpha2_true),
         sd(alpha2_est_scad),
         mean(alpha3_est_scad - alpha3_true),
         sd(alpha3_est_scad))


mcp = c(mean(K_est_mcp),
        median(K_est_mcp),
        sd(K_est_mcp),
        sum(K_est_mcp == K_true)/NN,
        mean(RI_mcp),
        mean(alpha1_est_mcp - alpha1_true),
        sd(alpha1_est_mcp),
        mean(alpha2_est_mcp - alpha2_true),
        sd(alpha2_est_mcp),
        mean(alpha3_est_mcp - alpha3_true),
        sd(alpha3_est_mcp))


ma = c(mean(K_est_ma),
       median(K_est_ma),
       sd(K_est_ma),
       sum(K_est_ma == K_true)/NN,
       mean(RI_ma),
       mean(alpha1_est_ma - alpha1_true),
       sd(alpha1_est_ma),
       mean(alpha2_est_ma - alpha2_true),
       sd(alpha2_est_ma),
       mean(alpha3_est_ma - alpha3_true),
       sd(alpha3_est_ma))


df2 = cbind(scad, mcp, ma)
rownames(df2) = c('K_mean', 'K_median', 'K_std', 'prop', 'RI',
                  'alpha1_bias', 'alpha1_std', 'alpha2_bias', 'alpha2_std', 'alpha3_bias', 'alpha3_std')
# ====================== Table 1 and Table 2 End ========================================



# ====================== Figure 1 Start ========================================
load('./simulation/results/example2.RData')

Alpha1Mse <- function(alpha){
  sqrt((alpha - alpha1_true)^2) / sqrt(p)
}

Alpha2Mse <- function(alpha){
  sqrt((alpha - alpha2_true)^2) / sqrt(p)
}

Alpha3Mse <- function(alpha){
  sqrt((alpha - alpha3_true)^2) / sqrt(p)
}

EtaMse <- function(eta) {
  sqrt(sum((eta - eta_true)^2)) / sqrt(q)
}

alpha1_or <- Alpha1Mse(alpha1_oracle)
alpha1_scad <- Alpha1Mse(na.omit(alpha1_est_scad))
alpha1_mcp <- Alpha1Mse(na.omit(alpha1_est_mcp))
alpha1_ma <- Alpha1Mse(alpha1_est_ma)


alpha2_or <- Alpha2Mse(alpha2_oracle)
alpha2_scad <- Alpha2Mse(na.omit(alpha2_est_scad))
alpha2_mcp <- Alpha2Mse(na.omit(alpha2_est_mcp))
alpha2_ma <- Alpha2Mse(alpha2_est_ma)


alpha3_or <- Alpha3Mse(alpha3_oracle)
alpha3_scad <- Alpha3Mse(na.omit(alpha3_est_scad))
alpha3_mcp <- Alpha3Mse(na.omit(alpha3_est_mcp))
alpha3_ma <- Alpha3Mse(alpha3_est_ma)



alpha1_bias_or = alpha1_oracle - alpha1_true
alpha1_bias_mcp = alpha1_est_mcp - alpha1_true
alpha1_bias_scad = alpha1_est_scad - alpha1_true
alpha1_bias_ma = alpha1_est_ma - alpha1_true


alpha2_bias_or = alpha2_oracle - alpha2_true
alpha2_bias_mcp = alpha2_est_mcp - alpha2_true
alpha2_bias_scad = alpha2_est_scad - alpha2_true
alpha2_bias_ma = alpha2_est_ma - alpha2_true


alpha3_bias_or = alpha3_oracle - alpha3_true
alpha3_bias_mcp = alpha3_est_mcp - alpha3_true
alpha3_bias_scad = alpha3_est_scad - alpha3_true
alpha3_bias_ma = alpha3_est_ma - alpha3_true


oldpar <- par(mfrow = c(2, 3))
boxplot(alpha1_bias_scad, alpha1_bias_mcp, alpha1_bias_ma, names = c('SCAD', 'MCP','Ma'),
        ylab = expression(paste('Bias for ', hat(symbol(a))[1])), ylim = c(-1, 1), las=1, cex.lab=1, cex.axis=1
)
boxplot(alpha2_bias_scad, alpha2_bias_mcp, alpha2_bias_ma, names = c('SCAD', 'MCP','Ma'),
        ylab = expression(paste('Bias for ', hat(symbol(a))[2])), ylim = c(-1, 1), las=1, cex.lab=1, cex.axis=1
)
boxplot(alpha3_bias_scad, alpha3_bias_mcp, alpha3_bias_ma, names = c('SCAD', 'MCP','Ma'),
        ylab = expression(paste('Bias for ', hat(symbol(a))[3])), ylim = c(-1, 1), las=1, cex.lab=1, cex.axis=1
)

boxplot(alpha1_scad, alpha1_mcp, alpha1_ma,  names = c('SCAD', 'MCP','Ma'),
        ylab = expression(paste('MSE for ', hat(symbol(a))[1])), ylim = c(0, 1), cex.lab=1, cex.axis=1
)
boxplot(alpha2_scad, alpha2_mcp, alpha2_ma, names = c('SCAD', 'MCP', 'Ma'),
        ylab = expression(paste('MSE for ', hat(symbol(a))[2])), ylim = c(0, 1), cex.lab=1,cex.axis=1
)
boxplot(alpha3_scad, alpha3_mcp, alpha3_ma, names = c('SCAD', 'MCP', 'Ma'),
        ylab = expression(paste('MSE for ', hat(symbol(a))[3])), ylim = c(0, 1), cex.lab=1,cex.axis=1
)

par(oldpar)
# ====================== Figure 1 End ========================================




# ====================== Table 3 Start ========================================
load('./simulation/results/example3.RData')

bias = apply(eta_oracle, 1, function(eta){eta - c(eta_true)})
apply(bias, 1, mean)
apply(eta_oracle, 2, sd)
apply(sd_eta_oracle, 2, mean)

bias = apply(eta_est_scad, 1, function(eta){eta - c(eta_true)})
apply(bias, 1, mean)
apply(eta_est_scad, 2, sd)
apply(sd_eta_scad, 2, mean)

bias = apply(eta_est_mcp, 1, function(eta){eta - c(eta_true)})
apply(bias, 1, mean)
apply(eta_est_mcp, 2, sd)
apply(sd_eta_mcp, 2, mean)

eta = rbind(apply(bias, 1, mean),
           apply(eta_oracle, 2, sd),
           apply(sd_eta_oracle, 2, mean),
           apply(bias, 1, mean),
           apply(eta_est_scad, 2, sd),
           apply(sd_eta_scad, 2, mean),
           apply(bias, 1, mean),
           apply(eta_est_mcp, 2, sd),
           apply(sd_eta_mcp, 2, mean))

beta = c(mean(beta_oracle - beta_true),
         sd(beta_oracle),
         mean(sd_beta_oracle),
         mean(beta_est_scad - beta_true),
         sd(beta_est_scad),
         mean(sd_beta_scad),
         mean(beta_est_mcp - beta_true),
         sd(beta_est_mcp),
         mean(sd_beta_mcp))

df3 = cbind(beta, eta)
colnames(df3) = c('beta', 'eta1', 'eta2', 'eta3')
rownames(df3) = c('oracle_bias', 'oracle_ese', 'oracle_ase',
                  'scad_bias', 'scad_ese', 'scad_ase',
                  'mcp_bias', 'mcp_ese', 'mcp_ase')
# ====================== Table 3 End ========================================



# ====================== Figure 2 Start ========================================
load('./real_data/real_data.RData')

## fit a homogeneous model
m = lm(BMI_change ~ Age_cal + factor(race) + factor(college) + factor(Group) + SBP1 + DBP1, data = nine_dat)
# summary(m)
# plot(density(m$residuals), xlab = 'residuals', main = 'homogeneous model')

## heterogeneous model =======================
x = as.matrix(nine_dat[, c(7, 9, 10, 2, 4, 6)])
x = cbind(rep(1, n), x)
theta1 = matrix(c(eta_star, alpha[1]), nrow = 7)
theta2 = matrix(c(eta_star, alpha[2]), nrow = 7)
c1 = y[group_index[[1]]] - x[group_index[[1]], ] %*% theta1 
c2 = y[group_index[[2]]] - x[group_index[[2]], ] %*% theta2
c = c(c1, c2)
# plot(density(c), xlab = 'residuals',main = '')

oldpar <- par(mfrow=c(1, 2))
plot(density(m$residuals), xlab = 'Residuals',main = '(a)', las=1)
plot(density(c), xlab = 'Residuals',main = '(b)', las=1)
par(oldpar)
# ====================== Figure 2 End ========================================



# ====================== Figure 3 Start ========================================
load('./real_data/real_data.RData')

results = results_mcp

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
# ====================== Figure 3 End ========================================



# ====================== Figure 4 Start ========================================
load('./real_data/real_data.RData')

group_index_est <- vector(length = n)
group_index_est[1:n %in% group_index[[1]]] <- 0
group_index_est[1:n %in% group_index[[2]]] <- 1


nine_dat$grouping <- group_index_est


scatterplotMatrix(~SBP1+DBP1+Age_cal+BMI_change | grouping, data=nine_dat,
                  regLine=FALSE,
                  smooth=FALSE,
                  diagonal=TRUE,
                  var.labels=c('SBP', 'DBP', 'Age', 'BMI change'),
                  cex.axis=2,
                  cex.labels=3)
# main="Scatter plot matrix for the data.",
# cex.main=2)
# ====================== Figure 4 End ========================================



# ====================== Table 4 Start ========================================
load('./real_data/real_data.RData')

# 计算系数估计的方差和 P-value
group_index_est <- vector(length = n)
group_index_est[1:n %in% group_index[[1]]] <- 0
group_index_est[1:n %in% group_index[[2]]] <- 1

H_qta = dummy(group_index_est)
H = kronecker(H_qta, diag(p))
# -------------------------------------------------------
M = solve(P_Z - t(WZ) %*% H %*% solve(t(H) %*% W_qta %*% H) %*% t(H) %*% WZ) %*% (t(2 * Z_a) - t(WZ) %*% H %*% solve(t(H) %*% W_qta %*% H) %*% t(2 * W_a %*% H))
N = solve(t(H) %*% W_qta %*% H - t(H) %*% WZ %*% solve(P_Z) %*% t(WZ) %*% H) %*% (t(2 * W_a %*% H) - t(H) %*% WZ %*% solve(P_Z) %*% t(2 * Z_a))

eta_oracle = M %*% y
alpha_oracle = N %*% y
## compute the sigma_hat
sigma_sqr = 1 / (n - q - K * p) * sum((y - Z_a %*% eta_oracle - W_a %*% H %*% alpha_oracle)^2)
sigma_hat = 1 / (n - q - K * p) * sum((y - Z_a %*% eta_star - W_a %*% H %*% c(alpha))^2)
## 
a11 = c(1, 0)
a12 = c(0, 1)
c12 = c(1, -1)

sd11 = sqrt(sigma_hat  * t(a11) %*% N %*% t(N) %*% a11)
sd12 = sqrt(sigma_hat  * t(a12) %*% N %*% t(N) %*% a12)


alpha/c(sd11,sd12)
pt(alpha/c(sd11,sd12), n-1)*2
(1 - pt(alpha/c(sd11,sd12), n-1))*2

t_z = eta_star/sqrt(diag(sigma_hat * M %*% t(M)))
t_z
(1 - pt(t_z[1:3], n-1))*2
pt(t_z[4:5], n-1)*2
(1 - pt(t_z[6], n-1))*2


pf((c12%*%alpha/(sigma_hat * t(c12) %*% N %*% t(N) %*% c12))/p,p, n-K*p-q)*2

pt(n-K*p-q)
pf(p, n-K*p-q)
# ====================== Table 4 End ========================================










# ===================== Appendix =========================================

## ============= Table 1 and Figure 3 Start ========================
load('./simulation/results/example1.RData')

EtaMse <- function(eta) {
  sqrt(sum((eta - eta_true)^2)) / sqrt(q)
}

EtaBias <- function(eta) {
  eta - eta_true
}


eta_or_mse <- apply(eta_or,1,EtaMse)
eta_scad_mse <- apply(eta_est_scad,1,EtaMse)
eta_mcp_mse <- apply(eta_est_mcp,1,EtaMse)
eta_ma_mse <- apply(eta_est_ma,1,EtaMse)



eta_or_bias <- apply(eta_or, 1, EtaBias)
eta_scad_bias <- apply(eta_est_scad, 1, EtaBias)
eta_mcp_bias <- apply(eta_est_mcp, 1, EtaBias)
eta_ma_bias <- apply(eta_est_ma, 1, EtaBias)


apply(eta_scad_bias, 1, mean)
apply(eta_scad_bias, 1, sd)
apply(eta_mcp_bias, 1, mean)
apply(eta_mcp_bias, 1, sd)
apply(eta_ma_bias, 1, mean)
apply(eta_ma_bias, 1, sd)


# boxplot(apply(eta_or_bias, 2, mean), apply(eta_scad_bias, 2, mean), apply(eta_mcp_bias, 2, mean), apply(eta_ma_bias, 2, mean))

oldpar <- par(mfrow = c(1, 3))
boxplot(eta_scad_bias[1, ], eta_mcp_bias[1, ], eta_ma_bias[1, ], names = c('SCAD', 'MCP','Ma'),
        ylab = expression(paste('Bias for ', hat(symbol(eta))[1])), ylim = c(-1, 1), las=1
)
boxplot(eta_scad_bias[2, ], eta_mcp_bias[2, ], eta_ma_bias[2, ], names = c('SCAD', 'MCP','Ma'),
        ylab = expression(paste('Bias for ', hat(symbol(eta))[2])), ylim = c(-1, 1), las=1
)
boxplot(eta_scad_bias[3, ], eta_mcp_bias[3, ], eta_ma_bias[3, ], names = c('SCAD', 'MCP', 'Ma'),
        ylab = expression(paste('Bias for ', hat(symbol(eta))[3])), ylim = c(-1, 1), las=1
)
par(oldpar)
title(main = 'Example 1')


load('./simulation/results/example2.RData')

eta_or_mse <- apply(eta_or,1,EtaMse)
eta_scad_mse <- apply(eta_est_scad,1,EtaMse)
eta_mcp_mse <- apply(eta_est_mcp,1,EtaMse)
eta_ma_mse <- apply(eta_est_ma,1,EtaMse)



eta_or_bias <- apply(eta_or, 1, EtaBias)
eta_scad_bias <- apply(eta_est_scad, 1, EtaBias)
eta_mcp_bias <- apply(eta_est_mcp, 1, EtaBias)
eta_ma_bias <- apply(eta_est_ma, 1, EtaBias)


apply(eta_scad_bias, 1, mean)
apply(eta_scad_bias, 1, sd)
apply(eta_mcp_bias, 1, mean)
apply(eta_mcp_bias, 1, sd)
apply(eta_ma_bias, 1, mean)
apply(eta_ma_bias, 1, sd)


# boxplot(apply(eta_or_bias, 2, mean), apply(eta_scad_bias, 2, mean), apply(eta_mcp_bias, 2, mean), apply(eta_ma_bias, 2, mean))

oldpar <- par(mfrow = c(1, 3))
boxplot(eta_scad_bias[1, ], eta_mcp_bias[1, ], eta_ma_bias[1, ], names = c('SCAD', 'MCP','Ma'),
        ylab = expression(paste('Bias for ', hat(symbol(eta))[1])), ylim = c(-1, 1), las=1
)
boxplot(eta_scad_bias[2, ], eta_mcp_bias[2, ], eta_ma_bias[2, ], names = c('SCAD', 'MCP','Ma'),
        ylab = expression(paste('Bias for ', hat(symbol(eta))[2])), ylim = c(-1, 1), las=1
)
boxplot(eta_scad_bias[3, ], eta_mcp_bias[3, ], eta_ma_bias[3, ], names = c('SCAD', 'MCP', 'Ma'),
        ylab = expression(paste('Bias for ', hat(symbol(eta))[3])), ylim = c(-1, 1), las=1
)
par(oldpar)
title(main = 'Example 2')
## ============= Table 1 and Figure 3 End ========================



## ============= Figure 1 Start ========================

### use example 1 for illustration, similar procedure for example 2 and 3
source("UpdateDelta.R")
source("iteration_ma.R")
source("iteration_procedures.R")

library(MASS)
library(snow)
library(fossil)
library(rlecuyer)

n = 100
q = 3
p = 1
rho_z = 0.3
sigma = 1
sigma_me_sqr = 0.25
corr_Z = matrix(rho_z, nrow = q-1, ncol = q-1)
diag(corr_Z) = 1
eta_true = matrix(c(1, 1, 1), nrow = q)
K_true = 2
alpha1_true = 2
alpha2_true = -2
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
Cn = log(n * p + q)

z = mvrnorm(n, rep(1, q-1), corr_Z)
Z = cbind(rep(1, n), z)

w = matrix(rnorm(n, 1, 1))
# w = cbind(rep(1, n), rnorm(n, 0, sigma))
# u = rnorm(n)


xi_z1 = rnorm(n, mean = 0, sd = sqrt(sigma_me_sqr))
xi_z2 = rnorm(n, mean = 0, sd = sqrt(sigma_me_sqr))
xi_z3 = rnorm(n, mean = 0, sd = sqrt(sigma_me_sqr))
xi_z4 = rnorm(n, mean = 0, sd = sqrt(sigma_me_sqr))
xi_x1 = rnorm(n, mean = 0, sd = sqrt(sigma_me_sqr))
xi_x2 = rnorm(n, mean = 0, sd = sqrt(sigma_me_sqr))



Z1 = cbind(Z[, 1], Z[, 2] + xi_z1, Z[, 3] + xi_z2)
Z2 = cbind(Z[, 1], Z[, 2] + xi_z3, Z[, 3] + xi_z4)

w1 = w + xi_x1
w2 = w + xi_x2

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

a = c(rep(0,n/2), rep(1, n/2))
b = sample(1:n, n, replace = FALSE)
group_index_true = a[b]
# group_index_true = matrix((Z[, 1]^2 + u < 1) * 1, nrow = n)
real_tag = rep(group_index_true, each = p)
beta_true = matrix(alpha1_true * (real_tag == 0) + alpha2_true * (real_tag == 1), nrow = n)
y = Z %*% eta_true + W %*% beta_true + rnorm(n, 0, 0.5)

# generate matrix Z_a
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

## Solution Path
lambda_min = 1/n + 0.02
# lambda_max = max(abs(t(cbind(Z_a, W_a)) %*% y) / n)
lambda_max = 1.5
k = 20
lambda_list = exp(seq(log(lambda_min), log(lambda_max), length.out = k))

beta_hat = matrix(nrow = n * p, ncol = k)
eta_hat = matrix(nrow = q, ncol = k)
delta_hat = matrix(nrow = choose(n, 2) * p, ncol = k)
K_hat = matrix(nrow = 1, ncol = k)
BIC = matrix(nrow = 1, ncol = k)

cl <- makeCluster(6)
clusterExport(cl, c('n', 'p', 'q', 'A'))
clusterExport(cl,c('L2Norm','S','delta_lasso', 'delta_SCAD','delta_MCP','UpdateDelta'))
set.seed(10)
clusterSetupRNG(cl, type = "RNGstream", seed=c(1, 22, 333, 444, 55, 6))
clusterExport(cl, c('Z_a', "W_a", "P_Z", "W_qta", "WZ", "V_ZW", "U_ZW", "Q_Z"))

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
## ============= Figure 1 End ========================


## ============= Figure 2 Start ========================
load('./simulation/results/example1.RData')

Alpha1Mse <- function(alpha){
  sqrt((alpha - alpha1_true)^2) / sqrt(p)
}

Alpha2Mse <- function(alpha){
  sqrt((alpha - alpha2_true)^2) / sqrt(p)
}


EtaMse <- function(eta) {
  sqrt(sum((eta - eta_true)^2)) / sqrt(q)
}

alpha1_or <- Alpha1Mse(alpha1_oracle)
alpha1_scad <- Alpha1Mse(na.omit(alpha1_est_scad))
alpha1_mcp <- Alpha1Mse(na.omit(alpha1_est_mcp))
alpha1_ma <- Alpha1Mse(alpha1_est_ma)


alpha2_or <- Alpha2Mse(alpha2_oracle)
alpha2_scad <- Alpha2Mse(na.omit(alpha2_est_scad))
alpha2_mcp <- Alpha2Mse(na.omit(alpha2_est_mcp))
alpha2_ma <- Alpha2Mse(alpha2_est_ma)


alpha1_bias_or = alpha1_oracle - alpha1_true
alpha1_bias_mcp = alpha1_est_mcp - alpha1_true
alpha1_bias_scad = alpha1_est_scad - alpha1_true
alpha1_bias_ma = alpha1_est_ma - alpha1_true


alpha2_bias_or = alpha2_oracle - alpha2_true
alpha2_bias_mcp = alpha2_est_mcp - alpha2_true
alpha2_bias_scad = alpha2_est_scad - alpha2_true
alpha2_bias_ma = alpha2_est_ma - alpha2_true



oldpar <- par(mfrow = c(2, 2))
boxplot(alpha1_bias_scad, alpha1_bias_mcp, alpha1_bias_ma, names = c('SCAD', 'MCP','Ma'),
        ylab = expression(paste('Bias for ', hat(symbol(a))[1])), ylim = c(-1, 1), las=1, cex.lab=1, cex.axis=1
)
boxplot(alpha2_bias_scad, alpha2_bias_mcp, alpha2_bias_ma, names = c('SCAD', 'MCP','Ma'),
        ylab = expression(paste('Bias for ', hat(symbol(a))[2])), ylim = c(-1, 1), las=1, cex.lab=1, cex.axis=1
)

boxplot(alpha1_scad, alpha1_mcp, alpha1_ma,  names = c('SCAD', 'MCP','Ma'),
        ylab = expression(paste('MSE for ', hat(symbol(a))[1])), ylim = c(0, 1), cex.lab=1, cex.axis=1
)
boxplot(alpha2_scad, alpha2_mcp, alpha2_ma, names = c('SCAD', 'MCP', 'Ma'),
        ylab = expression(paste('MSE for ', hat(symbol(a))[2])), ylim = c(0, 1), cex.lab=1,cex.axis=1
)

par(oldpar)

## ============= Figure 2 End ========================




## ============= Figure 4 and Table 2 Start ========================
load('./real_data/real_data.RData')

# results = results_mcp
results = results_scad
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

# 计算系数估计的方差和 P-value
group_index_est <- vector(length = n)
group_index_est[1:n %in% group_index[[1]]] <- 0
group_index_est[1:n %in% group_index[[2]]] <- 1

H_qta = dummy(group_index_est)
H = kronecker(H_qta, diag(p))
# -------------------------------------------------------
M = solve(P_Z - t(WZ) %*% H %*% solve(t(H) %*% W_qta %*% H) %*% t(H) %*% WZ) %*% (t(2 * Z_a) - t(WZ) %*% H %*% solve(t(H) %*% W_qta %*% H) %*% t(2 * W_a %*% H))
N = solve(t(H) %*% W_qta %*% H - t(H) %*% WZ %*% solve(P_Z) %*% t(WZ) %*% H) %*% (t(2 * W_a %*% H) - t(H) %*% WZ %*% solve(P_Z) %*% t(2 * Z_a))

eta_oracle = M %*% y
alpha_oracle = N %*% y
## compute the sigma_hat
sigma_sqr = 1 / (n - q - K * p) * sum((y - Z_a %*% eta_oracle - W_a %*% H %*% alpha_oracle)^2)
sigma_hat = 1 / (n - q - K * p) * sum((y - Z_a %*% eta_star - W_a %*% H %*% c(alpha))^2)
## 
a11 = c(1, 0)
a12 = c(0, 1)
c12 = c(1, -1)

sd11 = sqrt(sigma_hat  * t(a11) %*% N %*% t(N) %*% a11)
sd12 = sqrt(sigma_hat  * t(a12) %*% N %*% t(N) %*% a12)


alpha/c(sd11,sd12)
pt(alpha/c(sd11,sd12), n-1)*2
(1 - pt(alpha/c(sd11,sd12), n-1))*2

t_z = eta_star/sqrt(diag(sigma_hat * M %*% t(M)))
t_z
(1 - pt(t_z[1:3], n-1))*2
pt(t_z[4:5], n-1)*2
(1 - pt(t_z[6], n-1))*2


pf((c12%*%alpha/(sigma_hat * t(c12) %*% N %*% t(N) %*% c12))/p,p, n-K*p-q)*2

pt(n-K*p-q)
pf(p, n-K*p-q)
## ============= Figure 4 and Table 2 End ========================







