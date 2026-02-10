rm(list=ls())
library(nleqslv)
library(xtable)

#----------------------------------------------------------
#------------------self-defined functions------------------
#----------------------------------------------------------

# randomized systematic sampling
syspps = function(N_f,x){ 
  U = sample(N_f,N_f)
  xx = x[U]
  cx = cumsum(xx)
  r = runif(1)
  s = numeric()
  for(i in 1:N_f){
    if(cx[i] >= r){
      s = c(s,U[i])
      r = r + 1
    }
  }
  return(s)
}

# Estimate theta with proposed method (Maximum pseudo-likelihood)
hat_theta_f = function(xsa,xsb,wsb){
  col_xsa = colSums(xsa)
  theta_0 = rep(0,length(col_xsa))
  theta_1 = solve(t(xsb) %*% (wsb*xsb), col_xsa-t(xsb) %*% wsb)
  while (abs(max(theta_1-theta_0)) > 10^(-8)){
    theta_0 = theta_1
    ps = 1/(1 + exp(-xsb %*% theta_0))
    theta_1 = theta_0 + solve(t(xsb) %*% (c(wsb*ps*(1-ps))*xsb), col_xsa - t(xsb) %*% c(wsb*ps))
  }
  return(theta_1)
}

# Estimate theta with incorrect pooling method (mu_C1 and mu_C2)
hat_theta_pooled = function(xsa,xsb){
  # Create pooled indicator: R_tilde = 1 for SA, 0 for SB
  nsa = nrow(xsa)
  nsb = nrow(xsb)
  
  X_pooled = rbind(xsa, xsb)
  R_pooled = c(rep(1, nsa), rep(0, nsb))
  
  # Unweighted logistic regression
  theta_0 = rep(0, ncol(xsa))
  theta_1 = theta_0 + 0.1
  
  while (abs(max(theta_1-theta_0)) > 10^(-8)){
    theta_0 = theta_1
    ps = 1/(1 + exp(-X_pooled %*% theta_0))
    
    # Score function
    score = t(X_pooled) %*% (R_pooled - ps)
    
    # Hessian
    hess = -t(X_pooled) %*% (c(ps*(1-ps))*X_pooled)
    
    theta_1 = theta_0 - solve(hess, score)
  }
  
  return(theta_1)
}

# Estimate theta with KH method
# Estimate theta and beta with KH method (CORRECT IMPLEMENTATION)
hk_TT_f = function(xsa, xsb, wsb, ysa, xsa_reg){
  k1 = ncol(xsa)      # dimension of propensity score model
  k2 = ncol(xsa_reg)  # dimension of outcome regression model
  
  if (k1 != k2){
    return(list(theta=rep(NA,k1), beta=rep(NA,k2)))
  }
  
  nsa = nrow(xsa)
  nsb = nrow(xsb)
  
  # Initial values
  theta_0 = rep(0, k1)
  beta_0 = solve(t(xsa_reg) %*% xsa_reg, t(xsa_reg) %*% ysa)
  
  max_iter = 100
  tol = 10^(-8)
  
  for(iter in 1:max_iter){
    ps = 1/(1 + exp(-xsa %*% theta_0))
    
    # Equation (16): sum_{i in SA} (1/pi(theta) - 1) * (yi - m(xi, beta)) * xi = 0
    eq16 = t(xsa) %*% c((1/ps - 1) * (ysa - xsa_reg %*% beta_0))
    
    # Equation (17): sum_{i in SA} (1/pi(theta)) * xi = sum_{i in SB} wsb * xi
    # Rearranged: sum_{i in SA} (1/pi(theta)) * xi - sum_{i in SB} wsb * xi = 0
    eq17 = t(xsa_reg) %*% c(1/ps) - t(xsb) %*% wsb
    
    # Combined score function
    score = c(eq16, eq17)
    
    # Check convergence
    if(max(abs(score)) < tol){
      break
    }
    
    # Jacobian matrix (approximation for Newton-Raphson)
    # d(eq16)/d(theta)
    J11 = -t(xsa) %*% (c((1/ps - 1) * ps * (1-ps)) * xsa) - 
      t(xsa) %*% (c((1/ps^2) * ps * (1-ps) * (ysa - xsa_reg %*% beta_0)) * xsa)
    
    # d(eq16)/d(beta)
    J12 = -t(xsa) %*% (c((1/ps - 1)) * xsa_reg)
    
    # d(eq17)/d(theta)
    J21 = -t(xsa_reg) %*% (c((1/ps^2) * ps * (1-ps)) * xsa)
    
    # d(eq17)/d(beta)
    J22 = matrix(0, k2, k2)
    
    # Jacobian matrix
    J = rbind(cbind(J11, J12), cbind(J21, J22))
    
    # Newton-Raphson update
    update = solve(J, -score)
    
    theta_0 = theta_0 + update[1:k1]
    beta_0 = beta_0 + update[(k1+1):(k1+k2)]
    
    if(iter == max_iter){
      warning("KH method did not converge")
    }
  }
  
  return(list(theta=theta_0, beta=beta_0))
}

# Bootstrap variance
boot_var_f = function(B_f,ysa,xsa_res,xsb_res,xsa_reg,xsb_reg,wsb,nsa,nsb){
  bdr2 = numeric(B_f)
  for (i in 1:B_f){
    bsa_index = sample(nsa,nsa,replace=TRUE)
    bxsa_reg = xsa_reg[bsa_index,]
    bxsa_res = xsa_res[bsa_index,]
    bysa = ysa[bsa_index]
    
    bbeta_ls = solve(t(bxsa_reg) %*% (bxsa_reg), t(bxsa_reg) %*% bysa)
    
    bsb_index = sample(nsb,nsb,replace=TRUE)
    bxsb_reg = xsb_reg[bsb_index,]
    bxsb_res = xsb_res[bsb_index,]
    bwsb = wsb[bsb_index]
    
    btheta = hat_theta_f(bxsa_res,bxsb_res,bwsb)
    bscore = 1/(1 + exp(-bxsa_res %*% btheta))
    
    bdr2[i] = sum(bwsb*c(bxsb_reg %*% bbeta_ls))/sum(bwsb) + 
      sum((1/bscore)*c(bysa-bxsa_reg %*% bbeta_ls))/sum(1/bscore)
  }
  
  return(var(bdr2))
}

################## Parameter settings ##############

set.seed(123)
B = 5000                    # Bootstrap sample size
N = 20000                   # Population total
iteration = 500             # Simulation iterations (Table uses 5000)
beta = c(2,1,1,1,1)         # True beta for y=X*beta + e
theta = c(0.1,0.2,0.1,0.2)  # True theta for PS (no intercept)

# Generate covariates
x1 = rbinom(N,1,0.5) 
x2 = runif(N,min=0,max=2) + 0.3*x1                                                   
x3 = rexp(N,1) + 0.2*(x1+x2)                                                    
x4 = rchisq(N,4) + 0.1*(x1+x2+x3)

# Design matrices (these will be modified for different scenarios)
X_full = cbind(rep(1,N), x1, x2, x3, x4)

# Simulation parameters
rho = 0.3               # Correlation: try 0.3, 0.5, 0.8
mean_n_SA = 500         # Sample size of S_A: try 500, 1000
mean_n_SB = 1000        # Sample size of S_B: try 500, 1000

# Model scenarios: "TT", "FT", "TF", "FF"
scenario = "TT"

# Set working models based on scenario
if(scenario == "TT"){
  # Both correct
  X_reg = X_full
  X_res = X_full
} else if(scenario == "FT"){
  # Outcome model wrong, PS model correct
  X_reg = cbind(rep(1,N), x1, x2, x3)  # x4 omitted
  X_res = X_full
} else if(scenario == "TF"){
  # Outcome model correct, PS model wrong
  X_reg = X_full
  X_res = cbind(rep(1,N), x1, x2, x3)  # x4 omitted
} else { # "FF"
  # Both wrong
  X_reg = cbind(rep(1,N), x1, x2, x3)
  X_res = cbind(rep(1,N), x1, x2, x3)
}

# Find sigma for desired correlation
sigma = sqrt((1/rho^2-1)*var(c(X_full%*%beta)))

# Generate y
e = rnorm(N, mean=0, sd=1)
Y = X_full %*% beta + e * sigma

# Check correlation
cat("Actual correlation:", cor(Y, X_full %*% beta), "\n")
mu = mean(Y)

# Find intercept for propensity scores
etap = as.vector(X_full[,-1] %*% theta)
dif = 1
L = -10
R = -1

while(dif > 0){
  M = (L + R)/2
  ps = exp(M + etap)/(1 + exp(M + etap))
  if(sum(ps) <= mean_n_SA) L = M
  if(sum(ps) >= mean_n_SA) R = M
  if(abs(sum(ps) - mean_n_SA) <= 0.5) dif = 0
}

# True propensity scores
score = exp(M + etap)/(1 + exp(M + etap))   
cat("Min score:", min(score), "Max score:", max(score), "Sum:", sum(score), "\n")

# Set up PPS sampling weights (max/min = 50)
z = x3
c_const = (max(z) - 50*min(z))/49
s = sum(z + c_const)
w = s/(z + c_const)/mean_n_SB
z1 = (z + c_const)/s

#----------------------------------------------------------
#------------------ Simulation starts here ----------------
#----------------------------------------------------------

# Storage for all estimators
estimator_names = c("mu_A", "mu_C1", "mu_C2", "mu_IPW1", "mu_IPW2", 
                    "mu_REG", "mu_DR1", "mu_DR2")
n_est = length(estimator_names)

Bias = rep(0, n_est)
MSE = rep(0, n_est)

# For variance estimation (indices for IPW1, IPW2, DR2_plug, DR2_boot, KH)
var_est_names = c("vIPW1", "vIPW2", "vDR2_PLUG", "vDR2_BOOT", "vKH")
n_var = length(var_est_names)
est_var = rep(0, n_var)
c_rate = rep(0, n_var)

for (g in 1:iteration){
  
  # Draw non-probability sample (Poisson sampling)
  sam = rbinom(N, 1, score)
  n_SA = sum(sam)
  index_SA = which(sam == 1)
  
  Y_SA = Y[index_SA]
  X_SA_reg = X_reg[index_SA,]
  X_SA_res = X_res[index_SA,]
  
  # Draw probability sample with systematic PPS
  index_SB = syspps(N, z1*mean_n_SB)
  n_SB = length(index_SB)
  
  X_SB_res = X_res[index_SB,]
  X_SB_reg = X_reg[index_SB,]
  w_SB = w[index_SB]
  
  N_B = sum(w_SB)
  
  #-----------------------------------------------------
  #---------------- Propensity Score Estimation --------
  #-----------------------------------------------------
  
  # Proposed method (correct weighting)
  hat_theta = hat_theta_f(X_SA_res, X_SB_res, w_SB)
  hat_score_SA = 1/(1 + exp(-X_SA_res %*% hat_theta))
  hat_score_SB = 1/(1 + exp(-X_SB_res %*% hat_theta))
  N_A = sum(1/hat_score_SA)
  
  # Pooled method (incorrect - for mu_C1 and mu_C2)
  hat_theta_pool = hat_theta_pooled(X_SA_res, X_SB_res)
  hat_score_SA_pool = 1/(1 + exp(-X_SA_res %*% hat_theta_pool))
  N_A_pool = sum(1/hat_score_SA_pool)
  
  # Regression coefficient estimation
  ls_beta = solve(t(X_SA_reg) %*% X_SA_reg, t(X_SA_reg) %*% Y_SA)
  
  # Kim-Haziza method (only if dimensions match)
  k1 = ncol(X_SA_res)
  k2 = ncol(X_SA_reg)
  
  if (k1 == k2){
    hk_result = hk_TT_f(X_SA_res, X_SB_res, w_SB, Y_SA, X_SA_reg)
    hk_theta = hk_result$theta
    hk_beta = hk_result$beta
    hk_score_SA = 1/(1 + exp(-X_SA_res %*% hk_theta))
  } else {
    hk_theta = rep(NA, k1)
    hk_beta = rep(NA, k2)
    hk_score_SA = NA
  }
  
  #-----------------------------------------------------
  #---------------- Point Estimates --------------------
  #-----------------------------------------------------
  
  mu_A = mean(Y_SA)
  mu_C1 = sum(Y_SA/hat_score_SA_pool)/N
  mu_C2 = sum(Y_SA/hat_score_SA_pool)/N_A_pool
  mu_IPW1 = sum(Y_SA/hat_score_SA)/N
  mu_IPW2 = sum(Y_SA/hat_score_SA)/N_A
  mu_REG = sum(w_SB * c(X_SB_reg %*% ls_beta))/N_B
  mu_DR1 = sum(w_SB * c(X_SB_reg %*% ls_beta))/N + 
    sum((1/hat_score_SA) * c(Y_SA - X_SA_reg %*% ls_beta))/N
  mu_DR2 = sum(w_SB * c(X_SB_reg %*% ls_beta))/N_B + 
    sum((1/hat_score_SA) * c(Y_SA - X_SA_reg %*% ls_beta))/N_A
  
  muhat = c(mu_A, mu_C1, mu_C2, mu_IPW1, mu_IPW2, mu_REG, mu_DR1, mu_DR2)
  
  Bias = Bias + (muhat - mu)/iteration
  MSE = MSE + (muhat - mu)^2/iteration
  
  #-----------------------------------------------------
  #---------------- Variance Estimation ----------------
  #-----------------------------------------------------
  
  # Plug-in variance estimators
  a = t(X_SB_res) %*% (c(w_SB*hat_score_SB*(1-hat_score_SB))*X_SB_res)
  
  # For IPW1
  hat_b1 = solve(a, t(X_SA_res) %*% c((1/hat_score_SA-1)*Y_SA))
  tt1 = c(X_SB_res %*% hat_b1) * c(hat_score_SB)
  st1 = sum(w_SB*tt1)
  est_vd1 = 1/N^2 * sum((w_SB*tt1 - st1/n_SB)^2)
  est_vr1 = (1/N^2)*sum((1-hat_score_SA)*(Y_SA/hat_score_SA - X_SA_res %*% hat_b1)^2)
  est_var_ipw1 = est_vd1 + est_vr1
  
  # For IPW2
  hat_b2 = solve(a, t(X_SA_res) %*% c((1/hat_score_SA-1)*(Y_SA-mu_IPW2)))
  tt2 = c(X_SB_res %*% hat_b2) * c(hat_score_SB)
  st2 = sum(w_SB*tt2)
  est_vd2 = 1/N^2 * sum((w_SB*tt2 - st2/n_SB)^2)
  est_vr2 = (1/N^2)*sum((1-hat_score_SA)*((Y_SA-mu_IPW2)/hat_score_SA - X_SA_res %*% hat_b2)^2)
  est_var_ipw2 = est_vd2 + est_vr2
  
  # For DR2 (plug-in)
  hat_d0 = mu_DR2 - mu_REG
  hat_b4 = solve(a, t(X_SA_res) %*% c((1/hat_score_SA-1)*(Y_SA-X_SA_reg %*% ls_beta-hat_d0)))
  tt4 = c(X_SB_reg %*% ls_beta) + c(hat_score_SB)*c(X_SB_res %*% hat_b4) - mu_REG
  st4 = sum(w_SB*tt4)
  est_vd4 = 1/N^2 * sum((w_SB*tt4 - st4/n_SB)^2)
  est_vr4 = (1/N^2)*sum((1-hat_score_SA)*((Y_SA-X_SA_reg %*% ls_beta-hat_d0)/hat_score_SA - X_SA_res %*% hat_b4)^2)
  est_var_dr2_plug = est_vd4 + est_vr4
  
  # Bootstrap variance for DR2
  boot_var_dr2 = boot_var_f(B, Y_SA, X_SA_res, X_SB_res, X_SA_reg, X_SB_reg, w_SB, n_SA, n_SB)
  
  # Kim-Haziza DR variance
  if(k1 == k2 && !any(is.na(hk_score_SA))){
    tt_hk = c(X_SB_reg %*% hk_beta)
    st_hk = sum(w_SB*tt_hk)
    est_vdhk = 1/N^2 * sum((w_SB*tt_hk - st_hk/n_SB)^2)
    est_vrhk = (1/N^2)*sum((1-hk_score_SA)*((Y_SA-X_SA_reg %*% hk_beta)/hk_score_SA)^2)
    hat_sigma2 = mean((Y_SA - X_SA_reg %*% hk_beta)^2)
    est_var_hk = est_vdhk + est_vrhk - (1/N^2)*(sum(1/hk_score_SA) - N_B) * hat_sigma2
  } else {
    est_var_hk = NA
  }
  
  # Accumulate variance estimates
  var_estimates = c(est_var_ipw1, est_var_ipw2, est_var_dr2_plug, boot_var_dr2, est_var_hk)
  est_var = est_var + var_estimates/iteration
  
  # Confidence intervals and coverage
  sd_ipw1 = sqrt(est_var_ipw1)
  sd_ipw2 = sqrt(est_var_ipw2)
  sd_dr2_plug = sqrt(est_var_dr2_plug)
  sd_dr2_boot = sqrt(boot_var_dr2)
  
  CI_ipw1_c = (mu >= mu_IPW1 - 1.96*sd_ipw1) * (mu <= mu_IPW1 + 1.96*sd_ipw1)
  CI_ipw2_c = (mu >= mu_IPW2 - 1.96*sd_ipw2) * (mu <= mu_IPW2 + 1.96*sd_ipw2)
  CI_dr2_plug_c = (mu >= mu_DR2 - 1.96*sd_dr2_plug) * (mu <= mu_DR2 + 1.96*sd_dr2_plug)
  CI_dr2_boot_c = (mu >= mu_DR2 - 1.96*sd_dr2_boot) * (mu <= mu_DR2 + 1.96*sd_dr2_boot)
  
  if(!is.na(est_var_hk)){
    sd_hk = sqrt(est_var_hk)
    mu_KH = sum(w_SB*c(X_SB_reg %*% hk_beta))/N + sum((1/hk_score_SA)*c(Y_SA-X_SA_reg %*% hk_beta))/N
    CI_hk_c = (mu >= mu_KH - 1.96*sd_hk) * (mu <= mu_KH + 1.96*sd_hk)
  } else {
    CI_hk_c = NA
  }
  
  coverage = c(CI_ipw1_c, CI_ipw2_c, CI_dr2_plug_c, CI_dr2_boot_c, CI_hk_c)
  c_rate = c_rate + coverage/iteration
  
  if(floor(g/100) == g/100) cat("Iteration:", g, "\n")
}

#################  Final Results  ###############

MC_variance = MSE - Bias^2
RB = 100*Bias/mu

# Variance estimator relative bias
var_indices = c(4, 5, 8, 8, NA)  # Indices for IPW1, IPW2, DR2, DR2, KH
MC_var_subset = MC_variance[var_indices[!is.na(var_indices)]]
RB_var = (est_var[!is.na(est_var)] - MC_var_subset)/MC_var_subset * 100

################### Output Tables ###################

# Table 1: Point estimators
output1 = data.frame(
  Estimator = estimator_names,
  RB = round(RB, 2),
  MSE = round(MSE, 2)
)

table_title1 = paste0("Scenario=", scenario, ", nA=", mean_n_SA, ", nB=", mean_n_SB, 
                      ", rho=", rho)
cat("\n", table_title1, "\n")
print(xtable(output1, caption = table_title1, digits = 2),
      caption.placement = 'top', include.rownames = FALSE)

# Table 2: Variance estimators
output2 = data.frame(
  Estimator = var_est_names,
  RB_var = round(c(RB_var, ifelse(is.na(est_var[5]), NA, 
                                  (est_var[5]-MC_variance[8])/MC_variance[8]*100)), 2),
  CP = round(c_rate*100, 2)
)

table_title2 = paste0("Variance estimators - Scenario=", scenario, 
                      ", nA=", mean_n_SA, ", nB=", mean_n_SB, ", rho=", rho)
cat("\n", table_title2, "\n")
print(xtable(output2, caption = table_title2, digits = 2),
      caption.placement = 'top', include.rownames = FALSE)

# Summary output
cat("\n=== Summary ===\n")
cat("True mean mu:", round(mu, 3), "\n")
cat("Number of iterations:", iteration, "\n")
cat("Scenario:", scenario, "\n")