rm(list=ls())
# install.packages("nleqslv") # Décommenter si pas encore installé
library(nleqslv)
library(xtable)

#----------------------------------------------------------
#------------------ Fonctions Définies --------------------
#----------------------------------------------------------

# 1. Randomized Systematic Sampling
syspps = function(N_f,x){ 
  ##Population is first randomized!
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

# 2. Estimate theta with proposed method
hat_theta_f = function(xsa,xsb,wsb){
  col_xsa = colSums(xsa)
  theta_0 = rep(0,length(col_xsa))
  # Solution initiale via moindres carrés pondérés approximatifs ou solve simple
  # Note: Ajout d'une protection pour la matrice singulière si nécessaire
  theta_1 = solve(t(xsb) %*% (wsb*xsb), col_xsa-t(xsb) %*% wsb)
  
  while (abs(max(theta_1-theta_0)) > 10^(-8)){
    theta_0 = theta_1
    ps = 1/(1 + exp(-xsb %*% theta_0))
    # Newton-Raphson update
    theta_1 = theta_0 + solve(t(xsb) %*% (c(wsb*ps*(1-ps))*xsb), col_xsa - t(xsb) %*% c(wsb*ps))
  }
  return(theta_1)
}

# 3. Estimate theta with KH method
# J'ai ajouté 'theta_init' comme argument pour éviter les erreurs de portée de variable
hk_TT_f = function(xsa,xsb,wsb, theta_init){
  col_xsa = colSums(xsa)
  theta_0 = rep(0,length(col_xsa))
  theta_1 = theta_init
  
  while (abs(max(theta_1-theta_0)) > 10^(-8)){
    theta_0 = theta_1
    ps = 1/(1 + exp(-xsa %*% theta_0))
    # Note: Dans la méthode KH, on utilise xsa pour la mise à jour
    theta_1 = theta_0 + solve(t(xsa) %*% (c((1-ps)/ps)*xsa), t(xsa) %*% c(1/ps) - t(xsb) %*% c(wsb))
  }
  return(theta_1)
}

# 4. Bootstrap variance
boot_var_f = function(B_f,ysa,xsa_res,xsb_res,xsa_reg,xsb_reg,wsb,nsa,nsb){
  bdr2 = numeric(0)
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
    
    bdr2[i] = sum(bwsb*c(bxsb_reg %*% bbeta_ls))/sum(bwsb) + sum((1/bscore)*c(bysa-bxsa_reg %*% bbeta_ls))/sum(1/bscore)
  }
  return(var(bdr2))
}


#----------------------------------------------------------
#------------------ Paramètres et Génération --------------
#----------------------------------------------------------

set.seed(123)
B = 1000                    # bootstrap sample size (peut être réduit pour tester plus vite)
N = 20000                   # population total
iteration = 500             # simulation iteration                        
beta = c(2,1,1,1,1)         # True beta for "y=x beta + e"
theta = c(0.1,0.2,0.1,0.2)  # True theta for PS, no intercept       

# Covariate Generation
x1 = rbinom(N,1,0.5) 
x2 = runif(N,min=0,max=2)+0.3*x1                                     
x3 = rexp(N,1)+0.2*(x1+x2)                                    
x4 = rchisq(N,4)+0.1*(x1+x2+x3)

# Full Design matrix (Used for generation and TRUE models)
X_full = cbind(rep(1,N),x1,x2,x3,x4)

# Misspecified Design matrix (Without x4)
X_miss = cbind(rep(1,N),x1,x2,x3)

# Error term
e = rnorm(N,mean=0,sd=1)

rho = 0.3             # cor(Y,XB) target
mean_n_SA = 500       # sample size of S_A 
mean_n_SB = 1000      # sample size of S_B 

# Find sigma to match rho
sigma = sqrt((1/rho^2-1)*var(c(X_full%*%beta)))

# Generate Y
Y = X_full %*% beta + e * sigma
# Check cor(Y,XB)
# cor(Y,X_full %*% beta)
mu = mean(Y)

# Find intercept for Propensity Score
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

# True score (propensity)
score = exp(M + etap)/(1 + exp(M + etap))   

# Probability Sample Parameters (SB)
z = x3
c_val = (max(z) - 50*min(z))/49 # renamed 'c' to 'c_val' to avoid conflict with c() function
s = sum(z + c_val)
w = s/(z + c_val)/mean_n_SB
z1 = (z + c_val)/s


#----------------------------------------------------------
#------------------ Boucle des Scénarios ------------------
#----------------------------------------------------------

# Liste des scénarios à tester (Comme dans le Tableau 1 de l'article)
scenarios = c("TT", "FT", "TF", "FF")
results_summary = data.frame()

for (scen in scenarios) {
  
  cat("\nTraitement du scénario :", scen, "\n")
  
  # 1. Configuration des matrices selon le scénario
  # TT: Propensity (Res) Correct, Outcome (Reg) Correct
  # FT: Propensity (Res) Correct, Outcome (Reg) Faux (misspecified)
  # TF: Propensity (Res) Faux, Outcome (Reg) Correct
  # FF: Propensity (Res) Faux, Outcome (Reg) Faux
  
  if (scen == "TT") {
    X_res = X_full
    X_reg = X_full
  } else if (scen == "FT") {
    X_res = X_full
    X_reg = X_miss 
  } else if (scen == "TF") {
    X_res = X_miss
    X_reg = X_full
  } else if (scen == "FF") {
    X_res = X_miss
    X_reg = X_miss
  }
  
  Bias_scen = 0
  MSE_scen = 0
  est_var = 0 # Si vous voulez traquer la variance aussi
  
  # 2. Itération Monte Carlo
  for (g in 1:iteration){
    
    # --- Tirage des échantillons (Le mécanisme de sélection reste le VRAI) ---
    sam = rbinom(N,1,score)
    n_SA = sum(sam)
    index_SA = order(-sam)[1:n_SA] 
    Y_SA = Y[index_SA]
    
    # Selection des covariables "visibles" pour l'analyste selon le scénario
    X_SA_reg = X_reg[index_SA,]
    X_SA_res = X_res[index_SA,]
    
    index_SB = syspps(N,z1*mean_n_SB)
    n_SB = mean_n_SB
    
    # Selection des covariables "visibles" pour SB
    X_SB_res = X_res[index_SB,]
    X_SB_reg = X_reg[index_SB,]
    w_SB = w[index_SB]
    
    N_B = sum(w_SB)
    
    
    # --- Estimations ---
    
    # 1. Estimated theta and scores
    hat_theta = hat_theta_f(X_SA_res,X_SB_res,w_SB)
    hat_score_SA = 1/(1+exp(-X_SA_res%*%hat_theta))
    # hat_score_SB = 1/(1+exp(-X_SB_res%*%hat_theta)) # Pas strictement nécessaire pour les estimateurs ponctuels sauf var
    
    N_A = sum(1/hat_score_SA)
    
    # 2. Estimated LS beta (Outcome model)
    ls_beta = solve(t(X_SA_reg) %*% (X_SA_reg), t(X_SA_reg) %*% Y_SA)
    
    # 3. Estimated HK theta/beta
    # Calculable uniquement si les dimensions correspondent
    k1 = dim(X_SA_res)[2]
    k2 = dim(X_SA_reg)[2]
    
    if (k1 == k2){
      hk_theta = hk_TT_f(X_SA_res,X_SB_res,w_SB, hat_theta) # Pass hat_theta as init
      hk_score_SA = 1/(1+exp(-X_SA_res%*%hk_theta))
      hk_beta = solve(t(X_SA_reg) %*% (c(1/hk_score_SA-1) * (X_SA_reg)),
                      t(X_SA_reg) %*% c((1/hk_score_SA-1) * Y_SA))
      
      mu_dr_hk = sum(w_SB*c(X_SB_reg %*% hk_beta))/N + sum((1/hk_score_SA)*c(Y_SA-X_SA_reg %*% hk_beta))/N
    } else {
      mu_dr_hk = NA
    }
    
    
    # --- Calcul des Estimateurs Ponctuels ---
    mu_naive = mean(Y_SA)
    mu_ipw1 = sum(Y_SA/hat_score_SA)/N
    mu_ipw2 = sum(Y_SA/hat_score_SA)/N_A
    mu_reg = sum(w_SB*c(X_SB_reg%*%ls_beta))/N_B
    mu_dr1 = sum(w_SB*c(X_SB_reg %*% ls_beta))/N + sum((1/hat_score_SA)*c(Y_SA-X_SA_reg %*% ls_beta))/N
    mu_dr2 = sum(w_SB*c(X_SB_reg %*% ls_beta))/N_B + sum((1/hat_score_SA)*c(Y_SA-X_SA_reg %*% ls_beta))/N_A
    
    muhat = c(mu_naive, mu_ipw1, mu_ipw2, mu_reg, mu_dr1, mu_dr2, mu_dr_hk)
    
    # Accumulation
    Bias_scen = Bias_scen + (muhat - mu)/iteration
    MSE_scen = MSE_scen + (muhat - mu)^2/iteration
    
    if(g %% 50 == 0) cat(".")
  }
  
  # Stockage des résultats pour le scénario courant
  # Ordre: Naive, IPW1, IPW2, REG, DR1, DR2, DR_HK
  row_res = data.frame(
    Scenario = rep(scen, 7),
    Estimator = c("Naive", "IPW1", "IPW2", "REG", "DR1", "DR2", "DR_HK"),
    RB_percent = Bias_scen/mu * 100,
    MSE = MSE_scen
  )
  results_summary = rbind(results_summary, row_res)
}

#----------------------------------------------------------
#------------------ Affichage Final -----------------------
#----------------------------------------------------------

cat("\n\n################ RÉSULTATS FINAUX ################\n")
print(xtable(results_summary, digits = 4), include.rownames=FALSE)

# Pour une vue plus lisible dans la console R directement:
print(results_summary, digits=4)

