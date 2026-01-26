rm(list=ls())
install.packages("nleqslv")
library(nleqslv)
library(xtable)
#----------------------------------------------------------
#------------------self-defined functions------------------
#----------------------------------------------------------

# randomized systematic sampling
# Echantillon de probabiliste (S_B)

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


# Estimate theta with proposed method
# Estimation avec score de propension (estimation de theta du modèle de regression logistique)

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



# Estimate theta with KH method

hk_TT_f = function(xsa,xsb,wsb){
  col_xsa = colSums(xsa)
  theta_0 = rep(0,length(col_xsa))
  theta_1 = hat_theta
  while (abs(max(theta_1-theta_0)) > 10^(-8)){
    theta_0 = theta_1
    ps = 1/(1 + exp(-xsa %*% theta_0))
    theta_1 = theta_0 + solve(t(xsa) %*% (c((1-ps)/ps)*xsa), t(xsa) %*% c(1/ps) - t(xsb) %*% c(wsb))
  }
  return(theta_1)
}




#bootstrap variance

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



##################parameter settings##############



set.seed(123)
B = 1000                            #bootstrap sample size
N = 20000                           #population total
iteration = 500                     #simulation iteration                        
beta = c(2,1,1,1,1)                 #True beta for "y=x beta + e"
theta = c(0.1,0.2,0.1,0.2)          #True theta for PS, no intercept       
#Covariate
x1 = rbinom(N,1,0.5) 
x2 = runif(N,min=0,max=2)+0.3*x1                                                   
x3 = rexp(N,1)+0.2*(x1+x2)                                                    
x4 = rchisq(N,4)+0.1*(x1+x2+x3)

#Design matrix
X = cbind(rep(1,N),x1,x2,x3,x4)
#working matrices
X_reg = cbind(rep(1,N),x1,x2,x3,x4)
X_res = cbind(rep(1,N),x1,x2,x3,x4)
e = rnorm(N,mean=0,sd=1)

rho = 0.3              #cor(Y,XB)
mean_n_SA = 500        #sample size of S_A 
mean_n_SB = 1000       #sample size of S_B 

#find sigma
sigma = sqrt((1/rho^2-1)*var(c(X%*%beta)))

#Genreate y
Y = X %*% beta + e * sigma
#Check cor(Y,XB)
cor(Y,X %*% beta)
mu = mean(Y)

#find intercept
etap = as.vector(X[,-1] %*% theta)
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


#true score
score = exp(M + etap)/(1 + exp(M + etap))   
(min(score))
(max(score))
(sum(score))

#set max(z1)/min(z1)=50
#find c

z = x3
c = (max(z) - 50*min(z))/49
s = sum(z + c)
w = s/(z + c)/mean_n_SB
z1 = (z + c)/s

#----------------------------------------------------------
#------------------iteration starts here-------------------
#----------------------------------------------------------
Bias = 0
MSE = 0
est_var = 0
c_rate = 0
for (g in 1:iteration){
  
  #draw non-probability sample 
  sam = rbinom(N,1,score)
  n_SA = sum(sam)
  index_SA = order(-sam)[1:n_SA] 
  Y_SA = Y[index_SA]
  X_SA_reg = X_reg[index_SA,]
  X_SA_res = X_res[index_SA,]
  
  
  #draw probability sampe with sys method
  index_SB = syspps(N,z1*mean_n_SB)
  n_SB = mean_n_SB
  X_SB_res = X_res[index_SB,]
  X_SB_reg = X_reg[index_SB,]
  w_SB = w[index_SB]
  
  
  #check sum
  N_B = sum(w_SB)
  
  
  #estimated theta and scores
  hat_theta = hat_theta_f(X_SA_res,X_SB_res,w_SB)
  hat_score_SA = 1/(1+exp(-X_SA_res%*%hat_theta))
  hat_score_SB = 1/(1+exp(-X_SB_res%*%hat_theta))
  
  N_A = sum(1/hat_score_SA)
  
  #estimated ls beta 
  ls_beta = solve(t(X_SA_reg) %*% (X_SA_reg), t(X_SA_reg) %*% Y_SA)
  
  #estimated hk theta/beta
  
  #attain hk estimators if two equations have the same dim
  #set NA ow
  k1 = dim(X_SA_res)[2]
  k2 = dim(X_SA_reg)[2]
  
  if (k1 == k2){
    hk_theta=hk_TT_f(X_SA_res,X_SB_res,w_SB)
    hk_score_SA = 1/(1+exp(-X_SA_res%*%hk_theta))
    hk_beta = solve(t(X_SA_reg) %*% (c(1/hk_score_SA-1) * (X_SA_reg)),
                    t(X_SA_reg) %*% c((1/hk_score_SA-1) * Y_SA))
  }else{
    hk_theta = rep(NA,k1)
    hk_beta = rep(NA,k2)
    hk_score_SA = 1/(1 + exp(-X_SA_res %*% hk_theta))
  }
  
  
  #-----------------------------------------------------
  #--------------------point estimates------------------
  #-----------------------------------------------------
  
  mu_naive = mean(Y_SA)
  mu_ipw1 = sum(Y_SA/hat_score_SA)/N
  mu_ipw2 = sum(Y_SA/hat_score_SA)/N_A
  mu_reg = sum(w_SB*c(X_SB_reg%*%ls_beta))/N_B
  mu_dr1 = sum(w_SB*c(X_SB_reg %*% ls_beta))/N + sum((1/hat_score_SA)*c(Y_SA-X_SA_reg %*% ls_beta))/N
  mu_dr2 = sum(w_SB*c(X_SB_reg %*% ls_beta))/N_B + sum((1/hat_score_SA)*c(Y_SA-X_SA_reg %*% ls_beta))/N_A
  mu_dr_hk = sum(w_SB*c(X_SB_reg %*% hk_beta))/N + sum((1/hk_score_SA)*c(Y_SA-X_SA_reg %*% hk_beta))/N
  
  muhat = c(mu_naive,mu_ipw1,mu_ipw2,mu_reg,mu_dr1,mu_dr2,mu_dr_hk)
  Bias = Bias + (muhat - mu)/iteration
  MSE = MSE + (muhat - mu)^2/iteration
  
  
  
  #-----------------------------------------------------
  #-----------------------Variance-----------------------
  #-----------------------------------------------------
  a = t(X_SB_res) %*% (c(w_SB*hat_score_SB*(1-hat_score_SB))*X_SB_res)
  hat_b1 = solve(a,t(X_SA_res) %*% c((1/hat_score_SA-1)*Y_SA))
  hat_b2 = solve(a,t(X_SA_res) %*% c((1/hat_score_SA-1)*(Y_SA-mu_ipw2)))
  hat_b3 = solve(a,t(X_SA_res) %*% c((1/hat_score_SA-1)*(Y_SA-X_SA_reg %*% ls_beta)))
  hat_d0 = mu_dr2 - mu_reg
  hat_b4 = solve(a,t(X_SA_res) %*% c((1/hat_score_SA-1)*(Y_SA-X_SA_reg %*% ls_beta-hat_d0)))
  
  tt = c(X_SB_res %*% hat_b1)*c(hat_score_SB)
  st = sum(w_SB*tt)
  est_vd1 = 1/N^2 * sum((w_SB*tt-st/n_SB)^2)
  
  tt = c(X_SB_res %*% hat_b2)*c(hat_score_SB)
  st = sum(w_SB*tt)
  est_vd2 = 1/N^2 * sum((w_SB*tt-st/n_SB)^2)
  
  tt = c(X_SB_reg %*% ls_beta) + c(hat_score_SB)*c(X_SB_res %*% hat_b3)
  st = sum(w_SB*tt)
  est_vd3 = 1/N^2 * sum((w_SB*tt-st/n_SB)^2)
  
  
  tt = c(X_SB_reg %*% ls_beta) + c(hat_score_SB)*c(X_SB_res %*% hat_b4)- mu_reg
  st = sum(w_SB*tt)
  est_vd4 = 1/N^2 * sum((w_SB*tt-st/n_SB)^2)
  
  tt = c(X_SB_reg %*% hk_beta)
  st = sum(w_SB*tt)
  est_vdhk = 1/N^2 * sum((w_SB*tt-st/n_SB)^2)
  
  
  est_vr1 = (1/N^2)*sum((1-hat_score_SA)*(Y_SA/hat_score_SA-X_SA_res %*% hat_b1)^2)
  est_vr2 = (1/N^2)*sum((1-hat_score_SA)*((Y_SA-mu_ipw2)/hat_score_SA-X_SA_res %*% hat_b2)^2)
  est_vr3 = (1/N^2)*sum((1-hat_score_SA)*((Y_SA-X_SA_reg %*% ls_beta)/hat_score_SA-X_SA_res %*% hat_b3)^2)
  est_vr4 = (1/N^2)*sum((1-hat_score_SA)*((Y_SA-X_SA_reg %*% ls_beta-hat_d0)/hat_score_SA-X_SA_res %*% hat_b4)^2)
  est_vrhk = (1/N^2)*sum((1-hk_score_SA)*((Y_SA-X_SA_reg %*% hk_beta)/hk_score_SA)^2)
  
  est_var_ipw1 = est_vd1 + est_vr1
  est_var_ipw2 = est_vd2 + est_vr2
  est_var_dr1 = est_vd3 + est_vr3
  est_var_dr2 = est_vd4 + est_vr4
  hat_sigma2 = mean((Y_SA-X_SA_reg %*% hk_beta)^2)
  est_var_hk = est_vdhk + est_vrhk - (1/N^2)*(sum(1/hk_score_SA)-N_B) * hat_sigma2
  
  #bootstrap variance estimator.
  boot_var_dr2 = boot_var_f(B, Y_SA, X_SA_res, X_SB_res, X_SA_reg, X_SB_reg, w_SB, n_SA, n_SB)
  
  #estimated variance
  est_var = est_var + c(est_var_ipw1, est_var_ipw2, est_var_dr1, est_var_dr2, boot_var_dr2, est_var_hk)/iteration
  
  #sd
  sd_ipw1 = sqrt(est_var_ipw1)
  sd_ipw2 = sqrt(est_var_ipw2)
  sd_dr1 = sqrt(est_var_dr1)
  sd_dr2 = sqrt(est_var_dr2)
  sd_boot_dr2 = sqrt(boot_var_dr2)
  sd_dr_hk = sqrt(est_var_hk)
  
  
  #CI
  CI_ipw1 = c(mu_ipw1 - 1.96*sd_ipw1, mu_ipw1 + 1.96*sd_ipw1)
  CI_ipw2 = c(mu_ipw2 - 1.96*sd_ipw2, mu_ipw2 + 1.96*sd_ipw2)
  CI_dr1 = c(mu_dr1 - 1.96*sd_dr1, mu_dr1 + 1.96*sd_dr1)
  CI_dr2 = c(mu_dr2 - 1.96*sd_dr2, mu_dr2 + 1.96*sd_dr2)
  CI_boot_dr2 = c(mu_dr2 - 1.96*sd_boot_dr2, mu_dr2 + 1.96*sd_boot_dr2)
  CI_dr_hk = c(mu_dr_hk - 1.96*sd_dr_hk, mu_dr_hk + 1.96*sd_dr_hk)
  
  
  #cp
  
  CI_ipw1_c = (mu >= CI_ipw1[1])*(mu <= CI_ipw1[2])
  CI_ipw2_c = (mu >= CI_ipw2[1])*(mu <= CI_ipw2[2])
  CI_dr1_c = (mu >= CI_dr1[1])*(mu <= CI_dr1[2])
  CI_dr2_c = (mu >= CI_dr2[1])*(mu <= CI_dr2[2])
  CI_boot_dr2_c = (mu >= CI_boot_dr2[1])*(mu <= CI_boot_dr2[2])
  CI_dr_hk_c = (mu >= CI_dr_hk[1])*(mu <= CI_dr_hk[2])
  
  
  
  c_rate = c_rate + c(CI_ipw1_c, CI_ipw2_c, CI_dr1_c, CI_dr2_c, CI_boot_dr2_c, CI_dr_hk_c)/iteration
  
  
  if(floor(g/10) == g/10) print(g)
  
}





#################final result###############

MC_variance = MSE - (Bias)^2
RB = 100*Bias/mu
est_var = est_var

muhat = c(mu_naive, mu_ipw1 ,mu_ipw2, mu_reg, mu_dr1, mu_dr2, mu_dr_hk)

RB_var = (est_var - MC_variance[c(2,3,5,6,6,7)])/(MC_variance[c(2,3,5,6,6,7)])*100


###################LATEX##############
output1 = data.frame("Estimator"=c("mu_naive","mu_ipw1","mu_ipw2","mu_reg","mu_dr1","mu_dr2","mu_dr_hk"),
                     "RB" = RB,"MSE" = MSE)

table_title = paste0("na=",mean_n_SA,",","nb=",mean_n_SB,", ","cor(Y,XB)=",rho*100)
print(xtable(output1, caption = table_title,digits = 2),
      caption.placement = 'top',include.colnames = T)

output2 = data.frame("Estimator"=c("ipw1","ipw2","dr1","dr2","boot_dr2","dr_hk"), 
                     "RB" = RB_var, "CR" = c_rate*100)

table_title2 = paste0("na=",mean_n_SA,",","nb=",mean_n_SB,", ","cor(Y,XB)=",rho*100)
print(xtable(output2, caption = table_title2,digits = 2),
      caption.placement = 'top',include.colnames = T)

