# Chargement du package nécessaire pour le tirage PPS
if(!require(sampling)) install.packages("sampling")
library(sampling)

# --- 1. PARAMÈTRES DE LA SIMULATION ---
set.seed(42)
N <- 20000
nA <- 500
nB <- 1000
B <- 1000 # L'article utilise 5000, réduit ici pour la vitesse
rhos <- c(0.3, 0.5, 0.8)

# --- 2. GÉNÉRATION DE LA POPULATION ---
# Variables auxiliaires (Section 5, p. 8) [cite: 824, 825]
z1 <- rbinom(N, 1, 0.5)
z2 <- runif(N, 0, 2)
z3 <- rexp(N, 1)
z4 <- rchisq(N, 4)

x1 <- z1
x2 <- z2 + 0.3 * x1
x3 <- z3 + 0.2 * (x1 + x2)
x4 <- z4 + 0.1 * (x1 + x2 + x3)

# Calcul du prédicteur linéaire pour Y
lp_y <- 2 + x1 + x2 + x3 + x4
var_lp <- var(lp_y)

# Fonction pour calculer sigma selon rho 
get_sigma <- function(rho) sqrt(var_lp * (1/rho^2 - 1))

# --- 3. MÉCANISMES DE SÉLECTION ---
# Probabilités d'inclusion pour S_B (PPS sur z = c + x3) [cite: 834, 835]
# c est choisi pour que max(z)/min(z) = 50
z_pps <- x3
c_val <- (max(z_pps) - 50 * min(z_pps)) / 49
z_pps <- z_pps + c_val
pi_B <- inclusionprobabilities(z_pps, nB)
dB <- 1/pi_B

# Scores de propension pour S_A (Modèle logistique q) [cite: 831, 832]
# Trouver theta0 pour atteindre nA en moyenne
f_theta0 <- function(t0) sum(1 / (1 + exp(-(t0 + 0.1*x1 + 0.2*x2 + 0.1*x3 + 0.2*x4)))) - nA
theta0 <- uniroot(f_theta0, c(-20, 20))$root
pi_A_true <- 1 / (1 + exp(-(theta0 + 0.1*x1 + 0.2*x2 + 0.1*x3 + 0.2*x4)))

# --- 4. FONCTION D'ESTIMATION DU SCORE (Pseudo-Log-Vraisemblance) ---
# Résolution de l'équation du score (Eq 6) [cite: 668, 695]
estimate_theta <- function(dataA, dataB, formula) {
  X_A <- model.matrix(formula, data = dataA)
  X_B <- model.matrix(formula, data = dataB)
  weights_B <- dataB$d
  
  score_eq <- function(theta) {
    pi_est <- 1 / (1 + exp(-as.matrix(X_B) %*% theta))
    colSums(X_A) - colSums(weights_B * as.vector(pi_est) * X_B)
  }
  
  # Utilisation de nleqslv ou équivalent pour résoudre
  res <- try(rootSolve::multiroot(score_eq, start = rep(0, ncol(X_A)))$root, silent=T)
  if(inherits(res, "try-error")) return(rep(NA, ncol(X_A)))
  return(res)
}

# --- 5. BOUCLE DE SIMULATION ---
run_sim <- function(rho_val) {
  sigma <- get_sigma(rho_val)
  y <- lp_y + rnorm(N, 0, sigma)
  mu_y_true <- mean(y)
  
  sim_results <- replicate(B, {
    # Tirages [cite: 833, 834]
    S_A_idx <- which(rbinom(N, 1, pi_A_true) == 1)
    S_B_idx <- which(UPsystematic(pi_B) == 1)
    
    # Données
    dfA <- data.frame(y = y[S_A_idx], x1 = x1[S_A_idx], x2 = x2[S_A_idx], x3 = x3[S_A_idx], x4 = x4[S_A_idx])
    dfB <- data.frame(x1 = x1[S_B_idx], x2 = x2[S_B_idx], x3 = x3[S_B_idx], x4 = x4[S_B_idx], d = dB[S_B_idx])
    
    # Scénarios (TT, FT, TF, FF) [cite: 836-840]
    # Modèles corrects : y ~ x1+x2+x3+x4 | logit(pi) ~ x1+x2+x3+x4
    # Modèles faux : omission de x4
    
    results <- list()
    for(scenario in c("TT", "FT", "TF", "FF")) {
      f_pi <- if(scenario %in% c("TT", "FT")) ~x1+x2+x3+x4 else ~x1+x2+x3
      f_reg <- if(scenario %in% c("TT", "TF")) y~x1+x2+x3+x4 else y~x1+x2+x3
      
      # Estimation de pi_A
      th <- estimate_theta(dfA, dfB, f_pi)
      X_A_pi <- model.matrix(f_pi, data = dfA)
      pi_A_est <- 1 / (1 + exp(-as.matrix(X_A_pi) %*% th))
      
      # Estimation de m(x, beta)
      mod_y <- lm(f_reg, data = dfA)
      y_hat_B <- predict(mod_y, newdata = dfB)
      y_hat_A <- predict(mod_y, newdata = dfA)
      
      # Estimateurs [cite: 711, 732, 740]
      N_hat_A <- sum(1/pi_A_est)
      N_hat_B <- sum(dfB$d)
      
      mu_IPW2 <- sum(dfA$y / pi_A_est) / N_hat_A
      mu_REG <- sum(dfB$d * y_hat_B) / N_hat_B
      mu_DR2 <- (sum((dfA$y - y_hat_A) / pi_A_est) / N_hat_A) + (sum(dfB$d * y_hat_B) / N_hat_B)
      
      results[[scenario]] <- c(IPW2 = mu_IPW2, REG = mu_REG, DR2 = mu_DR2)
    }
    unlist(results)
  })
  
  # Calcul du %RB et MSE
  means <- rowMeans(sim_results, na.rm=T)
  rb <- (means - mu_y_true) / mu_y_true * 100
  mse <- rowMeans((sim_results - mu_y_true)^2, na.rm=T)
  
  return(list(RB = rb, MSE = mse))
}

# --- 6. EXÉCUTION ET AFFICHAGE ---
final_table <- lapply(rhos, run_sim)
# (Le formatage final en tableau peut être fait avec 'data.frame')