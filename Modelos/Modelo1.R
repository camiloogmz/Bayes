library(readxl)
library(sqldf)
library(dplyr)
library(data.table)
library(tictoc)
library(coda)
library(tidyr)




setwd("C:/Users/camil/OneDrive - Universidad Nacional de Colombia/Documentos/Camilo UNAL/Estadística Bayesiana/Caso de Estudio/Modelo 1")
datos <- read.delim("datos.txt", sep = ";")


datos <- datos %>% arrange(cole_cod_depto_ubicacion, cole_cod_mcpio_ubicacion)



# Funciones para muestrear los parámetros ----

sample_Xi <- function(sig2_k, n_jk, n_k, Kappa2_jk, theta_k, ybar_jk, n_j){
  
  v <- 1/(1/rep(sig2_k,n_k) + n_jk/Kappa2_jk)
  
  mn <- (rep(theta_k/sig2_k, n_k) + n_jk*ybar_jk/Kappa2_jk)*v
  
  Xi <- rnorm(n_j, mean = mn, sd = sqrt(v))
}

sample_Kappa2_jk <- function(nu_k, n_jk, n_k, Kappa2_k, s2_jk, ybar_jk, Xi, n_j){
  
  a <- nu_k + n_jk
  
  b <- nu_k*rep(Kappa2_k, n_k) + (n_jk-1)*s2_jk + n_jk*(ybar_jk - Xi)^2
  
  Kappa2_jk <- 1/rgamma(n_j, shape = a/2, rate = b/2)
}

sample_theta_k <- function(Xi, deptos, tau2, n_k, sig2_k, mu, m){
  
  T_k <- as.numeric(rowsum(Xi, deptos))
  
  v <- 1/(1/tau2 + n_k/sig2_k)
  
  mn <- (mu/tau2 + T_k/sig2_k)*v
  
  theta_k <- rnorm(m, mean = mn, sd = sqrt(v))
}

sample_sig2_k <- function(Xi, theta_k, deptos, nu_sig, n_k, sig2, m){
  
  T_sig2 <- as.numeric(rowsum((Xi - rep(theta_k, n_k))^2, group = deptos))
  
  a <- nu_sig + n_k
  
  b <- nu_sig*sig2 + T_sig2
  
  sig2_k <- 1/rgamma(m, shape = a/2, rate = b/2)
}

sample_Kappa2_k <- function(Kappa2_jk, deptos, alpha_kappa, beta_kappa, n_k, nu_k, m){
  
  T_Kappa <- as.numeric(rowsum(1/Kappa2_jk, group = deptos))
  
  a <- alpha_kappa + n_k*nu_k
  
  b <- beta_kappa + nu_k*T_Kappa
  
  Kappa2_k <- rgamma(m, shape = a/2, rate = b/2)
}

sample_mu <- function(sig2_mu, m, tau2, mu_mu, theta_k){
  
  v <- 1/(1/sig2_mu + m/tau2)
  
  mn <- (mu_mu/sig2_mu + m*mean(theta_k)/tau2)*v
  
  mu <- rnorm(1, mean = mn, sd = sqrt(v))
}

sample_tau2 <- function(nu_tau, m, sig2_tau, theta_k, mu){
  
  a <- nu_tau + m
  
  b <- nu_tau*sig2_tau + sum((theta_k - mu)^2)
  
  tau2 <- 1/rgamma(1, shape = a/2, rate = b/2)
}

sample_s2 <- function(a_sig, nu_sig, b_sig, m, sig2_k){
  
  a <- a_sig + m*nu_sig
  
  b <- b_sig + nu_sig*sum(1/sig2_k)
  
  s2 <- rgamma(1, shape = a/2, rate = b/2)
}


lpost_alpha  <- function(alpha_kappa, m, beta_kappa, Kappa2_k, a_alpha_k, b_alpha_k){
  (m*alpha_kappa/2)*log(beta_kappa/2) - m*lgamma(alpha_kappa/2) + (alpha_kappa/2)*sum(log(Kappa2_k))
  + (a_alpha_k - 1)*log(alpha_kappa) - b_alpha_k*alpha_kappa 
}

sample_alpha_kappa <- function(alpha_kappa, delta, m, beta_kappa, 
                               Kappa2_k, a_alpha_k, b_alpha_k, ac){
  
  # Propuesta
  
  eta_b <- log(alpha_kappa)
  eta_ast <- rnorm(1, mean = eta_b, sd = sqrt(delta))
  alpha_ast <- exp(eta_ast)
  
  # Tasa de aceptación
  
  r <- exp(lpost_alpha(alpha_ast, m, beta_kappa, Kappa2_k, a_alpha_k, b_alpha_k) 
           - lpost_alpha(alpha_kappa, m, beta_kappa, Kappa2_k, a_alpha_k, b_alpha_k)
           + log(alpha_ast) - log(alpha_kappa))
  
  # Actualización
  
  if(runif(1) < r){
    alpha_kappa <- alpha_ast
    ac <- ac + 1
  }
  
  return(list(alpha_kappa = alpha_kappa, ac = ac, delta = delta))
}


sample_beta_kappa <- function(a_beta_k, alpha_kappa, b_beta_k, Kappa2_k, m){
  
  a <- 2*a_beta_k + m*alpha_kappa
  
  b <- 2*b_beta_k + sum(Kappa2_k)
  
  beta_kappa <- rgamma(1, shape = a/2, rate = b/2)
}




# Datos para la cadena ----

## Estadísticas por departamento ----

est_dpto <- datos %>%
  group_by(cole_cod_depto_ubicacion, cole_depto_ubicacion) %>%
  summarise(
    ybar_k = mean(punt_global),
    s2_k = var(punt_global),
    n_k = n_distinct(cole_cod_mcpio_ubicacion),
    n_ijk = n()
  ) %>%
  ungroup()

m <- nrow(est_dpto)

## Estadísticas por municipio ----

est_mun <- datos %>%
  group_by(cole_cod_mcpio_ubicacion, cole_mcpio_ubicacion, cole_cod_depto_ubicacion) %>%
  summarise(
    ybar_jk = mean(punt_global),
    s2_jk = var(punt_global),
    n_jk = n(),
    .groups = "drop"
  ) 

n_j <- nrow(est_mun)




# Cadena ----


MCMC1 <- function(y, datos, est_mun, est_dpto, n_burn, n_skip, n_iter, delta){
  
  
  y <- datos$punt_global
  deptos <- factor(est_mun$cole_cod_depto_ubicacion)
  
  # Estadísticos Suficientes 
  
  ybar_jk <- est_mun$ybar_jk # Municipio
  s2_jk <- est_mun$s2_jk
  
  ybar_k <- est_dpto$ybar_k # Departamento
  s2_k <- est_dpto$s2_k
  var_ybar_jk <- with(est_mun,
                      tapply(ybar_jk,
                             cole_cod_depto_ubicacion,
                             FUN = var))
  var_ybar_jk[is.na(var_ybar_jk)] <- mean(var_ybar_jk, na.rm = TRUE)
  
  mean.s2_jk <- with(est_mun,
                     tapply(s2_jk,
                            cole_cod_depto_ubicacion,
                            FUN = mean))
  
  
  # Tamaños
  
  n_j <- nrow(est_mun) # Número de municipios
  m <- nrow(est_dpto) # Número de departamentos
  n_k <- est_dpto$n_k # Municipios por Departamento
  n_jk <- est_mun$n_jk # Número de estudiantes por municipio
  
  
  
  # Almacenamiento
  
  XI <- matrix(NA, nrow = n_iter, ncol = n_j)
  KAPPA2_JK  <- matrix(NA, nrow = n_iter, ncol = n_j)
  THETA_K <- matrix(NA, nrow = n_iter, ncol = m)
  SIG2_K <- matrix(NA, nrow = n_iter, ncol = m)
  KAPPA2_K <- matrix(NA, nrow = n_iter, ncol = m)
  MU <- matrix(NA, nrow = n_iter, ncol = 1)
  TAU2 <- matrix(NA, nrow = n_iter, ncol = 1)
  SIG2 <- matrix(NA, nrow = n_iter, ncol = 1)
  ALPHA_KAPPA <- matrix(NA, nrow = n_iter, ncol = 1)
  BETA_KAPPA <- matrix(NA, nrow = n_iter, ncol = 1)
  LL <- matrix(NA, nrow = n_iter, ncol = 1)
  
  # Hiperparámetros
  
  nu_k <- 1
  nu_sig <- 1
  
  mu_mu <- 250
  sig2_mu <- 50^2
  
  nu_tau <- 1
  sig2_tau <- 50^2
  
  a_sig <- 1
  b_sig <- 1/50^2
  
  a_alpha_k <- 1
  b_alpha_k <- 1
  
  a_beta_k <- 1
  b_beta_k <- 50^2
  
  
  # Valores Iniciales
  
  Xi <- ybar_jk
  Kappa2_jk <- s2_jk
  
  theta_k <- ybar_k
  sig2_k <- var_ybar_jk
  
  alpha_kappa <- 1
  beta_kappa <- 1/50^2
  
  Kappa2_k <- mean.s2_jk
  
  mu <- mean(theta_k)
  tau2 <- var(theta_k)
  
  sig2 <- 50^2
  
  
  # Número de iteraciones
  
  B <- n_burn + n_skip * n_iter
  
  ncat <- floor(B / 10)
  
  
  # Guardar las muestras
  
  for (b in 1:B) {
    
    # Muestrear Parámetros
    
    Xi <- sample_Xi(sig2_k, n_jk, n_k, Kappa2_jk, theta_k, ybar_jk, n_j)
    
    Kappa2_jk <- sample_Kappa2_jk(nu_k, n_jk, n_k, Kappa2_k, s2_jk, ybar_jk, Xi, n_j)
    
    theta_k <- sample_theta_k(Xi, deptos, tau2, n_k, sig2_k, mu, m)
    
    sig2_k <- sample_sig2_k(Xi, theta_k, deptos, nu_sig, n_k, sig2, m)
    
    Kappa2_k <- sample_Kappa2_k(Kappa2_jk, deptos, alpha_kappa, beta_kappa, n_k, nu_k, m)
    
    mu <- sample_mu(sig2_mu, m, tau2, mu_mu, theta_k)
    
    tau2 <- sample_tau2(nu_tau, m, sig2_tau, theta_k, mu)
    
    sig2 <- sample_s2(a_sig, nu_sig, b_sig, m, sig2_k)
    
    t <- sample_alpha_kappa(alpha_kappa, delta, m, beta_kappa, 
                            Kappa2_k, a_alpha_k, b_alpha_k, ac)
    
    alpha_kappa <- t$alpha_kappa
    
    beta_kappa <- sample_beta_kappa(a_beta_k, alpha_kappa, b_beta_k, Kappa2_k, m)
    
    
    ac <- t$ac
    
    
    ll <- sum(dnorm(x = y, mean = rep(Xi, n_jk), sd = sqrt(rep(Kappa2_jk, n_jk)), log = TRUE))
    
    # Almacenar
    
    if ((b > n_burn) & (b %% n_skip == 0)) {
      
      i <- (b - n_burn) / n_skip
      
      XI[i,] <- Xi
      KAPPA2_JK[i,]  <- Kappa2_jk
      THETA_K[i,] <- theta_k
      SIG2_K[i,] <- sig2_k
      KAPPA2_K[i,] <- Kappa2_k
      MU[i] <- mu
      TAU2[i] <- tau2
      SIG2[i] <- sig2
      ALPHA_KAPPA[i] <- alpha_kappa
      BETA_KAPPA[i] <- beta_kappa
      LL[i] <- ll
      
    }
    
    # Progreso
    
    if (b %% ncat == 0) {
      cat(sprintf("Iteración: %d, Tasa de aceptación: %.2f%%, delta: %.3f ",
                  b, (ac / b) * 100, delta))
    }
    
  }
  
  # Salida final
  
  return(list(XI = XI,
              KAPPA2_JK = KAPPA2_JK,
              THETA_K = THETA_K,
              SIG2_K = SIG2_K,
              KAPPA2_K = KAPPA2_K,
              MU = MU,
              TAU2 =TAU2,
              SIG2 = SIG2,
              ALPHA_KAPPA = ALPHA_KAPPA,
              BETA_KAPPA = BETA_KAPPA,
              LL = LL))
  
}







# Modelo 1 ----


ac <- 0
delta <- 5.5
n_burn = 1000; n_iter = 1000; n_skip = 10

set.seed(1)

tictoc::tic()

M1 <- MCMC1(y, datos, est_mun, est_dpto, n_burn, n_skip, n_iter, delta)

tictoc::toc()

M1 <- readRDS("M1_data.rds")

## save(M1, file = "M1.RData")


neff <- lapply(M1, effectiveSize)

lapply(neff, summary)

lapply(M1, colMeans)


EEMC_XI <- apply(X = M1$XI, MARGIN = 2, FUN = sd)/sqrt(neff$XI)
CVMC_XI <- 100*abs(EEMC_XI/colMeans(M1$XI))
round(summary(CVMC_CSI), 4)
hist(CVMC_XI)

# KAPPA2_jk 

EEMC_KAPPA2_jk  <- apply(X = M1$KAPPA2_JK, MARGIN = 2, FUN = sd)/sqrt(neff$KAPPA2_JK)
CVMC_KAPPA2_jk  <- 100*abs(EEMC_KAPPA2_jk/colMeans(M1$KAPPA2_JK))
round(summary(CVMC_KAPPA2_jk), 4)

# THETA_k

EEMC_THETA_k <- apply(X = M1$THETA_K, MARGIN = 2, FUN = sd)/sqrt(neff$THETA_K)
CVMC_THETA_k <- 100*abs(EEMC_THETA_k/colMeans(M1$THETA_K))
round(summary(CVMC_THETA_k), 4)
hist(CVMC_THETA_k)

# S2_k

EEMC_S2_k <- apply(X = M1$SIG2_K, MARGIN = 2, FUN = sd)/sqrt(neff$SIG2_K)
CVMC_S2_k <- 100*abs(EEMC_S2_k/colMeans(M1$SIG2_K))
round(summary(CVMC_S2_k), 4)
hist(CVMC_S2_k)

# KAPPA2_k 

EEMC_KAPPA2_k <- apply(X = M1$KAPPA2_K, MARGIN = 2, FUN = sd)/sqrt(neff$KAPPA2_K)
CVMC_KAPPA2_k <- 100*abs(EEMC_KAPPA2_k/colMeans(M1$KAPPA2_K))
round(summary(CVMC_KAPPA2_k), 4)
hist(CVMC_KAPPA2_k)


# MU

EEMC_MU <- sd(unlist(M1$MU))/sqrt(neff$MU)
CVMC_MU <- 100*abs(EEMC_MU/mean(unlist(M1$MU)))
CVMC_MU


# TAU2

EEMC_TAU2 <- sd(unlist(M1$TAU2))/sqrt(neff$TAU2)
CVMC_TAU2 <- 100*abs(EEMC_TAU2/mean(unlist(M1$TAU2)))
CVMC_TAU2

# S2

EEMC_S2 <- sd(unlist(M1$SIG2))/sqrt(neff$SIG2)
CVMC_S2 <- 100*abs(EEMC_S2/mean(unlist(M1$SIG2)))
CVMC_S2

# BETA_KAPPA

EEMC_BETA_KAPPA <- sd(unlist(M1$BETA_KAPPA))/sqrt(neff$BETA_KAPPA)
CVMC_BETA_KAPPA <- 100*abs(EEMC_BETA_KAPPA/mean(unlist(M1$BETA_KAPPA)))
CVMC_BETA_KAPPA


# ALPHA_KAPPA

EEMC_ALPHA_KAPPA <- sd(unlist(M1$ALPHA_KAPPA))/sqrt(neff$ALPHA_KAPPA)
CVMC_ALPHA_KAPPA <- 100*abs(EEMC_ALPHA_KAPPA/mean(unlist(M1$ALPHA_KAPPA)))
CVMC_ALPHA_KAPPA


