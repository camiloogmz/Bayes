library(sqldf)
library(dplyr)
library(tictoc)
library(coda)
library(mvtnorm)
library(Matrix)


setwd("C:/Users/camil/OneDrive - Universidad Nacional de Colombia/Documentos/Camilo UNAL/Estadística Bayesiana/Caso de Estudio/M2")

load("exp_covariables.RData")

datos <- dataM2 %>% arrange(cole_cod_depto_ubicacion, cole_cod_mcpio_ubicacion)

rm(dataM2)


fit <- lm(punt_global ~ edu_madre + computador + internet 
          + libros + estrato + etnia
          + nbi.m + doc_est.m + riesgo.m
          + PIB.d + pob_rural.d + riesgo.d, data = datos)


# Funciones para muestrear los parámetros ----

sample_Beta <- function(s2_beta, s2_E, s2_M, s2_D, p_E, p_M, p_D, Kappa2_jk, n_jk,
                        X, y, Beta0){
  
  d0   <- c(s2_beta, 
            rep(s2_E, p_E), 
            rep(s2_M, p_M), 
            rep(s2_D, p_D))
  i0   <- 1 / d0              # Diag(Sigma_0^{-1})  
  
  Sigma.inv <- 1/rep(Kappa2_jk, n_jk)
  
  X.S <- sweep(X, 1, Sigma.inv, FUN = "*")    # X^t Sigma^{-1}
  Xt.Sinv.X <- crossprod(X, X.S)     # X^t Sigma^{-1} X
  Xt.Sinv.y <- crossprod(X, Sigma.inv * y)     # X^t Sigma^{-1} y
  
  P <- diag(i0) + Xt.Sinv.X
  M <- i0 * Beta0 + Xt.Sinv.y
  
  V <- chol2inv(chol(P))
  mn <- V %*% M
  
  Beta <- rmvnorm(1, mean = as.vector(mn), sigma = V)
  
  return(as.vector(Beta))
}


sample_Kappa2_jk <- function(nu_k, n_jk, Kappa2_k, n_k, Xi, y, mcpio, n_j){
  
  a <- nu_k + n_jk
  
  b <- nu_k*rep(Kappa2_k, n_k) + as.numeric(rowsum((y - Xi)^2, mcpio))
  
  kappa2_jk <- 1/rgamma(n_j, shape = a/2, rate = b/2)
}


sample_s2_beta <- function(nu_beta, gamma_beta, Beta, mu_beta){
  
  a <- nu_beta + 1
  
  b <- nu_beta*gamma_beta + sum((Beta - mu_beta)^2)
  
  s2_beta <- 1/rgamma(1, shape = a/2, rate = b/2)
}


sample_s2_E <- function(Beta, nu_E, p_E, gamma_E, mu_E){
  
  Beta_E <- Beta[2:14]
  
  a <- nu_E + p_E
  
  b <- nu_E*gamma_E + as.numeric(crossprod(Beta_E - mu_E))
  
  s2_E <- 1/rgamma(1, shape = a/2, rate = b/2)
}


sample_s2_M <- function(Beta, nu_M, p_M, gamma_M, mu_M){
  
  Beta_M <- Beta[15:17]
  
  a <- nu_M + p_M

  b <- nu_M*gamma_M + as.numeric(crossprod(Beta_M - mu_M))
  
  s2_M <- 1/rgamma(1, shape = a/2, rate = b/2)
}



sample_s2_D <- function(Beta, nu_D, p_D, gamma_D, mu_D){
  
  Beta_D <- Beta[18:20]
  
  a <- nu_D + p_D 
  
  b <- nu_D*gamma_D + as.numeric(crossprod(Beta_D - mu_D))
  
  s2_D <- 1/rgamma(1, shape = a/2, rate = b/2)
}


sample_Kappa2_k <- function(alpha_kappa, n_k, nu_k, beta_kappa, Kappa2_jk, deptos, m){
  
  a <- alpha_kappa + n_k*nu_k
  
  b <- beta_kappa + nu_k*as.numeric(rowsum(1/Kappa2_jk, deptos))
  
  Kappa2_k <- rgamma(m, shape = a/2, rate = b/2)
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
  
  b <- 2*b_beta_k + sum(1/Kappa2_k)
  
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


MCMC2 <- function(y, fit, datos, est_mun, est_dpto, n_burn, n_skip, n_iter, delta){
  
  
  y <- datos$punt_global
  X <- model.matrix(fit)
  deptos <- factor(est_mun$cole_cod_depto_ubicacion)
  mcpio <- factor(datos$cole_cod_mcpio_ubicacion)
  
  Beta_OLS <- as.vector(fit$coefficients)
  s2_OLS <- sum(fit$residuals^2)/fit$df.residual 
  s_OLS <- sqrt(s2_OLS)
  
  
  # Estadísticos Suficientes 
  
  ybar_jk <- est_mun$ybar_jk # Municipio
  s2_jk <- est_mun$s2_jk
  
  ybar_k <- est_dpto$ybar_k # Departamento
  s2_k <- est_dpto$s2_k
  
  mean.s2_jk <- with(est_mun,
                     tapply(s2_jk,
                            cole_cod_depto_ubicacion,
                            FUN = mean))

  # Tamaños
  
  n_j <- nrow(est_mun) # Número de municipios
  m <- nrow(est_dpto) # Número de departamentos
  n_k <- est_dpto$n_k # Municipios por Departamento
  n_jk <- est_mun$n_jk # Número de estudiantes por municipio
  
  p_E <- ncol(X[,2:14]) # Parámetros por estudiante
  p_M <- ncol(X[,15:17]) # Parámetros por municipio
  p_D <- ncol(X[,18:20]) # Parámetros por departamento
  p <- fit$rank
  
  
  
  
  # Almacenamiento
  
  BETA <- matrix(NA, nrow = n_iter, ncol = p)
  KAPPA2_JK  <- matrix(NA, nrow = n_iter, ncol = n_j)
  S2_beta <- matrix(NA, nrow = n_iter, ncol = 1)
  S2_E <- matrix(NA, nrow = n_iter, ncol = 1)
  S2_M <- matrix(NA, nrow = n_iter, ncol = 1)
  S2_D <- matrix(NA, nrow = n_iter, ncol = 1)
  KAPPA2_K <- matrix(NA, nrow = n_iter, ncol = m)
  ALPHA_KAPPA <- matrix(NA, nrow = n_iter, ncol = 1)
  BETA_KAPPA <- matrix(NA, nrow = n_iter, ncol = 1)
  LL <- matrix(NA, nrow = n_iter, ncol = 1)
  
  # Hiperparámetros
  
  mu_beta <- mu_E <- mu_M <- mu_D <- 0
  
  Beta0 <- c(mu_beta, 
             rep(mu_E, p_E),
             rep(mu_M, p_M),
             rep(mu_D, p_D))
  
  nu_beta <- nu_E <- nu_M <- nu_D <- nu_k <- 1
  
  gamma_beta <- gamma_E <- gamma_M <- gamma_D <- 100
  
  a_alpha_k <- 1
  b_alpha_k <- 1
  
  a_beta_k <- 1
  b_beta_k <- 1/50^2
  
  
  # Valores Iniciales
  
  Beta <- Beta_OLS
  
  Kappa2_jk <- s2_jk
  
  alpha_kappa <- 1
  beta_kappa <- 1
  
  s2_beta <- s2_E <- s2_M <- s2_D <- s2_OLS
  
  Kappa2_k <- mean.s2_jk
  
  
  # Número de iteraciones
  
  B <- n_burn + n_skip * n_iter
  
  ncat <- floor(B / 10)
  
  
  # Guardar las muestras
  
  for (b in 1:B) {
    
    # Muestrear Parámetros
    
    Beta <- sample_Beta(s2_beta, s2_E, s2_M, s2_D, p_E, p_M, p_D, Kappa2_jk, n_jk,
                        X, y, Beta0)
    Xi <- as.vector(X %*% Beta)
    
    Kappa2_jk <- sample_Kappa2_jk(nu_k, n_jk, Kappa2_k, n_k, Xi, y, mcpio, n_j)
    
    s2_beta <- sample_s2_beta(nu_beta, gamma_beta, Beta, mu_beta)
    
    s2_E <- sample_s2_E(Beta, nu_E, p_E, gamma_E, mu_E)
    
    s2_M <- sample_s2_M(Beta, nu_M, p_M, gamma_M, mu_M)
    
    s2_D <- sample_s2_D(Beta, nu_D, p_D, gamma_D, mu_D)
    
    Kappa2_k <- sample_Kappa2_k(alpha_kappa, n_k, nu_k, beta_kappa, Kappa2_jk, deptos, m)
    
    t <- sample_alpha_kappa(alpha_kappa, delta, m, beta_kappa, 
                            Kappa2_k, a_alpha_k, b_alpha_k, ac)
    
    alpha_kappa <- t$alpha_kappa
    
    beta_kappa <- sample_beta_kappa(a_beta_k, alpha_kappa, b_beta_k, Kappa2_k, m)
    
    
    ac <- t$ac
    
    
    ll <- sum(dnorm(x = y, mean = Xi, sd = sqrt(rep(Kappa2_jk, n_jk)), log = TRUE))
    
    # Almacenar
    
    if ((b > n_burn) & (b %% n_skip == 0)) {
      
      i <- (b - n_burn) / n_skip
      
      BETA[i,] <- Beta
      KAPPA2_JK[i,]  <- Kappa2_jk
      S2_beta[i] <- s2_beta
      S2_E[i] <- s2_E
      S2_M[i] <- s2_M
      S2_D[i] <- s2_D
      KAPPA2_K[i,] <- Kappa2_k
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
  
  return(list(BETA = BETA,
              KAPPA2_JK = KAPPA2_JK,
              S2_beta = S2_beta,
              S2_E = S2_E,
              S2_M = S2_M,
              S2_D = S2_D,
              KAPPA2_K = KAPPA2_K,
              ALPHA_KAPPA = ALPHA_KAPPA,
              BETA_KAPPA = BETA_KAPPA,
              LL = LL))
  
}





# Modelo 2 ----


ac <- 0
delta <- 5.5
n_burn = 5000; n_iter = 5000; n_skip = 10


set.seed(1)

tictoc::tic()

M2 <- MCMC2(y, fit, datos, est_mun, est_dpto, n_burn, n_skip, n_iter, delta)

tictoc::toc()


setwd("C:/Users/camil/OneDrive - Universidad Nacional de Colombia/Documentos/Camilo UNAL/Estadística Bayesiana/Caso de Estudio/Modelo 2")

save(M2, file = "M2.Rdata")
load("M2_2.RData")


neff <- lapply(M2, effectiveSize)

lapply(neff, summary)

lapply(M2, colMeans)




EEMC_KAPPA2_jk  <- apply(X = M2$KAPPA2_JK, MARGIN = 2, FUN = sd)/sqrt(neff$KAPPA2_JK)
CVMC_KAPPA2_jk  <- 100*abs(EEMC_KAPPA2_jk/colMeans(M2$KAPPA2_JK))
round(summary(CVMC_KAPPA2_jk), 4)

EEMC_KAPPA2_k <- apply(X = M2$KAPPA2_K, MARGIN = 2, FUN = sd)/sqrt(neff$KAPPA2_K)
CVMC_KAPPA2_k <- 100*abs(EEMC_KAPPA2_k/colMeans(M2$KAPPA2_K))
round(summary(CVMC_KAPPA2_k), 4)
hist(CVMC_KAPPA2_k)

EEMC_BETA_KAPPA <- sd(unlist(M2$BETA_KAPPA))/sqrt(neff$BETA_KAPPA)
CVMC_BETA_KAPPA <- 100*abs(EEMC_BETA_KAPPA/mean(unlist(M2$BETA_KAPPA)))
CVMC_BETA_KAPPA

