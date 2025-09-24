library(dplyr)
library(coda)

setwd("C:/Users/camil/OneDrive - Universidad Nacional de Colombia/Documentos/Camilo UNAL/Estadística Bayesiana/Caso de Estudio/Modelo 1")

datos <- read.delim("datos.txt", sep = ";")


datos <- datos %>% arrange(cole_cod_depto_ubicacion, cole_cod_mcpio_ubicacion)

# Modelo 1

M1 <- readRDS("M1.rds")

# Modelo 2

setwd("C:/Users/camil/OneDrive - Universidad Nacional de Colombia/Documentos/Camilo UNAL/Estadística Bayesiana/Caso de Estudio/Modelo 2")

load("M2.RData")

load("exp_covariables.RData")

fit <- lm(punt_global ~ edu_madre + computador + internet 
          + libros + estrato + etnia
          + nbi.m + doc_est.m + riesgo.m
          + PIB.d + pob_rural.d + riesgo.d, data = dataM2)

coef(fit)

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

n_jk <- est_mun$n_jk

# Resolución de Puntos ----

setwd("C:/Users/camil/OneDrive - Universidad Nacional de Colombia/Documentos/Camilo UNAL/Estadística Bayesiana/Caso de Estudio/Solución Puntos")

## 7. LL ----

### Por separado ----

pdf("cadena1.pdf", width = 7, height = 5) 

par(mar = c(2.75,2.75,1.5,0.5), mgp = c(1.7,0.7,0))
plot(M1$LL,   type = "p", pch = ".", cex = 2, col = "purple", 
     xlab = "Iteracion", ylab = "Log-verosimilitud", 
     main = expression(italic("Cadena del Modelo 1")))
abline(h = mean(M1$LL), col = "purple3", lwd = 3)

dev.off()

pdf("cadena2.pdf", width = 7, height = 5) 

par(mar = c(2.75,2.75,1.5,0.5), mgp = c(1.7,0.7,0))
plot(M1$LL,   type = "p", pch = ".", cex = 2, col = "dodgerblue", 
     xlab = "Iteracion", ylab = "Log-verosimilitud", 
     main = expression(italic("Cadena del Modelo 2")))
abline(h = mean(M1$LL), col = "dodgerblue2", lwd = 3)

dev.off()

### Conjunto ----

pdf("cadena_conjunta.pdf", width = 7, height = 5) 

yrange <- range(M1$LL, M2$LL)

par(mar = c(3.2, 3.2, 2, 0.8), mgp = c(2, 0.7, 0))


plot(M1$LL, type = "p", pch = ".", cex = 1.7, col = "purple",
     xlab = "Iteración", ylab = "Log-verosimilitud",
     main = expression(italic("Cadenas del Modelo 1 y 2")),
     ylim = yrange)


points(M2$LL, type = "p", pch = ".", cex = 1.7, col = "dodgerblue")

legend(x = 4000, y = -1243000,
       legend = c("Modelo 1", "Modelo 2"),
       col = c("purple", "dodgerblue"),
       pch = 20, pt.cex = 1.5, bty = "n", text.col = "black")

dev.off()

## 8. Tamaños efectivos de muestra ----

### Modelo 1 ----

neffM1 <- lapply(M1, effectiveSize)

lapply(neffM1, summary)

### Modelo 2 ----

neffM2 <- lapply(M2, effectiveSize)

lapply(neffM2, summary)


## 9. CVMC ----

### Modelo 1 ----

EEMC_XI <- apply(X = M1$XI, MARGIN = 2, FUN = sd)/sqrt(neffM1$XI)
round(summary(EEMC_XI), 4)
CVMC_XI <- 100*abs(EEMC_XI/colMeans(M1$XI))
round(summary(CVMC_XI), 4)
hist(CVMC_XI)

# KAPPA2_jk 

EEMC_KAPPA2_jk  <- apply(X = M1$KAPPA2_JK, MARGIN = 2, FUN = sd)/sqrt(neffM1$KAPPA2_JK)
round(summary(EEMC_KAPPA2_jk), 4)
CVMC_KAPPA2_jk  <- 100*abs(EEMC_KAPPA2_jk/colMeans(M1$KAPPA2_JK))
round(summary(CVMC_KAPPA2_jk), 4)

# THETA_k

EEMC_THETA_k <- apply(X = M1$THETA_K, MARGIN = 2, FUN = sd)/sqrt(neffM1$THETA_K)
round(summary(EEMC_THETA_k), 4)
CVMC_THETA_k <- 100*abs(EEMC_THETA_k/colMeans(M1$THETA_K))
round(summary(CVMC_THETA_k), 4)


# S2_k

EEMC_S2_k <- apply(X = M1$SIG2_K, MARGIN = 2, FUN = sd)/sqrt(neffM1$SIG2_K)
round(summary(EEMC_S2_k), 4)
CVMC_S2_k <- 100*abs(EEMC_S2_k/colMeans(M1$SIG2_K))
round(summary(CVMC_S2_k), 4)
hist(CVMC_S2_k)

# KAPPA2_k 

EEMC_KAPPA2_k <- apply(X = M1$KAPPA2_K, MARGIN = 2, FUN = sd)/sqrt(neffM1$KAPPA2_K)
round(summary(EEMC_KAPPA2_k), 4)
CVMC_KAPPA2_k <- 100*abs(EEMC_KAPPA2_k/colMeans(M1$KAPPA2_K))
round(summary(CVMC_KAPPA2_k), 4)
hist(CVMC_KAPPA2_k)


# MU

EEMC_MU <- sd(unlist(M1$MU))/sqrt(neffM1$MU)
CVMC_MU <- 100*abs(EEMC_MU/mean(unlist(M1$MU)))
CVMC_MU


# TAU2

EEMC_TAU2 <- sd(unlist(M1$TAU2))/sqrt(neffM1$TAU2)
CVMC_TAU2 <- 100*abs(EEMC_TAU2/mean(unlist(M1$TAU2)))
CVMC_TAU2

# S2

EEMC_S2 <- sd(unlist(M1$SIG2))/sqrt(neffM1$SIG2)
CVMC_S2 <- 100*abs(EEMC_S2/mean(unlist(M1$SIG2)))
CVMC_S2

# BETA_KAPPA

EEMC_BETA_KAPPA <- sd(unlist(M1$BETA_KAPPA))/sqrt(neffM1$BETA_KAPPA)
CVMC_BETA_KAPPA <- 100*abs(EEMC_BETA_KAPPA/mean(unlist(M1$BETA_KAPPA)))
CVMC_BETA_KAPPA


# ALPHA_KAPPA

EEMC_ALPHA_KAPPA <- sd(unlist(M1$ALPHA_KAPPA))/sqrt(neffM1$ALPHA_KAPPA)
CVMC_ALPHA_KAPPA <- 100*abs(EEMC_ALPHA_KAPPA/mean(unlist(M1$ALPHA_KAPPA)))
CVMC_ALPHA_KAPPA



### Modelo 2 ----

CVMC <- function(param, modelo, neff) {

  param_samples <- modelo[[param]]
  param_neff <- neff[[param]]
  
  eemc <- apply(param_samples, 2, sd) / sqrt(param_neff)
  cvmc <- 100 * abs(eemc / colMeans(param_samples))

  return(round(summary(eemc), 4))
}

CVMC("BETA", M2, neffM2)
CVMC("KAPPA2_JK", M2, neffM2)
CVMC("KAPPA2_K", M2, neffM2)

CVMC("S2_beta", M2, neffM2)
CVMC("S2_E", M2, neffM2)
CVMC("S2_M", M2, neffM2)
CVMC("S2_D", M2, neffM2)
CVMC("ALPHA_KAPPA", M2, neffM2)
CVMC("BETA_KAPPA", M2, neffM2)

EEMC_KAPPA2_jk  <- apply(X = M2$BETA, MARGIN = 2, FUN = sd)/sqrt(neffM2$BETA)
CVMC_KAPPA2_jk  <- 100*abs(EEMC_KAPPA2_jk/colMeans(M2$BETA))
round(summary(CVMC_KAPPA2_jk), 4)






## 11. Ranking Bayesiano ----

pdf("ranking.pdf", width = 8, height = 10)

deptos <- est_dpto$cole_depto_ubicacion

THETA_k_post <- colMeans(M1$THETA_K)

ci <- apply(M1$THETA_K, 2, quantile, probs = c(0.025, 0.975))

ord <- order(THETA_k_post)
deptos_ord      <- deptos[ord]
theta_means_ord <- THETA_k_post[ord]
ci_ord          <- ci[, ord]

colors <- rep("darkgray", m)
colors[ ci_ord[1, ] > 250 ] <- "springgreen3"
colors[ ci_ord[2, ] < 250 ] <- "tomato2"

# 5. Dibujar el gráfico
par(mar = c(4, 10, 1.5, 1), mgp = c(2.5, 0.75, 0), las = 1)

plot(NA, NA,
     xlab  = "Puntaje",
     ylab  = "",
     main  = expression(italic("Ranking Bayesiano")),
     xlim  = c(190, 280),
     ylim  = c(1, m),
     cex.axis = 0.8,
     yaxt     = "n")

axis(side = 2, at = 1:m, labels = deptos_ord, las = 2, cex.axis = 0.7)
abline(v = 250, col = "gray", lwd = 2)          
abline(h = 1:m, col = "lightgray", lwd = 1, lty = "dotted")


for (j in seq_len(m)) {
  segments(x0 = ci_ord[1, j], y0 = j,
           x1 = ci_ord[2, j], y1 = j,
           col = colors[j], lwd = 2)
  points(x = theta_means_ord[j], y = j,
         pch = 16, cex = 0.8,
         col = colors[j])
}

dev.off()

## 12 ----

B <- 5000

C <- matrix(0, nrow = m, ncol = m)
k_range <- 2:15
n_clusters <- length(k_range)


eps <- 0.1  



set.seed(123)
for (b in 1:B) {
  theta_b <- THETA[b, ]
  
  wss <- c()
  clusters_list <- list()
  
  for (i in seq_along(k_range)) {
    k  <- k_range[i]
    km <- kmeans(theta_b, centers = k)
    wss[i] <- km$tot.withinss
    clusters_list[[i]] <- km$cluster
  }
  
  # Calcular reducción relativa de WSS
  delta  <- abs(diff(wss)) / wss[-n_clusters]
  best_k <- which(delta < eps)
  
  if (length(best_k) == 0) {
    xi_b <- clusters_list[[1]]  # usar k = 2 si no hay gran mejora
  } else {
    xi_b <- clusters_list[[min(best_k) + 1]]
  }
  
  # Actualizar matriz de co-agrupamiento
  for (i in 1:(m-1)) {
    for (j in (i+1):m) {
      if (xi_b[i] == xi_b[j]) {
        C[i, j] <- C[i, j] + 1/B
      }
    }
  }
}

# Normalizar
C <- C + t(C)
diag(C) <- 1




## 16. Coeficientes de Regresión ----

mean(M2$BETA[,1])
quantile(M2$BETA[,1], probs = c(0.025, 0.975))


coef_names <- names(coef(fit))

colnames(M2$BETA) <- coef_names

### Individual ----

#### Intercepto ----

beta_intercepto <- M2$BETA[,1]
media <- mean(beta_intercepto)
LI <- quantile(beta_intercepto, probs = 0.025)
LS <- quantile(beta_intercepto, probs = 0.975)


color <- if (LI > 0 & media > 0) {
  "seagreen2"
} else if (LS < 0 & media < 0) {
  "tomato1"
} else {
  "black"
}

# Graficamos
pdf("efecto_intercepto.pdf", width = 6, height = 5)

par(mar = c(6, 5, 2.5, 1), mgp = c(2.5, 0.75, 0), las = 1)

# Seteamos el lienzo del plot
plot(NA, NA, 
     xlim = c(0.5, 1.5), 
     ylim = c(LI - 1, LS + 1),
     xlab = "", ylab = "Efecto del Intercepto", 
     xaxt = "n", cex.axis = 0.8)

# Ejes y guías
axis(side = 1, at = 1, labels = "Intercepto", cex.axis = 0.9)
abline(h = 0, col = "gray", lwd = 1.5)

# Segmento (intervalo de credibilidad)
segments(x0 = 1, y0 = LI, x1 = 1, y1 = LS, lwd = 3, col = color)

# Punto (media)
points(x = 1, y = media, pch = 16, cex = 1.2, col = color)

# Etiqueta de la media
text(x = 1, y = media, labels = round(media, 2), pos = 3, cex = 0.8)

dev.off()



#### Tablas ----

coef_ind <- colMeans(M2$BETA[,coef_names[c(2:14)]])

p_E <- length(coef_ind)

ci_ind <- apply(M2$BETA[,coef_names[c(2:14)]], 2, 
                quantile, probs = c(0.025, 0.975))


tabla_coef <- data.frame(
  Media = as.vector(coef_ind),
  LI         = ci_ind[1, ],
  LS         = ci_ind[2, ]
)


#### Variables binarias ----

colo <- rep("black", nrow(tabla_coef))
colo[(tabla_coef$LI > 0) & (tabla_coef$Media > 0)] <- "seagreen2"
colo[(tabla_coef$LS < 0) & (tabla_coef$Media < 0)] <- "tomato1"


pdf("efectos_binarios.pdf", width = 7, height = 5) 

#par(mar = c(6, 5, 2.5, 1), mgp = c(2.5, 0.75, 0), las = 1)

par(mar = c(4, 4, 1.5, 0.5), mgp = c(2, 0.5, 0), las = 1)

plot(NA, NA, 
     xlab = "", 
     ylab = "Efecto del Tratamiento", 
     xlim = c(0.5, 4.5), 
     ylim = range(tabla_coef[c(1:3,13), ]) + 3 * c(-1, 1),
     cex.axis = 0.75, xaxt = "n")

axis(side = 1, at = 1:4, labels = FALSE)

mtext(side = 1, at = 1:4, line = 1, cex = 0.7, 
      text = c("Madre con\neducación superior", 
               "Tiene\ncomputador", 
               "Tiene\ninternet", 
               "Pertenencia a\ngrupo étnico"))

abline(h = 0, col = "gray", lwd = 2)
abline(v = 1:4, col = "gray95", lwd = 1, lty = 2)

idx <- c(1:3, 13)

for (j in 1:4) {
  i <- idx[j]  # índice real
  
  segments(x0 = j, y0 = tabla_coef[i, 2], 
           x1 = j, y1 = tabla_coef[i, 3], 
           lwd = 3, col = colo[i])
  
  points(x = j, y = tabla_coef[i, 1], 
         pch = 16, cex = 0.8, col = colo[i])
  
  text(x = j, y = tabla_coef[i, 1], 
       labels = round(tabla_coef[i, 1], 2), 
       pos = 3, cex = 0.6)
}


dev.off()


#### Número de libros ----


pdf("efectos_libros.pdf", width = 7, height = 5) 

par(mar = c(4, 4, 1.5, 0.5), mgp = c(2, 0.5, 0), las = 1)

idx <- c(4:6)

plot(NA, NA, 
     xlab = "", 
     ylab = "Efecto del Tratamiento", 
     xlim = c(0.5, 3.5), ylim=range(tabla_coef[idx,])+3*c(-1, 1),
     cex.axis = 0.75, xaxt = "n")

axis(side = 1, at = 1:4, labels = FALSE)

mtext(side = 1, at = 1:4, line = 1, cex = 0.7, 
      text = c("11 a\n25 libros", 
               "26 a\n100 libros", 
               "Más de\n100 libros"))
abline(h = 0, col = "gray", lwd = 2)
abline(v = 1:4, col = "gray95", lwd = 1, lty = 2)

for (j in 1:3) {
  
  segments(x0 = j, y0 = tabla_coef[idx[j],2], 
           x1 = j, y1 = tabla_coef[idx[j],3], 
           lwd = 3, col = colo[j])
  
  
  lines(x = j, y = tabla_coef[idx[j],1], 
        type = "p", pch = 16, cex = 0.8, col = colo[j])
  
  
  text(x = j, y = tabla_coef[idx[j],1], 
       labels = round(tabla_coef[idx[j],1], 2), 
       pos = 3, cex = 0.6)
}


dev.off()


#### Estratos ----

pdf("efectos_estratos.pdf", width = 7, height = 5)
  
par(mar = c(4, 4, 1.5, 0.5), mgp = c(2, 0.5, 0), las = 1)


idx1 <- 7
idx2 <- 12
idx <- idx1:idx2
dif <- length(idx)

plot(NA, NA, 
     xlab = "", 
     ylab = "Efecto del Tratamiento", 
     xlim = c(0.5, dif + 0.5), 
     ylim = range(tabla_coef[idx, ]) + 3 * c(-1, 1),
     cex.axis = 0.75, xaxt = "n")

# Eje x
axis(side = 1, at = 1:dif, labels = FALSE)
mtext(side = 1, at = 1:dif, line = 1, cex = 0.7, 
      text = c(paste0("Estrato ", 2:6), "Sin Estrato"))
  # Líneas de referencia
abline(h = 0, col = "gray", lwd = 2)
abline(v = 1:dif, col = "gray95", lwd = 1, lty = 2)

# Dibujar puntos e intervalos con colores correctos
for (j in seq_len(dif)) {
  i <- idx[j]  # índice real dentro de tabla_coef
  
  segments(x0 = j, y0 = tabla_coef[i, 2], 
           x1 = j, y1 = tabla_coef[i, 3], 
           lwd = 3, col = colo[i])
  
  points(x = j, y = tabla_coef[i, 1], 
         pch = 16, cex = 0.8, col = colo[i])
  
  text(x = j, y = tabla_coef[i, 1], 
       labels = round(tabla_coef[i, 1], 2), 
       pos = 3, cex = 0.6)
}

dev.off()
  

### Municipal ----


coef_mun <- colMeans(M2$BETA[,coef_names[c(15:17)]])

p_M <- length(coef_mun)

ci_ind <- apply(M2$BETA[,coef_names[c(15:17)]], 2, 
                quantile, probs = c(0.025, 0.975))


tabla_coef <- data.frame(
  Media = as.vector(coef_mun),
  LI         = ci_ind[1, ],
  LS         = ci_ind[2, ]
)

colo <- rep("black", nrow(tabla_coef))
colo[(tabla_coef$LI > 0) & (tabla_coef$Media > 0)] <- "seagreen2"
colo[(tabla_coef$LS < 0) & (tabla_coef$Media < 0)] <- "tomato1"




pdf("efectos_mun.pdf", width = 7, height = 5)

par(mar = c(4, 4, 1.5, 0.5), mgp = c(2, 0.5, 0), las = 1)


idx1 <- 1
idx2 <- 3
idx <- idx1:idx2
dif <- length(idx)

plot(NA, NA, 
     xlab = "", 
     ylab = "Efecto del Tratamiento", 
     xlim = c(0.5, dif + 0.5), 
     ylim = range(tabla_coef[idx, ]) + 3 * c(-1, 1),
     cex.axis = 0.75, xaxt = "n")

# Eje x
axis(side = 1, at = 1:dif, labels = FALSE)
mtext(side = 1, at = 1:dif, line = 1, cex = 0.7, 
      text = c("NBI", 
               "Docentes\npor estudiante",
               "Índice de riesgo"))
# Líneas de referencia
abline(h = 0, col = "gray", lwd = 2)
abline(v = 1:dif, col = "gray95", lwd = 1, lty = 2)

# Dibujar puntos e intervalos con colores correctos
for (j in seq_len(dif)) {
  i <- idx[j]  # índice real dentro de tabla_coef
  
  segments(x0 = j, y0 = tabla_coef[i, 2], 
           x1 = j, y1 = tabla_coef[i, 3], 
           lwd = 3, col = colo[i])
  
  points(x = j, y = tabla_coef[i, 1], 
         pch = 16, cex = 0.8, col = colo[i])
  
  text(x = j, y = tabla_coef[i, 1], 
       labels = round(tabla_coef[i, 1], 2), 
       pos = 3, cex = 0.6)
}

dev.off()



### Departamental ----


coef_dep <- colMeans(M2$BETA[,coef_names[c(18:20)]])

p_D <- length(coef_mun)

ci_ind <- apply(M2$BETA[,coef_names[c(18:20)]], 2, 
                quantile, probs = c(0.025, 0.975))


tabla_coef <- data.frame(
  Media = as.vector(coef_dep),
  LI         = ci_ind[1, ],
  LS         = ci_ind[2, ]
)


colo <- rep("black", nrow(tabla_coef))
colo[(tabla_coef$LI > 0) & (tabla_coef$Media > 0)] <- "seagreen2"
colo[(tabla_coef$LS < 0) & (tabla_coef$Media < 0)] <- "tomato1"


pdf("efectos_dep.pdf", width = 7, height = 5)

par(mar = c(4, 4, 1.5, 0.5), mgp = c(2, 0.5, 0), las = 1)


idx1 <- 1
idx2 <- 3
idx <- idx1:idx2
dif <- length(idx)

plot(NA, NA, 
     xlab = "", 
     ylab = "Efecto del Tratamiento", 
     xlim = c(0.5, dif + 0.5), 
     ylim = range(tabla_coef[idx, ]) + 0.25*c(-1,1),
     cex.axis = 0.75, xaxt = "n")

# Eje x
axis(side = 1, at = 1:dif, labels = FALSE)
mtext(side = 1, at = 1:dif, line = 1, cex = 0.7, 
      text = c("PIB", 
               "Población\nRural",
               "Municipios en\nriesgo (%)"))
# Líneas de referencia
abline(h = 0, col = "gray", lwd = 2)
abline(v = 1:dif, col = "gray95", lwd = 1, lty = 2)

# Dibujar puntos e intervalos con colores correctos
for (j in seq_len(dif)) {
  i <- idx[j]  # índice real dentro de tabla_coef
  
  segments(x0 = j, y0 = tabla_coef[i, 2], 
           x1 = j, y1 = tabla_coef[i, 3], 
           lwd = 3, col = colo[i])
  
  points(x = j, y = tabla_coef[i, 1], 
         pch = 16, cex = 0.8, col = colo[i])
  
  text(x = j, y = tabla_coef[i, 1], 
       labels = round(tabla_coef[i, 1], 2), 
       pos = 3, cex = 0.6)
}

dev.off()







## 19. Predicciones ----

predict_data <- read.delim("predict_data.txt", sep = ";")

load("covariables_predict.RData")


### Modelo 1 ----

Xi_post <- colMeans(M1$XI)

Mun <- cbind(unique(datos$cole_cod_mcpio_ubicacion))

fitted_Mun <- data.frame(cod_mun = Mun, 
                          pred = Xi_post)

fitted_M1 <- (predict_data %>%
  left_join(fitted_Mun, by = c("cole_cod_mcpio_ubicacion"="cod_mun")))$pred

y <- predict_data$punt_global
ybar <- mean(y)

MSE_M1 <- mean((y - fitted_M1)^2)
MAE_M1 <- mean(abs(y - fitted_M1))
R2_M1 <- 1 -  sum((y - fitted_M1)^2)/sum((y - ybar)^2)


### Modelo 2 ----

fit_new <- lm(punt_global ~ edu_madre + computador + internet 
              + libros + estrato + etnia
              + nbi.m + doc_est.m + riesgo.m
              + PIB.d + pob_rural.d + riesgo.d, data = cov_predict)

X_new <- model.matrix(fit_new)

BETA_post <- colMeans(M2$BETA)

fitted_M2 <- as.numeric(X_new %*% BETA_post)

y <- cov_predict$punt_global
ybar <- mean(y)

MSE_M2 <- mean((y - fitted_M2)^2)
MAE_M2 <- mean(abs(y - fitted_M2))
R2_M2 <- 1 -  sum((y - fitted_M2)^2)/sum((y - ybar)^2)
