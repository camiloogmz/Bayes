library(dplyr)
library(readxl)
library(sqldf)

setwd("C:/Users/camil/OneDrive - Universidad Nacional de Colombia/Documentos/Camilo UNAL/Estadística Bayesiana/Caso de Estudio/Modelo 2")

datos <- read.delim("Examen_Saber_11_2022_2_Prueba.txt", sep = ";")

datos <- subset(datos, estu_nacionalidad == "COLOMBIA" & estu_pais_reside == "COLOMBIA"
                       & cole_depto_ubicacion != "SAN ANDRES")

datos <- datos[, c("punt_global",
                                 "fami_educacionmadre",
                                 "fami_tienecomputador",
                                 "fami_tieneinternet",
                                 "fami_numlibros",
                                 "fami_estratovivienda",
                                 "estu_tieneetnia",
                                 "estu_cod_reside_depto",
                                 "estu_depto_reside",
                                 "estu_cod_reside_mcpio",
                                 "estu_mcpio_reside",
                                 "cole_depto_ubicacion",
                                 "cole_cod_depto_ubicacion",
                                 "cole_mcpio_ubicacion",
                                 "cole_cod_mcpio_ubicacion")]

datos <- na.omit(datos)
write.table(datos, "predict_data.txt", sep = ";")





datos <- read.delim("predict_data.txt", sep = ";")

datos <- datos %>% arrange(cole_cod_depto_ubicacion, cole_mcpio_ubicacion)

## Covariables ##

# A nivel estudiante ----

setwd("C:/Users/camil/OneDrive - Universidad Nacional de Colombia/Documentos/Camilo UNAL/Estadística Bayesiana/Caso de Estudio")


datos$fami_educacionmadre <- ifelse(datos$fami_educacionmadre %in% 
                                      c("Postgrado", 
                                        "Educación profesional completa"),
                                    1, 0)
datos$fami_tienecomputador <- ifelse(datos$fami_tienecomputador=="Si", 1, 0)

datos$fami_tieneinternet <- ifelse(datos$fami_tieneinternet=="Si", 1, 0)


# A nivel municipal ----

## NBI municipal ----

car_gen <- read_xlsx("PANEL_CARACTERISTICAS_GENERALES(2022).xlsx")
car_gen <- subset(car_gen,
                    ano %in% c(2018)) # unicamente hay datos en 2018 (mas reciente)
car_gen <- sqldf("select coddepto, depto, codmpio, municipio, ano, nbi, pobl_tot, pobl_rur
                   from car_gen")

deptos <- car_gen[,1:4]

nbi <- sqldf("select codmpio, municipio, nbi
              from car_gen
              order by depto, municipio")

## Docentes por estudiante ----

load("edu.RData")
#edu <- read_xlsx("PANEL_DE_EDUCACION(2022).xlsx")
edu_mun <- subset(edu, ano == 2021)
edu_mun <- edu_mun[,c("codmpio", "docen_total", "alumn_total")]

#edu_mun <- sqldf("select e.*, municipio, depto
#                 from edu_mun e
#                  inner join deptos d
#                  on e.codmpio = d.codmpio")

edu_mun <- sqldf("select codmpio, docen_total/alumn_total as 'doc.est'
                  from edu_mun")

##########  edu_mun <- edu_mun %>% arrange(depto, municipio)

## Índice de riesgo por hechos victimizantes ----

LAFT <- read_excel("LAFT.xlsx",
                   sheet = "FT",
                   range = "A2:B1123",
                   col_names = c("DPTO_MUN","RISK_VICTIM_2022"))

LAFT <- LAFT %>% arrange(DPTO_MUN)

LAFT.m <- cbind(deptos %>% arrange(depto, municipio), riesgo.m = LAFT$RISK_VICTIM_2022)


# A nivel departamental ----


## PIB ----

pib_num <- read_excel("DANE - PIB.xlsx", 
                  sheet = 4,
                  range = "T11:T43",
                  col_names = "PIB")
pib_deptos <- read_excel("DANE - PIB.xlsx", 
                         sheet = 4,
                         range = "A11:B43",
                         col_names = c("Código", "Departamento"),
                         col_types = c("numeric", "guess"))

pib <- cbind(pib_deptos, pib_num); pib <- pib[pib$Código != 88,] 

pib <- pib %>% arrange(Departamento) %>% mutate(PIB=PIB/10^6)


## Población rural ----

pob_rural <- sqldf("select coddepto, depto, (sum(pobl_rur)/sum(pobl_tot))*100 as PobRural
                    from car_gen
                    group by coddepto, depto
                    order by depto")
pob_rural <- pob_rural[pob_rural$coddepto != 88,]

##  Porcentaje de municipios en riesgo de violencia ----


deptos_unique <- deptos %>% distinct(coddepto, depto)
deptos_unique <- deptos_unique[deptos_unique$coddepto != 88, ] # Sin San Andrés

df_pct <- data.frame(
  coddepto = c(5, 8, 11, 13, 15, 17, 18, 19, 20, 23,
               25, 27, 41, 44, 47, 50, 52, 54, 63, 66,
               68, 70, 73, 76, 81, 85, 86, 88, 91, 94,
               95, 97, 99),
  pct_riesgo = c(48.00, 17.39, 100.00, 45.65, 1.63,
                 3.70, 93.75, 76.19, 64.00, 53.33,
                 1.72, 90.00, 18.92, 60.00, 20.00,
                 62.07, 46.88, 47.50, 8.33, 28.57,
                 6.90, 69.23, 19.15, 47.62, 100.00,
                 21.05, 76.92, 0.00, 0.00, 55.56,
                 75.00, 0.00, 50.00),
  stringsAsFactors = FALSE
)




riesgo_dep <- df_pct %>%
                right_join(deptos_unique, by = "coddepto") %>%
                select(coddepto, depto, pct_riesgo) %>%
                arrange(depto)



# Arreglo Final ----

est_dpto <- datos %>%
  group_by(cole_depto_ubicacion, cole_cod_depto_ubicacion) %>%
  summarise(
    n_k = n_distinct(cole_cod_mcpio_ubicacion),
    n_ijk = n()
  ) %>%
  ungroup() %>%
  left_join(pib, by = c("cole_cod_depto_ubicacion" = "Código")) %>%
  left_join(pob_rural, by = c("cole_cod_depto_ubicacion" = "coddepto")) %>%
  left_join(riesgo_dep, by = c("cole_cod_depto_ubicacion" = "coddepto")) %>%
  select(cole_depto_ubicacion, cole_cod_depto_ubicacion, n_ijk, PIB, PobRural, pct_riesgo)


est_mun <- datos %>%
  group_by(cole_mcpio_ubicacion, cole_cod_mcpio_ubicacion,
           cole_depto_ubicacion, cole_cod_depto_ubicacion) %>%
  summarise(
    n_jk = n(),
    .groups = "drop"
  ) %>%
  arrange(cole_depto_ubicacion, cole_mcpio_ubicacion) %>%
  left_join(edu_mun, by = c("cole_cod_mcpio_ubicacion"="codmpio")) %>%
  left_join(nbi, by = c("cole_cod_mcpio_ubicacion" = "codmpio")) %>%
  left_join(LAFT.m, by = c("cole_cod_mcpio_ubicacion" = "codmpio")) %>%
  select(cole_cod_depto_ubicacion, cole_depto_ubicacion, cole_cod_mcpio_ubicacion,
         cole_mcpio_ubicacion, n_jk, nbi, doc.est, riesgo.m)

cov_predict <- datos %>% 
          left_join(est_mun,
                    by = c("cole_cod_mcpio_ubicacion"="cole_cod_mcpio_ubicacion")) %>%
          left_join(est_dpto, 
                    by = c("cole_cod_depto_ubicacion.x"="cole_cod_depto_ubicacion")) %>%
          select(punt_global, 
                 cole_depto_ubicacion = cole_depto_ubicacion, 
                 cole_cod_depto_ubicacion = cole_cod_depto_ubicacion.x,
                 cole_mcpio_ubicacion = cole_mcpio_ubicacion.x,
                 cole_cod_mcpio_ubicacion = cole_cod_mcpio_ubicacion,
                 edu_madre = fami_educacionmadre, 
                 computador = fami_tienecomputador,
                 internet = fami_tieneinternet,
                 libros = fami_numlibros,
                 estrato = fami_estratovivienda,
                 etnia = estu_tieneetnia,
                 nbi.m = nbi,
                 doc_est.m = doc.est,
                 riesgo.m = riesgo.m,
                 PIB.d = PIB,
                 pob_rural.d = PobRural,
                 riesgo.d = pct_riesgo) %>%
                 mutate(
                  across(
                    .cols = c(edu_madre, computador, internet, libros, estrato, etnia),
                    .fns  = as.factor
                  )
                 ) %>%
                na.omit()
  


setwd("C:/Users/camil/OneDrive - Universidad Nacional de Colombia/Documentos/Camilo UNAL/Estadística Bayesiana/Caso de Estudio/Solución Puntos")

save(cov_predict, file = "covariables_predict.RData")
