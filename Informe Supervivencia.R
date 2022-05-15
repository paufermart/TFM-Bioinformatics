#' ---
#' title: "Análisis de Supervivencia con algoritmos de Machine Learning"
#' author: "Paula Fernández Martínez"
#' date: '`r format(Sys.time(), "%d %B de %Y")`'
#' output:
#'   rmdformats::readthedown:
#'    self.contained: no
#'    highlight: tango
#'    code_folding: hide
#' toc_depth: 4
#' ---

#' <style>
#'   body {
#'     text-align: justify}
#' </style>

#+ include = FALSE
#+ message = FALSE, warning = FALSE

library(data.table)
library(readr)
library(ggplot2)
library(survival)
library(survminer)
library(caret)
library(gdata)
library(auxfun)
library(mice)
library(ggmice)
library(randomForestSRC)
library(ggRandomForests)
library(Hmisc)
library(scales)
library(dplyr)
library(e1071)
library(survivalsvm)
library(kableExtra)
library(compareC)

# Carga de datos
lake <- read_delim("~/Máster Bioinformática/TFM/DATOS/lake.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                   trim_ws = TRUE)

lake3 <- readxl::read_excel("~/Máster Bioinformática/TFM/DATOS/lake2.xlsx")
names(lake3)[32] <- "tiempo"

#' # Limpieza de datos
#' 
#' El primer paso antes de realizar el análisis es comprobar cuántos valores perdidos tenemos y si es necesario 
#' llevar a cabo alguna transformación.
#' 
#' Comenzamos visualizando el patrón general de NAs:

plot_pattern(lake, rotate = TRUE)
#' 
#' Calculamos el % en la base de datos:
#' 
# % NAs en toda la base de datos
porc_NA <- round((sum(is.na(lake)) / (ncol(lake) * nrow(lake))) *100,2)
porc_NA

#' Tenemos un `r porc_NA`% de valores perdidos.
#' 
#' A continuación eliminamos algunas variables no relevantes. Este paso dependerá de los datos con los que
#' trabajemos.
#' 
#' En este caso, eliminamos la procedencia de los datos, el nombre del paciente, el número de paciente, 
#' el número de visita, la fecha de nacimiento (ya tenemos una variabe edad) y los factores de riesgo individuales,
#' que ya que se encuentran en la variable factor_riesgo_total. También la fecha inicio Lake, las fechas 
#' correspondientes a las semanas de seguimiento y la fecha VIH porque ya tenemos el tiempo en meses.
#' 

lake <- lake[,-c(1:5,7,10:15,17,18,19,23:24,60:61,97:98,137:138,174:175, 174:175)]

#' Solucionamos las incongruencias con los tiempos de infección, que muestra valores negativos:
#' 
summary(lake$tpo_vih_meses)

lake$tpo_vih_meses <- ifelse(lake$tpo_vih_meses <=0, NA, lake$tpo_vih_meses)

#' Una vez solucionado nos queda así:
summary(lake$tpo_vih_meses)

#' Calculamos el porcentaje de NAs por variable y eliminamos aquellas que presentan un 20% o más.
#+ message = FALSE, warning = FALSE

porc_missing <- round(unlist(lapply(lake, function(x) sum(is.na(x))))/
                        nrow(lake)*100,2)
porc_missing

dd_limpio <- remove.vars(lake,names(porc_missing[porc_missing >=20]))

porc_NA <- round((sum(is.na(dd_limpio)) / (ncol(dd_limpio) * nrow(dd_limpio))) *100,2)
porc_NA

#' Ahora tenemos un `r porc_NA`% de NAs general y el patrón gráfico nos queda como sigue:
#' 
#+ message = FALSE
plot_pattern(dd_limpio, rotate = TRUE)
gg <- plot_pattern(dd_limpio, rotate = TRUE)
ggtoppt(gg)

#'
#' En este caso tuvimos que calcular el tiempo de seguimiento y determinar la variable éxito/fracaso aparte de la base
#' de datos, por lo que las incoportamos al data.frame:
#' 
#Incorporamos el tiempo de seguimiento y el status

dd_limpio$tiempo <- lake3$tiempo #contamos el tiempo de seguimiento total

dd_limpio$status <- lake3$status_sin0

dd_limpio <- dd_limpio[!is.na(dd_limpio$tiempo),]

dd_limpio$status <- ifelse(dd_limpio$status == 1, 0, 1)

dd_limpio$status <- ifelse(is.na(dd_limpio$status) == TRUE, 1, dd_limpio$status)

dd_limpio$edad <- ifelse(dd_limpio$edad == 0, NA, dd_limpio$edad)

table(dd_limpio$tiempo)

table(dd_limpio$status)

#' # Análisis descriptivo
#' 
#' ## Tranformación de los datos
#' 
#' Comenzamos modificando etiquetas para que las tablas y los gráficos sean más comprensibles.
#' 
lake2_2 <- dd_limpio

lake2_2$factor_riesgo_total <- factor(lake2_2$factor_riesgo_total, levels = 1:5, 
                                      labels = c("ADVP", "Heterosexual","Homosexual", "Hemofilia", "Otros"))

lake2_2$Grupo <- factor(lake2_2$Grupo, levels = c(-1,0), labels = c("EFV + Kivexa", "Kaletra + Kivexa"))

lake2_2$sexo <- factor(lake2_2$sexo, levels = c(1,2), labels = c("Hombre", "Mujer"))

lake2_2$tiempo <- as.factor(lake2_2$tiempo)

lake2_2$status <- as.factor(lake2_2$status)

lake2_2$edad <- ifelse(lake2_2$edad == 0, NA, lake2_2$edad)

#' Generamos un par de tablas con los datos más relevantes:
#' 
table1::table1( ~ sexo + edad + tpo_vih_meses | Grupo, lake2_2, overall = "Total")

table1::table1( ~ tiempo + status | Grupo, lake2_2, overall = "Total")

#' Y un gráfico con la evolución de la carga viral:
#' 
#+ message = FALSE
lake_cv50 <- lake[, c("cv50_0", "cv50_12", "cv50_24", "cv50_36", "cv50_48")]

lake_cv50 <- lapply(lake_cv50, factor)

lake_cv50$id <- 1:116

long_cv <- melt(setDT(lake_cv50), id.vars = c("id"), variable.name = "Carga_viral")

long_cv$value <- factor(long_cv$value, levels = 0:1,
                        labels = c("Superior al límite", "Indetectable"))

long_cv2 <- na.omit(long_cv) %>% 
  group_by(Carga_viral, value) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

ggplot(long_cv2, aes(x = Carga_viral, y = perc*100, fill = value)) +
  geom_bar(stat = "identity", position = 'dodge') +
  ylim(0,100) +
  labs(y = "%", x = "Semana", fill = "Carga viral")

gg <- ggplot(long_cv2, aes(x = Carga_viral, y = perc*100, fill = value)) +
  geom_bar(stat = "identity", position = 'dodge') +
  ylim(0,100) +
  labs(y = "%", x = "Semana", fill = "Carga viral")

ggtoppt(gg)

ggtoppt(export = TRUE, "graficos.pptx")

#' Generamos unas curvas de Kaplan Meier para ver cómo funcionan los modelos tradicionales.
#' 
#' Hacemos una con todas las variables:
#+ message = FALSE, warning = FALSE

sur1 <- survfit(Surv(tiempo, status) ~ 1, data = dd_limpio)

ggsurvplot(sur1,
           risk.table = T,
           surv.median.line = "hv",
           linetype = "strata",
           break.time.by = 12,
           xlab = "Semanas",
           title = "Todos")

#' Y otra por grupos:
#+ message = FALSE, warning = FALSE
lake_surv <- dd_limpio[, c("sexo", "edad", "Grupo", "tpo_vih_meses", "tiempo", "status")]

lake_surv$Grupo <- factor(lake_surv$Grupo, levels = c(-1,0), labels = c("EFV + Kivexa", "Kaletra + Kivexa"))
sur2 <- survfit(Surv(tiempo, status) ~ Grupo, data = lake_surv)
ggsurvplot(sur2, 
           risk.table = T, 
           surv.median.line = "hv", 
           pval = TRUE, 
           break.time.by = 12, 
           xlab = "Semanas",
           title = "")

#' # Algoritmos de Machine Learning
#' 
#' ## Preparación de los datos
#' 
#' Comenzamos dividiendo la base de datos en dos partes: entrenamiento (67%) y test (33%).
#' 
dd_limpio_ml <- dd_limpio
dd_limpio_ml$tiempo <- ifelse(dd_limpio_ml$tiempo == 0 , 0.0, dd_limpio_ml$tiempo)
dd_limpio_ml$Grupo <- ifelse(dd_limpio_ml$Grupo == -1, 1, 0)

#Implantamos la semilla
set.seed(123)

#Separamos las dos partes creando un objeto que contenga el 67% de los datos.
objeto1 <- sample(nrow(dd_limpio), nrow(dd_limpio)*0.67)

#Al subset train le añadimos el 67% de los datos
train <- dd_limpio[objeto1, ]

#Al subset test le añadimos el resto (33%)
test <- dd_limpio[-objeto1, ]

#' ## Árbol de supervivencia
#'Primero realizamos un modelo con 100 árboles.
#' 
train_forest <- train

test_forest <- test

rf_1 <- rfsrc(Surv(tiempo, status) ~ ., 
              data = train_forest,
              ntree = 100,
              na.action = "na.impute",
              tree.err = TRUE,importance = TRUE)

rf_1


#' Representamos la estimación del error de predicción en función del número de árboles en el bosque:
#' 
plot(gg_error(rf_1)) 

#'
#' Hacemos las predicciones y volvemos a representar gráficamente el error:

rf_test <- predict(rf_1, newdata = test_forest, na.action = "na.impute", importance = TRUE)

rf_test

plot(gg_error(rf_test)) 

c_100rf <- compareC::estC(test_forest$tiempo, test_forest$status, rf_test$predicted)

#' 
#' También podemos representar qué variables son las que más influyen en los resultados: en azul aparecen las variables que reducen
#' el error y en rojo las que lo aumentan.
#' 
plot(gg_vimp(rf_test))
#'
#' El índice C de Harrell para este modelo es: `r round(c_100rf,4)`
#' 
#' Repetimos el modelo pero aumentando el número de árboles:
#' 
rf_2 <- rfsrc(Surv(tiempo, status) ~ ., 
              data = train_forest,
              ntree = 300,
              na.action = "na.impute",
              tree.err = TRUE,importance = TRUE)

rf_2

plot(gg_error(rf_2))

rf_test2 <- predict(rf_2, newdata = test_forest, na.action = "na.impute", importance = TRUE)

rf_test2

plot(gg_error(rf_test2))

plot(gg_vimp(rf_test2))

c_300rf <- compareC::estC(test_forest$tiempo, test_forest$status, rf_test2$predicted)

#' 
#' Ahora el índice C es `r round(c_300rf,4)`.
#' 
#' ## Naive Bayes
#' 
#' De este algoritmo también vamos a hacer dos versiones: una sin la corrección de Laplace y otra con ella.
#' 
#' ###  Laplace = 0
bayes_lap0_model<-naiveBayes(Surv(tiempo,status) ~ ., data=train, laplace=0, na.action = na.pass)

summary(bayes_lap0_model)

bayes_lap0_pred <- predict(bayes_lap0_model, test)

train_labels <- Surv(train$tiempo, train$status)
test_labels <- as.factor(Surv(test$tiempo,test$status))

c_nb0 <- compareC::estC(test$tiempo, test$status, bayes_lap0_pred)

#' Obtenemos un índice C de `r round(c_nb0,4)`.
#' 
#' ### Laplace 1
#' 
bayes_lap1_model <- naiveBayes(Surv(tiempo,status) ~ ., data=train, laplace=1, na.action = na.pass)

summary(bayes_lap1_model)
bayes_lap1_pred <- predict(bayes_lap1_model, test)

c_nb1 <- compareC::estC(test$tiempo, test$status, bayes_lap1_pred)

#' Obtenemos un índice C de `r round(c_nb1,4)`.
#' 
#' 
#' ## Support Vector Machine 
#' 
#' En el caso del SVM, debemos eliminar los NAs para que no den problemas con los modelos.
#' 
#Implantamos la semilla
set.seed(123)

dd_svm <- na.omit(dd_limpio_ml)

objeto_svm <- sample(nrow(dd_svm), nrow(dd_svm)*0.67)

train_svm <- dd_svm[objeto_svm,]

test_svm <- dd_svm[-objeto_svm,]
#' 
#' Vamos a probar tres variantes:
#' 
#' ### Kernel lineal
surv_lin <- survivalsvm(Surv(tiempo, status) ~.,
                        data = train_svm,
                        type = "regression", gamma.mu = 1,
                        kernel = "lin_kernel")

print(surv_lin)

pred_surv_lin <- predict(surv_lin, newdata = test_svm)

print(pred_surv_lin)

c_svmlin <- compareC::estC(test_svm$tiempo, test_svm$status, pred_surv_lin$predicted)

#' Obtenemos una C de `r round(c_svmlin,4)`
#'
#' ### Kernel radial

surv_rbf <- survivalsvm(Surv(tiempo, status) ~.,
                        data = train_svm,
                        type = "regression", gamma.mu = 1,
                        opt.meth = "quadprog", kernel = "rbf_kernel")

print(surv_rbf)

pred_surv_rbf <- predict(surv_rbf, test_svm)

print(pred_surv_rbf)

c_svmrbf <- compareC::estC(test_svm$tiempo, test_svm$status, pred_surv_rbf$predicted)

#' Obtenemos una C de `r round(c_svmrbf,4)`
#' 
#' ### Kernel aditivo
#' 

surv_add <- survivalsvm(Surv(tiempo, status) ~.,
                        data = train_svm[,2:39],
                        type = "regression", gamma.mu = 1,
                        opt.meth = "quadprog", kernel = "add_kernel")

print(surv_add)

pred_surv_add <- predict(surv_add, newdata = test_svm)

print(pred_surv_add)

conindex(pred_surv_add, test_svm$tiempo)

c_svmadd <- compareC::estC(test_svm$tiempo, test_svm$status, pred_surv_add$predicted)

#' Obtenemos una C de `r round(c_svmadd,4)`
#' 
#' # Resumen de los resultados
#' 
#' En la siguiente tabla se muestran los valores de C para los diferentes algoritmos:
#' 
c_values <- data.frame(names = c("Árbol de Supervivencia n = 100", "Árbol de supervivencia n = 300",
                                 "Naive Bayes L = 0", "Naive Bayes L = 1",
                                 "SVM Lineal", "SVM RBF", "SVM Aditivo"),
                       values = c(c_100rf, c_300rf, c_nb0, c_nb1, c_svmlin, c_svmrbf, c_svmadd))

kable(c_values, caption = "Resumen de los valores de C para todos los algoritmos", 
    col.names = c("Algoritmo", "C-Index"),
    digits = 4,
    format = "html", table.attr = "style='width:100%;'", align = "c") %>%
  kable_styling(bootstrap_options = c("striped", "hover","condensed"))

max_c <- max(c_values$values)

#' El algoritmo con mayor valor de C es `r c_values[c_values$values == max_c,][[1]]`, con un índice C de `r round(max_c,4)`.
#' 
#' 
