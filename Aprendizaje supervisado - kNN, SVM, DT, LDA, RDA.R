'''
https://depmap.org/portal/data_page/?tab=allData

OmicsExpressionProteinCodingGenesTPMLogp1.csv
Model.csv


Los datos que compartes pertenecen a un conjunto de datos de expresión génica de la plataforma DepMap (Dependency Map), lanzado en su versión 24Q2. Esta plataforma, mantenida por el Broad Institute, proporciona datos sobre la expresión génica y otras características de las líneas celulares, basados en mediciones de RNA-seq. En este caso, los valores de expresión génica están representados en TPM (Transcripts Per Million), y luego transformados con el logaritmo en base 2, utilizando un valor de pseudo-cuenta de 1 para evitar valores cero en los cálculos (es decir, log2(TPM + 1)).

Detalles clave :

--> Datos de expresión génica: Miden los niveles de expresión de genes codificantes de proteínas en las líneas celulares de DepMap. Esta información se obtiene a partir de datos de secuenciación de RNA (RNA-seq).

--> Log2 transformación: La transformación logarítmica en base 2 se usa comúnmente en los datos de expresión génica para normalizar los datos y hacer que las diferencias de expresión sean más manejables, especialmente cuando los valores varían en órdenes de magnitud.

--> Uso de valores TPM: El valor TPM (Transcripts Per Million) se usa para cuantificar la expresión de los genes. Se normaliza por millón de lecturas, lo que facilita la comparación entre diferentes genes y muestras.

--> CRISPR y otras herramientas: Además de los datos de expresión, la plataforma DepMap incluye información sobre la dependencia genética de las células en función de la perturbación genética, como las pantallas CRISPR, pantallas de drogas PRISM, y características relacionadas con el número de copias, mutaciones y fusiones de genes.
Referencias y recursos adicionales: Este conjunto de datos es parte de una liberación pública de DepMap que está disponible para investigadores interesados en estudiar la biología celular, las interacciones genéticas y cómo los fármacos pueden afectar la expresión génica de diversas líneas celulares. Más detalles sobre los datos y el proceso de análisis pueden encontrarse en el repositorio de GitHub de DepMap.

'''

#### ---- Preparamos la base de datos para modelos de clasificación ----

rm(list=ls())
getwd()
setwd("~/R/Algoritmo")
df <- read.csv("dataset_expresiongenes_cancer.csv")
str(df)
table(df$primaryormetastasis)
data <- subset(df, select = 3:101)
str(data)

install.packages("glmnet")
install.packages("rpart.plot")
install.packages("rattle")
install.packages("pROC")
install.packages("PRROC")
install.packages("MASS")
install.packages("klard")
library(glmnet) # ElasticNet
library(tidyverse)
library(caret) # ML
library(rpart) # DT
library(rpart.plot) # DT plot
library(rattle) # DT plot
library(pROC) # ROC
library(PRROC) # PR-Curve
library(MASS) # LDA
library(klaR) # RDA
library(gridExtra) # juntar los gráficos


# Machine learning methods: http://topepo.github.io/caret/train-models-by-tag.html
# in method="XXXX"
# in metric="XXXX" -> "RMSE" para regression y "Accuracy" para classification
# Apoyo: https://rpubs.com/nomarpicasso/1150167
names(getModelInfo())
modelLookup(model = "knn")
modelLookup(model = "svmLinear")
modelLookup(model = "svmRadial")
modelLookup(model = "svmPoly")
modelLookup(model = "rpart")








genes <- names(df[4:19197])


# Preparar los datos para el modelo LASSO - ridge (mantiene todas las varibles, incluso =0),lasso(si
#elimina las =0), elasticNet (en funcion de los datos el mantiene un intermedio entre las dos anteriores)
x <- as.matrix(df[, genes])
y <- factor(df$primaryormetastasis)

# Ajustar el modelo LASSO
set.seed(1995)
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)
selected_genes <- coef(lasso_model, s = "lambda.min")
selected_genes <- as.matrix(selected_genes) # Convertir a matriz densa si es necesario (esto convierte el formato disperso a un formato manejable)
selected_genes_df <- as.data.frame(selected_genes) # Convertir la matriz a data frame para facilitar su manipulación
non_zero_indices <- selected_genes_df[selected_genes_df$s1 != 0, , drop = FALSE] # Filtrar los coeficientes no cero
dim(non_zero_indices)
non_zero_indices

names <- rownames(non_zero_indices)[2:180]
names

data <- df %>% dplyr:: select(primaryormetastasis, names)
rows <- df$code
rownames(data) <- rows
str(data)
#extraer bases de datos en formato excel
write_csv(data, "/Users/vic/Library/CloudStorage/GoogleDrive-vdelaopascual@gmail.com/Mi unidad/MU en Bioinformática (UNIR 2023)/Presentaciones Algoritmos e Inteligencia Artificial (AAAA-MM-DD)/stata.csv")
names <- colnames(data)[-1]
data <- data %>% dplyr::select(names, primaryormetastasis)


# Dividir el conjunto de datos en conjuntos de entrenamiento y prueba
set.seed(1995)
trainIndex <- createDataPartition(data$primaryormetastasis, p = 0.8, list = FALSE)
data$primaryormetastasis <- as.factor(data$primaryormetastasis)
trainData <- data[trainIndex,]
testData <- data[-trainIndex,]





#### ---- kNN ----
# Crear un modelo de k-NN utilizando el paquete caret (el . hace referencia a toda la base de datos, cv validacion
# )
knnModel <- train(primaryormetastasis ~ .,
                  data = trainData,
                  method = "knn",
                  trControl = trainControl(method = "cv", number = 10),
                  preProcess = c("center", "scale"),
                  tuneLength = 15)

knnModel

plot(knnModel)

# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
predictions <- predict(knnModel, newdata = testData )
predictions

# Evaluar la precisión del modelo utilizando la matriz de confusión
confusionMatrix(predictions, testData$primaryormetastasis)


# Obtener probabilidades
probabilities_knn <- predict(knnModel, newdata = testData, type = "prob")
probabilities_knn




#### ---- Support Vector Machine ----
# Crear un modelo de SVM lineal utilizando el paquete caret
# parámetro C "cost" por defecto es 1, pero puedes tunearlo. Controla la flexibilidad del modelo para encontrar un equilibrio entre un margen amplio y la clasificación correcta de las muestras
svmModelLineal <- train(primaryormetastasis ~.,
                        data = trainData,
                        method = "svmLinear",
                        trControl = trainControl(method = "cv", number = 10),
                        preProcess = c("center", "scale"),
                        tuneGrid = expand.grid(C = seq(0, 2, length = 20)), #C grande lleva al sobreajuste, C pequeño al infraajuste
                        prob.model = TRUE) 
#expand.grid para expandir la busqueda entre 0 y 2 y 20 veces
svmModelLineal

plot(svmModelLineal)

# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
predictions <- predict(svmModelLineal, newdata = testData )
predictions

# Evaluar la precisión del modelo utilizando la matriz de confusión
confusionMatrix(predictions, testData$primaryormetastasis)

# SVM lineal
probabilities_svm_linear <- predict(svmModelLineal, newdata = testData, type = "prob")
probabilities_svm_linear

#Knn daba o 0 a 1, el lineal es más fino, no redondea como el knn


# Crear un modelo de SVM tipo kernel utilizando el paquete caret
# no hace falta tunear el parámetro C "cost" 
svmModelKernel <- train(primaryormetastasis ~.,
                        data = trainData,
                        method = "svmRadial",
                        trControl = trainControl(method = "cv", number = 10),
                        preProcess = c("center", "scale"),
                        tuneLength = 10, #numero de repeticiones
                        prob.model = TRUE) 
svmModelKernel

plot(svmModelKernel)


# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
predictions <- predict(svmModelKernel, newdata = testData )
predictions

# Evaluar la precisión del modelo utilizando la matriz de confusión
confusionMatrix(predictions, testData$primaryormetastasis)

# SVM kernel
probabilities_svm_kernel <- predict(svmModelKernel, newdata = testData, type = "prob")
probabilities_svm_kernel



# Crear un modelo de SVM tipo kernel polynomial utilizando el paquete caret
# no hace falta tunear el parámetro C "cost" 
svmModelKernelPolynomial <- train(primaryormetastasis ~.,
                                  data = trainData,
                                  method = "svmPoly",
                                  trControl = trainControl(method = "cv", number = 10),
                                  preProcess = c("center", "scale"),
                                  tuneLength = 5,
                                  prob.model = TRUE) 
svmModelKernelPolynomial

plot(svmModelKernelPolynomial)


# Realizar predicciones en el conjunto de prueba utilizando el modelo entrenado
predictions <- predict(svmModelKernelPolynomial, newdata = testData )
predictions

# Evaluar la precisión del modelo utilizando la matriz de confusión
confusionMatrix(predictions, testData$primaryormetastasis)

# SVM kernel
probabilities_svm_kernelpol <- predict(svmModelKernelPolynomial, newdata = testData, type = "prob")
probabilities_svm_kernelpol



#### ---- Decission Tree ----
# Crear un modelo de DT utilizando el paquete caret
dtModel <- train(primaryormetastasis ~.,
                 data = trainData,
                 method = "rpart",
                 trControl = trainControl(method = "cv", number = 10),
                 preProcess = c("center", "scale"),
                 tuneLength = 10)
dtModel
plot(dtModel)

fancyRpartPlot(dtModel$finalModel, type=4)


# Evaluar el modelo con el conjunto de prueba
predictions_raw <- predict(dtModel, newdata = testData, type = "raw") # raw = clases
predictions_raw


# Evaluar la precisión del modelo utilizando la matriz de confusión
confusionMatrix(predictions_raw, testData$primaryormetastasis)

# Obtener probabilidades
probabilities_dt <- predict(dtModel, newdata = testData, type = "prob")







#### ---- LDA (lineal) ----
formula <- as.formula(paste("primaryormetastasis ~", paste(names, collapse = "+")))
formula

# Ajustar el modelo LDA en el entrenamiento
lda_model <- lda(formula, data = trainData)
lda_model$scaling # contribuciones/coeficientes

# Realizar predicciones sobre el conjunto de prueba
lda_predictions <- predict(lda_model, newdata = testData)
lda_predictions$x

predicted_classes <- lda_predictions$class # Obtener la predicción (predicciones de la clase)
true_classes <- as.factor(testData$primaryormetastasis) # Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
confusionMatrix(predicted_classes, true_classes) # Crear la matriz de confusión (tumor predicho testing vs. tumor real testing)
probabilities_lda <- predict(lda_model, newdata = testData, type = "prob") # Obtener probabilidades






#### ---- RDA (regularizado) ----
# Ajustar el modelo RDA en el entrenamiento
rda_model <- rda(formula, data = trainData)

# Realizar predicciones sobre el conjunto de prueba
rda_predictions <- predict(rda_model, newdata = testData)

predicted_classes <- rda_predictions$class # Obtener la predicción (predicciones de la clase)
true_classes <- as.factor(testData$primaryormetastasis) # Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
confusionMatrix(predicted_classes, true_classes) # Crear la matriz de confusión (tumor predicho testing vs. tumor real testing)
probabilities_rda <- predict(rda_model, newdata = testData, type = "prob") # Obtener probabilidades







#### ---- Curvas ROC ----
roc_knn <- roc(testData$primaryormetastasis, probabilities_knn[,2]) # Cambia [,2] según la clase positiva
auc_knn <- auc(roc_knn)
cat("AUC k-NN:", auc_knn, "\n")

roc_svm_linear <- roc(testData$primaryormetastasis, probabilities_svm_linear[,2])
auc_svm_linear <- auc(roc_svm_linear)
cat("AUC SVM Lineal:", auc_svm_linear, "\n")

roc_svm_kernel <- roc(testData$primaryormetastasis, probabilities_svm_kernel[,2])
auc_svm_kernel <- auc(roc_svm_kernel)
cat("AUC SVM Kernel:", auc_svm_kernel, "\n")


roc_svm_kernelpol <- roc(testData$primaryormetastasis, probabilities_svm_kernelpol[,2])
auc_svm_kernelpol <- auc(roc_svm_kernelpol)
cat("AUC SVM Kernel Polynomial:", auc_svm_kernelpol, "\n")

roc_dt <- roc(testData$primaryormetastasis, probabilities_dt[,2])
auc_dt <- auc(roc_dt)
cat("AUC Árbol de decisión:", auc_dt, "\n")

roc_lda <- roc(testData$primaryormetastasis, probabilities_lda$posterior[, 2])
auc_lda <- auc(roc_lda)
cat("AUC LDA:", auc_lda, "\n")

roc_rda <- roc(testData$primaryormetastasis, probabilities_rda$posterior[, 2])
auc_rda <- auc(roc_rda)
cat("AUC QDA:", auc_rda, "\n")





# Asegúrate de que la primera curva ROC sea la base
plot(roc_knn, col = "blue", main = "Curvas ROC", lwd = 2)
plot(roc_svm_linear, col = "red", add = TRUE, lwd = 2)
plot(roc_svm_kernel, col = "green", add = TRUE, lwd = 2)
plot(roc_svm_kernelpol, col = "orange", add = TRUE, lwd = 2)
plot(roc_dt, col = "purple", add = TRUE, lwd = 2)
plot(roc_lda, col = "pink", add = TRUE, lwd = 2)
plot(roc_rda, col = "yellow", add = TRUE, lwd = 2)


# Agregar leyenda
knn_legend <- paste("AUC k-NN:", round(auc_knn, 2))  # Redondeamos a 2 decimales, si es necesario
svm_legend <- paste("AUC SVM Lineal:", round(auc_svm_linear, 2))  # Redondeamos a 2 decimales, si es necesario
svmk_legend <- paste("AUC SVM Kernel:", round(auc_svm_kernel, 2))  # Redondeamos a 2 decimales, si es necesario
svmp_legend <- paste("AUC SVM Kernel Polynomial:", round(auc_svm_kernelpol, 2))  # Redondeamos a 2 decimales, si es necesario
dt_legend <- paste("AUC Decission Tree:", round(auc_dt, 2))  # Redondeamos a 2 decimales, si es necesario
lda_legend <- paste("AUC LDA:", round(auc_lda, 2))  # Redondeamos a 2 decimales, si es necesario
rda_legend <- paste("AUC RDA:", round(auc_rda, 2))  # Redondeamos a 2 decimales, si es necesario

legend("bottomright", legend = c(knn_legend, svm_legend, svmk_legend, svmp_legend, dt_legend, lda_legend, rda_legend),
       col = c("blue", "red", "green", "orange", "purple" ,"pink" ,"yellow" ), lwd = 2)




table(testData$primaryormetastasis)

# Calcular curvas PR para cada modelo
pr_knn <- pr.curve(scores.class0 = probabilities_knn[,2], weights.class0 = testData$primaryormetastasis == "Primary", curve = TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
pr_svm_linear <- pr.curve(scores.class0 = probabilities_svm_linear[,2], weights.class0 = testData$primaryormetastasis == "Primary", curve = TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
pr_svm_kernel <- pr.curve(scores.class0 = probabilities_svm_kernel[,2], weights.class0 = testData$primaryormetastasis == "Primary", curve = TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
pr_svm_kernelpol <- pr.curve(scores.class0 = probabilities_svm_kernelpol[,2], weights.class0 = testData$primaryormetastasis == "Primary", curve = TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
pr_dt <- pr.curve(scores.class0 = probabilities_dt[,2], weights.class0 = testData$primaryormetastasis == "Primary", curve = TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
pr_lda <- pr.curve(scores.class0 = probabilities_lda$posterior[, 2], weights.class0 = testData$primaryormetastasis == "Primary", curve = TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
valid_indices <- !is.na(probabilities_rda$posterior[, 2]) & !is.na(testData$primaryormetastasis)
pr_rda <- pr.curve(scores.class0 = probabilities_rda$posterior[valid_indices, 2], weights.class0 = testData$primaryormetastasis[valid_indices] == "Primary", curve = TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)

plot(pr_knn, col = "blue", lwd = 2, rand.plot = TRUE, fill.area = TRUE)
plot(pr_svm_linear, col = "red", add = TRUE, lwd = 2, rand.plot = TRUE, fill.area = TRUE)
plot(pr_svm_kernel, col = "green", add = TRUE, lwd = 2, rand.plot = TRUE, fill.area = TRUE)
plot(pr_svm_kernelpol, col = "orange", add = TRUE, lwd = 2, rand.plot = TRUE, fill.area = TRUE)
plot(pr_dt, col = "purple", add = TRUE, lwd = 2, rand.plot = TRUE, fill.area = TRUE)
plot(pr_lda, col = "pink", add = TRUE, lwd = 2, rand.plot = TRUE, fill.area = TRUE)
plot(pr_rda, col = "yellow", add = TRUE, lwd = 2, rand.plot = TRUE, fill.area = TRUE)


# Agregar leyenda
knn_legend <- paste("PR-Curve k-NN:", round(pr_knn$auc.integral, 2))  # Redondeamos a 2 decimales, si es necesario
svm_legend <- paste("PR-Curve SVM Lineal:", round(pr_svm_linear$auc.integral, 2))  # Redondeamos a 2 decimales, si es necesario
svmk_legend <- paste("PR-Curve SVM Kernel:", round(pr_svm_kernel$auc.integral, 2))  # Redondeamos a 2 decimales, si es necesario
svmp_legend <- paste("PR-Curve SVM Kernel Polynomial:", round(pr_svm_kernelpol$auc.integral, 2))  # Redondeamos a 2 decimales, si es necesario
dt_legend <- paste("PR-Curve Decission Tree:", round(pr_dt$auc.integral, 2))  # Redondeamos a 2 decimales, si es necesario
lda_legend <- paste("PR-Curve LDA:", round(pr_lda$auc.integral, 2))  # Redondeamos a 2 decimales, si es necesario
rda_legend <- paste("PR-Curve RDA:", round(pr_rda$auc.integral, 2))  # Redondeamos a 2 decimales, si es necesario

legend("bottomright", legend = c(knn_legend, svm_legend, svmk_legend, svmp_legend, dt_legend, lda_legend, rda_legend),
       col = c("blue", "red", "green", "orange", "purple" ,"pink" ,"yellow" ), lwd = 2)














######################## MODELOS DE REGRESION LINEAL

#### ---- Preparamos la base de datos para modelos de regresión ----

path2 <- '/Users/vic/Library/CloudStorage/GoogleDrive-vdelaopascual@gmail.com/Mi unidad/MU en Bioinformática (UNIR 2023)/Presentaciones Algoritmos e Inteligencia Artificial (AAAA-MM-DD)/Otras bases de datos/Framingham Dataset_dep.csv'
df2 <- read.csv(path2)

str(df2)
anyNA(df2)
df2 <- na.omit(df2)


# Eliminar 'code' y configurar las variables
df2 <- df2[, !names(df2) %in% c("code")]
df2 <- df2[df2$antihypertensive_use != "", ]
df2$sex <- as.factor(df2$sex)
df2$current_smoker <- as.factor(df2$current_smoker)
df2$diabetes_status <- as.factor(df2$diabetes_status)
df2$antihypertensive_use <- as.factor(df2$antihypertensive_use)
df2$prev_coronary_hd <- as.factor(df2$prev_coronary_hd)
df2$prev_angina <- as.factor(df2$prev_angina)
df2$prev_myocardial_infarction <- as.factor(df2$prev_myocardial_infarction)
df2$prev_stroke <- as.factor(df2$prev_stroke)
df2$prev_hypertension <- as.factor(df2$prev_hypertension)
df2$exam_cycle <- as.factor(df2$exam_cycle)
df2$death_status <- as.factor(df2$death_status)
df2$followup_angina <- as.factor(df2$followup_angina)
df2$hosp_myocardial_infarction <- as.factor(df2$hosp_myocardial_infarction)
df2$hosp_mi_or_fatal_chd <- as.factor(df2$hosp_mi_or_fatal_chd)
df2$followup_chd <- as.factor(df2$followup_chd)
df2$followup_stroke <- as.factor(df2$followup_stroke)
df2$followup_cvd <- as.factor(df2$followup_cvd)
df2$followup_hypertension <- as.factor(df2$followup_hypertension)
str(df2)
dim(df2)
df2 <- df2[, -c(18, 29:36)]


# Separar conjuntos de entrenamiento y prueba
set.seed(1995)
trainIndex <- createDataPartition(df2$total_cholesterol, p = 0.8, list = FALSE)
trainData_2 <- df2[trainIndex, ]
testData_2 <- df2[-trainIndex, ]




#### ---- kNN regression ----
knnModel_reg <- train(total_cholesterol ~.,
                      data = trainData_2,
                      method = "knn",
                      trControl = trainControl(method = "cv", number = 10),
                      preProcess = c("center", "scale"),
                      tuneLength = 30)
knnModel_reg
plot(knnModel_reg)

# Predicciones y métricas
knn_preds <- predict(knnModel_reg, newdata = testData_2)
postResample(knn_preds, testData_2$total_cholesterol)
knn_rmse <- sqrt(mean((knn_preds - testData_2$total_cholesterol)^2))  # RMSE para kNN

knn_preds_df <- data.frame(Real = testData_2$total_cholesterol,
                           Predicted = knn_preds)

knnModel_graph <- ggplot(knn_preds_df, aes(x = Real, y = Predicted)) +
  geom_point(alpha = 0.6) +  # Puntos de la dispersión
  geom_abline(slope = 1, intercept = 0, color = "red") +  # Línea de igualdad
  labs(title = paste("Modelo kNN - RMSE:", round(knn_rmse, 2)),
       x = "Valores Reales", 
       y = "Predicciones") +
  theme_minimal()

knnModel_graph



#### ---- SVM regression ----
svmModelLineal_reg <- train(total_cholesterol ~.,
                            data = trainData_2,
                            method = "svmLinear",
                            trControl = trainControl(method = "cv", number = 10),
                            preProcess = c("center", "scale"),
                            tuneGrid = expand.grid(C = seq(0, 2, length = 10))) #C grande lleva al sobreajuste, C pequeño al infraajuste
svmModelLineal_reg
plot(svmModelLineal_reg)

# Predicciones y métricas
svmModelLineal_preds <- predict(svmModelLineal_reg, newdata = testData_2)
postResample(svmModelLineal_preds, testData_2$total_cholesterol)
svm_rmse <- sqrt(mean((svmModelLineal_preds - testData_2$total_cholesterol)^2))  # RMSE para SVM

svmModelLineal_preds_df <- data.frame(Real = testData_2$total_cholesterol,
                                      Predicted = svmModelLineal_preds)

svmModelLineal_reg_graph <- ggplot(svmModelLineal_preds_df, aes(x = Real, y = Predicted)) +
  geom_point(alpha = 0.6) +  # Puntos de la dispersión
  geom_abline(slope = 1, intercept = 0, color = "red") +  # Línea de igualdad
  labs(title = paste("Modelo SVM lineal - RMSE:", round(svm_rmse, 2)),
       x = "Valores Reales", 
       y = "Predicciones") +
  theme_minimal()






svmModelRadial_reg <- train(total_cholesterol ~.,
                            data = trainData_2,
                            method = "svmRadial",
                            trControl = trainControl(method = "cv", number = 10),
                            preProcess = c("center", "scale"),
                            tuneLength = 10)
svmModelRadial_reg
plot(svmModelRadial_reg)

# Predicciones y métricas
svmModelRadial_preds <- predict(svmModelRadial_reg, newdata = testData_2)
postResample(svmModelRadial_preds, testData_2$total_cholesterol)
svm_kernel_rmse <- sqrt(mean((svmModelRadial_preds - testData_2$total_cholesterol)^2))  # RMSE para SVM Kernel

svmModelRadial_preds_df <- data.frame(Real = testData_2$total_cholesterol,
                                      Predicted = svmModelRadial_preds)

svmModelRadial_reg_graph <- ggplot(svmModelRadial_preds_df, aes(x = Real, y = Predicted)) +
  geom_point(alpha = 0.6) +  # Puntos de la dispersión
  geom_abline(slope = 1, intercept = 0, color = "red") +  # Línea de igualdad
  labs(title = paste("Modelo SVM Kernel - RMSE:", round(svm_kernel_rmse, 2)),
       x = "Valores Reales", 
       y = "Predicciones") +
  theme_minimal()



#### ---- Decission Tree regression ----
dtModel_reg <- train(total_cholesterol ~.,
                     data = trainData_2,
                     method = "rpart",
                     trControl = trainControl(method = "cv", number = 10),
                     preProcess = c("center", "scale"),
                     tuneGrid = expand.grid(cp = seq(0.0001, 0.5, by = 0.005)))  # Más valores para cp
dtModel_reg
plot(dtModel_reg)

# Predicciones y métricas
dtModel_preds <- predict(dtModel_reg, newdata = testData_2)
postResample(dtModel_preds, testData_2$total_cholesterol)
dt_rmse <- sqrt(mean((dtModel_preds - testData_2$total_cholesterol)^2))  # RMSE para Decision Tree

dtModel_preds_df <- data.frame(Real = testData_2$total_cholesterol,
                               Predicted = dtModel_preds)

dt_reg_graph <- ggplot(dtModel_preds_df, aes(x = Real, y = Predicted)) +
  geom_point(alpha = 0.6) +  # Puntos de la dispersión
  geom_abline(slope = 1, intercept = 0, color = "red") +  # Línea de igualdad
  labs(title = paste("Modelo Decision Tree - RMSE:", round(dt_rmse, 2)),
       x = "Valores Reales", 
       y = "Predicciones") +
  theme_minimal()


graphs_ml_reg <- list(knnModel_graph, svmModelLineal_reg_graph, svmModelRadial_reg_graph, dt_reg_graph)
grid.arrange(grobs = graphs_ml_reg, ncol = 2 )

# Importance plots
plot(varImp(knnModel_reg, scale = TRUE)) # importancia de las variables
plot(varImp(svmModelLineal_reg, scale = TRUE)) # importancia de las variables
plot(varImp(svmModelRadial_reg, scale = TRUE)) # importancia de las variables
plot(varImp(dtModel_reg, scale = TRUE)) # importancia de las variables


