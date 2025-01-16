rm(list=ls())

path <- "/Users/vic/Library/CloudStorage/GoogleDrive-vdelaopascual@gmail.com/Mi unidad/MU en Bioinformática (UNIR 2023)/Actividades (AAAA-MM-DD)/Actividad 1_junio2024 update"

setwd(path)

df <- read.csv("Dataset expresión genes.csv")




#### ---- X = GENES Y= TUMOR; training testing ----
library(dplyr)
library(caret)

df_genes_tumor <- df %>% dplyr::select(starts_with("AQ_"), tumor) # lo pongo dplyr::select para que lo detecte dplyr ya que luego con MASS tiene la misma función
df_genes_tumor$tumor <- as.factor(df_genes_tumor$tumor)

set.seed(1995)
train_index <- createDataPartition(df_genes_tumor$tumor, p = 0.8, list = FALSE)
train_index

training_tumor <- df_genes_tumor[train_index, ]
testing_tumor <- df_genes_tumor[-train_index, ]

table(training_tumor$tumor)
table(testing_tumor$tumor)

str(training_tumor)


numerical_columns <- training_tumor[, sapply(training_tumor, is.numeric)]
scaled_data <- scale(numerical_columns)
training_tumor <- cbind(scaled_data, tumor = training_tumor$tumor)

numerical_columns <- testing_tumor[, sapply(testing_tumor, is.numeric)]
scaled_data <- scale(numerical_columns)
testing_tumor <- cbind(scaled_data, tumor = testing_tumor$tumor)

training_tumor <- as.data.frame(training_tumor)
testing_tumor <- as.data.frame(testing_tumor)

# Crear la formula sumando cada gen
genes <- colnames(training_tumor[1:46])
formula_tumor <- as.formula(paste("tumor ~", paste(genes, collapse = "+")))
formula_tumor
training_tumor$tumor <- as.factor(training_tumor$tumor)



#### ---- X = GENES Y= TRAT; training testing ----

df_genes_trat <- df %>% dplyr::select(starts_with("AQ_"), trat) # lo pongo dplyr::select para que lo detecte dplyr ya que luego con MASS tiene la misma función
df_genes_trat$trat <- as.factor(df_genes_trat$trat)

set.seed(1995)
train_index <- createDataPartition(df_genes_trat$trat, p = 0.8, list = FALSE)
train_index

training_trat <- df_genes_trat[train_index, ]
testing_trat <- df_genes_trat[-train_index, ]

table(training_trat$trat)
table(testing_trat$trat)

str(training_trat)


numerical_columns <- training_trat[, sapply(training_trat, is.numeric)]
scaled_data <- scale(numerical_columns)
training_trat <- cbind(scaled_data, trat = training_trat$trat)

numerical_columns <- testing_trat[, sapply(testing_trat, is.numeric)]
scaled_data <- scale(numerical_columns)
testing_trat <- cbind(scaled_data, trat = testing_trat$trat)

training_trat <- as.data.frame(training_trat)
testing_trat <- as.data.frame(testing_trat)

# Crear la formula sumando cada gen
genes <- colnames(training_trat[1:10])
formula_trat <- as.formula(paste("trat ~", paste(genes, collapse = "+")))
formula_trat
training_trat$trat <- as.factor(training_trat$trat)


#### ---- LDA (lineal) ----
library(MASS)

# Ajustar el modelo LDA en el entrenamiento
lda_model <- lda(formula_tumor, data = training_tumor)
lda_model$scaling # contribuciones/coeficientes

# LD1 y LD2 son las primeras dos funciones discriminantes generadas por el LDA. 
# --> LD1 es la combinación lineal de las variables que maximiza la separación entre las clases.
# --> LD2 es la segunda combinación lineal de las variables, que también busca separar las clases, pero es ortogonal (es decir, no correlacionada) a LD1.

# Si un coeficiente es positivo, significa que un aumento en esa variable aumentará el valor de la función discriminante.
# Si un coeficiente es negativo, significa que un aumento en esa variable disminuirá el valor de la función discriminante.

# AQ_ALOX5 tiene un coeficiente muy grande en LD1 (11.298), lo que significa que AQ_ALOX5 es una variable muy importante para separar las clases en la primera función discriminante.
# AQ_PPARG tiene un coeficiente negativo significativo en LD1 (-12.695), lo que indica que un aumento en AQ_PPARG tiende a desplazar la clasificación.

# En términos prácticos, las variables con coeficientes más grandes (en valor absoluto) tienen más poder para distinguir entre las clases. Por ejemplo:
# Variables como AQ_ALOX5 y AQ_PPARG (con coeficientes grandes) pueden estar desempeñando un papel importante en cómo el modelo clasifica las observaciones.


lda_pred <- predict(lda_model, newdata = training_tumor)
lda_pred$x

# Realizar predicciones sobre el conjunto de prueba
lda_predictions <- predict(lda_model, newdata = testing_tumor)
lda_predictions$x

# Obtener la predicción (predicciones de la clase)
predicted_classes <- lda_predictions$class
predicted_classes
length(predicted_classes)

# Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
true_classes <- as.factor(testing_tumor$tumor)
true_classes
length(true_classes)

# Crear la matriz de confusión (tumor predicho testing vs. tumor real testing)
confusion <- confusionMatrix(predicted_classes, true_classes)
print(confusion)

# Acuracy: El modelo tiene una precisión del 72.73%. Esto significa que el 72.73% de las predicciones fueron correctas (el modelo acertó en la clasificación).
# Sensitivity: La sensibilidad mide cuántas veces el modelo predijo correctamente la clase 1 entre las veces que realmente era clase 1. En este caso, es 1.0000, lo que significa que el modelo detectó todas las instancias de clase 1 correctamente (100%).
# Specificity: La especificidad mide cuántas veces el modelo predijo correctamente que no era clase 1. En este caso, 83.33% de las veces que no era clase 1, el modelo lo predijo correctamente.
# Pos Pred Value: Este valor indica el porcentaje de veces que el modelo predijo clase 1 correctamente cuando hizo esa predicción. Aquí es 83.33%.
# Neg Pred Value: Mide cuántas veces el modelo predijo correctamente que no era clase 1 cuando realmente no era clase 1. En este caso, el modelo predijo correctamente 100% de las veces que no era clase 1.
# Balanced Accuracy: es un promedio de la sensibilidad y especificidad. En este caso, es 91.67%, lo que sugiere que el modelo tiene un buen desempeño en la clase 1.





# grafico
predict(lda_model)

library(ggplot2)
lda.data <- cbind(training_tumor, predict(lda_model)$x)
ggplot(lda.data, aes(LD1, LD2)) + 
  geom_point(aes(color = tumor)) +
  theme_classic()


head(lda.data[,47:49])
# Observaciones con valores altos en LD1 están probablemente asociadas con una clase específica (en este caso, parece estar relacionado con tumor = 3 porque la observación 7 tiene un valor muy alto en LD1).
# Observaciones con valores altos en LD2 (como en las observaciones 1, 4 y 8) parecen estar asociadas con la clase tumor = 1.





#### ---- QDA (no lineal) ----

# Ajustar el modelo QDA en el entrenamiento
qda_model <- qda(formula_trat, data = training_trat)
qda_model


# Realizar predicciones sobre el conjunto de prueba
qda_predictions <- predict(qda_model, newdata = testing_trat)

# Obtener la predicción (predicciones de la clase)
predicted_classes <- qda_predictions$class
predicted_classes
length(predicted_classes)

# Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
true_classes <- as.factor(testing_trat$trat)
true_classes
length(true_classes)

# Crear la matriz de confusión (trat predicho testing vs. trat real testing)
confusion <- confusionMatrix(predicted_classes, true_classes)
print(confusion)

# Acuracy: El modelo tiene una precisión del 50%. Esto significa que el 50% de las predicciones fueron correctas (el modelo acertó en la clasificación).
# Sensitivity: La sensibilidad mide cuántas veces el modelo predijo correctamente la clase 1 entre las veces que realmente era clase 1. En este caso es 0.0, lo que significa que el modelo no detectó ninguna instancia de clase 1 correctamente (ningún verdadero positivo para la clase 1).
# Specificity: La especificidad mide cuántas veces el modelo predijo correctamente que no era clase 1. Aquí es 1.0, lo que significa que el modelo predijo correctamente todas las veces que no era clase 1 (es decir, todas las veces que realmente era clase 2).
# Pos Pred Value: Este valor indica el porcentaje de veces que el modelo predijo clase 1 correctamente cuando hizo esa predicción. Dado que no hay verdaderos positivos (TP = 0), este valor es NaN (no disponible).
# Neg Pred Value: Mide cuántas veces el modelo predijo correctamente que no era clase 1 cuando realmente no era clase 1. En este caso es 0.5, lo que significa que el modelo predijo correctamente 50% de las veces que no era clase 1.



#### ---- RDA (regularizado) ----
library(klaR)
formula_tumor

# Ajustar el modelo RDA en el entrenamiento
rda_model <- rda(formula_tumor, data = training_tumor)
rda_model

rda_pred <- predict(rda_model, newdata = training_tumor)


# Realizar predicciones sobre el conjunto de prueba
rda_predictions <- predict(rda_model, newdata = testing_tumor)

# Obtener la predicción (predicciones de la clase)
predicted_classes <- rda_predictions$class
predicted_classes
length(predicted_classes)

# Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
true_classes <- as.factor(testing_tumor$tumor)
true_classes
length(true_classes)


# Crear la matriz de confusión (trat predicho testing vs. trat real testing)
confusion <- confusionMatrix(predicted_classes, true_classes)
print(confusion)

# Acuracy: La precisión global del modelo es 72.73%, lo que significa que el modelo acertó el 72.73% de las predicciones.
# Sensitivity: La sensibilidad para la clase 1 es 1.0000, lo que significa que el modelo detectó correctamente todas las instancias de clase 1 (no hubo falsos negativos para clase 1).
# Specificity: La especificidad para la clase 1 es 0.6667, lo que significa que el modelo predijo correctamente 66.67% de las veces que no era clase 1 (cuando realmente era clase 2 o 3).
# Pos Pred Value: El valor predictivo positivo para la clase 1 es 0.7143, lo que indica que, cuando el modelo predice que es clase 1, hay un 71.43% de probabilidad de que sea realmente clase 1.
# Neg Pred Value: El valor predictivo negativo para la clase 1 es 1.0000, lo que significa que cuando el modelo predice que no es clase 1, la predicción es correcta el 100% de las veces.





#### ---- FDA (para datos no paramétricos, tarda mucho en correr...) ----
library(mda)

# Ajustar el modelo RDA en el entrenamiento
fda_model <- fda(formula_tumor, data = training_tumor, method = bruto)
fda_model

fda_pred <- predict(fda_model, newdata = training_tumor)


# Realizar predicciones sobre el conjunto de prueba
fda_predictions <- predict(fda_model, newdata = testing_tumor)

# Obtener la predicción (predicciones de la clase)
predicted_classes <- rda_predictions$class
predicted_classes
length(predicted_classes)

# Obtener las verdaderas etiquetas (las clases reales en el conjunto de prueba)
true_classes <- as.factor(testing_tumor$tumor)
true_classes
length(true_classes)


# Crear la matriz de confusión (trat predicho testing vs. trat real testing)
confusion <- confusionMatrix(predicted_classes, true_classes)
print(confusion)



























# simulamos más datos (no importa para los bioinfos)
library(synthpop)
library(dplyr)

## simulacion nº1
set.seed(1995) 
df_sintetico <- syn(df)  # Crea un objeto con los datos sintéticos
df_sintetico <- df_sintetico$syn # Obtener el dataframe de pacientes simulados

# Verifica los datos
str(df_sintetico)
str(df)
table(df_sintetico$tumor)
max(df$id)

df_sintetico$id <- 69:(69 + nrow(df_sintetico) - 1)
all(names(df) == names(df_sintetico))  # asegurar de que ambas bases tengan las mismas columnas = Esto debe devolver TRUE
df <- bind_rows(df, df_sintetico) # combinamos
any(duplicated(df$id))  # duplicados en id = debe ser FALSE



## simulacion nº2
set.seed(1995) 
df_sintetico <- syn(df)  # Crea un objeto con los datos sintéticos
df_sintetico <- df_sintetico$syn # Obtener el dataframe de pacientes simulados

# Verifica los datos
str(df_sintetico)
str(df)
table(df_sintetico$tumor)
max(df$id)

df_sintetico$id <- 134:(134 + nrow(df_sintetico) - 1)
all(names(df) == names(df_sintetico))  # asegurar de que ambas bases tengan las mismas columnas = Esto debe devolver TRUE
df <- bind_rows(df, df_sintetico) # combinamos
any(duplicated(df$id))  # duplicados en id = debe ser FALSE


## simulacion nº3
set.seed(1995) 
df_sintetico <- syn(df)  # Crea un objeto con los datos sintéticos
df_sintetico <- df_sintetico$syn # Obtener el dataframe de pacientes simulados

# Verifica los datos
str(df_sintetico)
str(df)
table(df_sintetico$tumor)
max(df$id)

df_sintetico$id <- 264:(264 + nrow(df_sintetico) - 1)
all(names(df) == names(df_sintetico))  # asegurar de que ambas bases tengan las mismas columnas = Esto debe devolver TRUE
df <- bind_rows(df, df_sintetico) # combinamos
any(duplicated(df$id))  # duplicados en id = debe ser FALSE




