# ============================================================================
# MODELAGEM DE DISTRIBUIÇÃO DE ESPÉCIES - Terminalia acuminata
# ============================================================================
# Baseado em: Fernandes et al. (2023) e Fernandes & Rego (2011)
# Data: 2025
# Versão: 30.0 
# Autoria - Antônio Lucas Barreira Rodrigues
# ============================================================================

# 1. CARREGAMENTO DE PACOTES --------------------------------------------------

# Lista de pacotes necessários
packages <- c("raster", "rgdal", "sp", "dismo", "vegan", 
              "ggplot2", "maps","rJava", "MASS", 
              "mgcv", "randomForest", "geodata", "terra",
              "corrplot", "pROC", "RColorBrewer", "sf")

# Função para instalar e carregar pacotes
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

install_and_load(packages)

# 2. CONFIGURAÇÃO INICIAL -----------------------------------------------------

# Definir diretório de trabalho
setwd("D:/Green_status/terminalia_aucuminata/analise")

# Criar pastas para organização
dir.create("dados_taucuminata", showWarnings = FALSE)
dir.create("resultados", showWarnings = FALSE)
dir.create("mapas", showWarnings = FALSE)

# 3. CARREGAMENTO DOS DADOS DE OCORRÊNCIA ------------------------------------

# Carregar dados do arquivo CSV
ocorrencias <- read.csv("D:/Green_status/terminalia_aucuminata/pontos_terminalia.csv", 
                        stringsAsFactors = FALSE)

# Verificar estrutura dos dados
str(ocorrencias)
head(ocorrencias)
cat("Número inicial de pontos:", nrow(ocorrencias), "\n")

# Remover duplicatas espaciais
#ocorrencias <- ocorrencias[!duplicated(ocorrencias[, c("longitude", "latitude")]), ]
#cat("Número de pontos após remover duplicatas:", nrow(ocorrencias), "\n")

# Criar objeto espacial
coordinates(ocorrencias) <- ~longitude + latitude
crs(ocorrencias) <- CRS("+proj=longlat +datum=WGS84")

# Visualizar pontos
plot(ocorrencias, pch = 19, col = "red", main = "Pontos de ocorrência - Terminalia acuminata")
map("world", add = TRUE)

# 4. DEFINIÇÃO DA ÁREA DE ESTUDO ---------------------------------------------

# Função para definir área de estudo
definir_area_estudo <- function(metodo = c("buffer", "estado_rj"), 
                                pontos_ocorrencia = ocorrencias,
                                tamanho_buffer = 2,
                                buffer_adicional_rj = 0.1) {
  
  metodo <- match.arg(metodo)
  
  if(metodo == "buffer") {
    cat("Definindo área de estudo por BUFFER\n")
    cat("Tamanho do buffer:", tamanho_buffer, "graus\n")
    
    # Criar extensão com buffer
    study_extent <- extent(pontos_ocorrencia) + tamanho_buffer
    
    # Criar polígono
    study_area <- as(study_extent, "SpatialPolygons")
    crs(study_area) <- CRS("+proj=longlat +datum=WGS84")
    
  } else {
    cat("Definindo área de estudo por LIMITES DO ESTADO RJ\n")
    
    # Verificar se geobr está instalado
    if(!require("geobr")) {
      install.packages("geobr")
      library(geobr)
    }
    
    # Baixar/carregar RJ
    rj_sf <- read_state(code_state = "RJ", year = 2020, showProgress = FALSE)
    rj_sp <- as(rj_sf, "Spatial")
    rj_sp <- spTransform(rj_sp, CRS("+proj=longlat +datum=WGS84"))
    
    # Definir extensão
    study_extent <- extent(rj_sp) + buffer_adicional_rj
    study_area <- rj_sp
    
    cat("Buffer adicional:", buffer_adicional_rj, "graus\n")
  }
  
  # Retornar lista com ambos objetos
  return(list(extent = study_extent, area = study_area))
}

# ESCOLHA O MÉTODO AQUI:
resultado <- definir_area_estudo(metodo = "buffer")  # ou "buffer" ou "estado_rj"

# Atribuir aos objetos originais
study_extent <- resultado$extent
study_area <- resultado$area

# Visualizar
plot(study_area, main = "Área de estudo selecionada")
points(ocorrencias, pch = 19, col = "red", cex = 0.8)

# 5. DOWNLOAD E PREPARAÇÃO DAS VARIÁVEIS AMBIENTAIS --------------------------

if (!require("geodata")) {
  install.packages("geodata")
  library(geodata)
}

# Alternativa 1: Usando geodata (se disponível)
tryCatch({
  bioclim_vars <- geodata::worldclim_global(var = "bio", res = 0.5, path = "dados/worldclim")
}, error = function(e) {
  cat("Erro com geodata::worldclim_global. Tentando método alternativo...\n")
})

# Download das variáveis bioclimáticas do WorldClim
bioclim_vars <- worldclim_global(var = "bio", res = 0.5, path = "dados/worldclim")

# Recortar para área de estudo
bioclim_crop <- crop(bioclim_vars, study_extent)

# Selecionar variáveis ambientais (conforme o estudo da Dimorphandra)
selected_indices <- c(1: 19)
env_vars <- bioclim_crop[[selected_indices]]

# Renomear variáveis para facilitar
names(env_vars) <- paste0("bio", 1:19)

# Definir paleta de cores verde-vermelho
green_red <- colorRampPalette(c("darkgreen", "green", "yellow", "orange", "red"))

# Visualizar variáveis com nova paleta
plot(env_vars, col = green_red(100), main = names(env_vars))



# 6. ANÁLISE DE CORRELAÇÃO E SELEÇÃO DE VARIÁVEIS ----------------------------

# Converter pontos para SpatVector
ocorrencias_terra <- terra::vect(ocorrencias)

# Extrair valores das variáveis ambientais nos pontos
env_values <- terra::extract(env_vars, ocorrencias_terra)

# Remover a coluna ID e NAs
env_values <- env_values[, -1]  # Remove coluna ID
env_values <- na.omit(env_values)

# Análise de correlação
cor_matrix <- cor(env_values)
print(round(cor_matrix, 2))

# Instalar e carregar corrplot se necessário
if (!require("corrplot")) {
  install.packages("corrplot")
  library(corrplot)
}

# Visualizar correlação (três métodos de visualização)
corrplot(cor_matrix, method = "circle", type = "upper", 
         order = "hclust", tl.col = "black", tl.srt = 45)

# Método alternativo: Usar ggplot2
    # Converter matriz para formato longo
cor_melted <- reshape2::melt(cor_matrix)
    

#CALCULO DE FIV E SELEÇÃO DE VARIÁVEIS

# Função customizada para calcular VIF
calculate_vif <- function(env_raster) {
  # Amostrar pontos aleatórios do raster
  set.seed(123)
  sample_points <- spatSample(env_raster, size = 5000, method = "random", 
                              na.rm = TRUE, as.points = TRUE)
  sample_df <- as.data.frame(values(sample_points))
  
  # Remover coluna de geometria se existir
  sample_df <- sample_df[, !names(sample_df) %in% c("x", "y", "geometry")]
  
  # Calcular VIF
  vif_values <- numeric(ncol(sample_df))
  names(vif_values) <- names(sample_df)
  
  for (i in 1:ncol(sample_df)) {
    formula <- as.formula(paste(names(sample_df)[i], "~ ."))
    model <- lm(formula, data = sample_df)
    r_squared <- summary(model)$r.squared
    vif_values[i] <- 1 / (1 - r_squared)
}
  
 return(vif_values)
}

# Método alternativo: Seleção iterativa de variáveis
# Remove uma variável por vez, começando pela de maior VIF
vif_iterative_selection <- function(env_raster, vif_threshold = 10) {
  current_vars <- names(env_raster)
  removed_vars <- character(0)
  
  repeat {
    # Calcular VIF para variáveis atuais
    current_env <- env_raster[[current_vars]]
    vif_current <- calculate_vif(current_env)
    
    # Verificar se alguma variável excede o threshold
    max_vif <- max(vif_current)
    if (max_vif <= vif_threshold) {
      break  # Todas as variáveis estão dentro do critério
    }
    
    # Remover a variável com maior VIF
    var_to_remove <- names(vif_current)[which.max(vif_current)]
    removed_vars <- c(removed_vars, var_to_remove)
    current_vars <- current_vars[current_vars != var_to_remove]
    
    cat("Removendo", var_to_remove, "- VIF =", round(max_vif, 2), "\n")
    
    # Verificação de segurança
    if (length(current_vars) < 2) {
      cat("ATENÇÃO: Menos de 2 variáveis restantes. Parando seleção.\n")
      break
    }
  }
  
  return(list(
    selected_vars = current_vars,
    removed_vars = removed_vars,
    final_vif = calculate_vif(env_raster[[current_vars]])
  ))
}

# Aplicar seleção iterativa

iterative_result <- vif_iterative_selection(env_vars, vif_threshold = 10)

cat("\nVariáveis finais selecionadas:", paste(iterative_result$selected_vars, collapse = ", "), "\n")
cat("VIF final das variáveis selecionadas:\n")
print(round(iterative_result$final_vif, 2))

# Criar raster final com variáveis selecionadas
env_vars_final <- env_vars[[iterative_result$selected_vars]]


# 7. GERAÇÃO DE PONTOS DE BACKGROUND/PSEUDO-AUSÊNCIAS ------------------------

# Número de pontos de background (10x o número de presenças)
n_bg <- nrow(ocorrencias) * 10

# Criar máscara para garantir que pontos sejam gerados apenas onde há dados
mask <- env_vars_final[[1]]
mask[!is.na(mask)] <- 1

# Converter para RasterLayer para usar randomPoints
mask_raster <- raster(mask)

# Gerar pontos aleatórios apenas onde há dados no raster
set.seed(123)
bg_points <- randomPoints(mask_raster, n = n_bg, p = coordinates(ocorrencias))

# Converter para dataframe
bg_points_df <- as.data.frame(bg_points)
names(bg_points_df) <- c("longitude", "latitude")

# Visualizar
plot(env_vars_final[[1]], col = green_red(100), 
     main = "Pontos de presença (preto) e background (azul)")
points(ocorrencias, pch = 19, col = "black", cex = 0.8)
points(bg_points_df, pch = 1, col = "blue", cex = 0.8)
legend("topright", legend = c("Presença", "Background"), 
       col = c("black", "blue"), pch = c(19, 1))


# 8. PREPARAÇÃO DOS DADOS PARA MODELAGEM -------------------------------------

# Extrair valores ambientais para presenças
pres_env <- terra::extract(env_vars_final, ocorrencias_terra)
pres_env <- pres_env[, -1]  # Remove ID
pres_data <- cbind(presence = 1, coordinates(ocorrencias), pres_env)

# Extrair valores ambientais para background
bg_points_terra <- vect(bg_points_df, geom = c("longitude", "latitude"), 
                        crs = "+proj=longlat +datum=WGS84")
bg_env <- terra::extract(env_vars_final, bg_points_terra)
bg_env <- bg_env[, -1]  # Remove ID
bg_data <- cbind(presence = 0, bg_points_df, bg_env)

# Combinar dados
model_data <- rbind(pres_data, bg_data)
model_data <- na.omit(model_data)

cat("\nDados para modelagem:")
cat("\n- Presenças:", sum(model_data$presence == 1))
cat("\n- Background:", sum(model_data$presence == 0))
cat("\n- Total:", nrow(model_data), "\n")

# 9. DIVISÃO DOS DADOS (TREINO/TESTE) ----------------------------------------

# 70% treino, 30% teste
set.seed(123)
train_index <- sample(1:nrow(model_data), size = 0.7 * nrow(model_data))
train_data <- model_data[train_index, ]
test_data <- model_data[-train_index, ]

# Separar coordenadas dos dados ambientais
train_coords <- train_data[, c("longitude", "latitude")]
test_coords <- test_data[, c("longitude", "latitude")]

# Remover coordenadas dos dados de treino/teste para modelagem
train_data_env <- train_data[, !(names(train_data) %in% c("longitude", "latitude"))]
test_data_env <- test_data[, !(names(test_data) %in% c("longitude", "latitude"))]

cat("\nDivisão dos dados:")
cat("\n- Treino:", nrow(train_data), "registros")
cat("\n- Teste:", nrow(test_data), "registros\n")

# 10. MODELAGEM COM DIFERENTES ALGORITMOS ------------------------------------

# 10.1 MAXENT ----------------------------------------------------------------

# Criar diretório para MaxEnt com caminho absoluto
maxent_dir <- file.path(getwd(), "resultados", "maxent")
dir.create(maxent_dir, recursive = TRUE, showWarnings = FALSE)

# Preparar dados para MaxEnt
# Pontos de presença do conjunto de treino
pres_train_idx <- which(train_data$presence == 1)
pres_train_coords <- train_coords[pres_train_idx, ]

# Background points do conjunto de treino
bg_train_idx <- which(train_data$presence == 0)
bg_train_coords <- train_coords[bg_train_idx, ]

# Converter SpatRaster para RasterStack
env_vars_raster <- raster::stack(env_vars_final)

# Verificar dados antes de executar
cat("\nVerificando dados para MaxEnt:\n")
cat("Pontos de presença:", nrow(pres_train_coords), "\n")
cat("Pontos de background:", nrow(bg_train_coords), "\n")
cat("Número de variáveis:", nlayers(env_vars_raster), "\n")
cat("Diretório de saída:", maxent_dir, "\n")

# Executar MaxEnt 
cat("\nExecutando MaxEnt...\n")
maxent_model <- dismo::maxent(
  x = env_vars_raster,
  p = pres_train_coords,
  a = bg_train_coords,
  path = maxent_dir,
  args = c("responsecurves=true", 
           "jackknife=true",
           "outputformat=logistic",
           "writeplotdata=false",
           "pictures=false",
           "askoverwrite=false",
           "maximumiterations=5000")
)

# Predição

maxent_pred_raster <- predict(maxent_model, env_vars_raster)

# Converter para SpatRaster
maxent_pred <- rast(maxent_pred_raster)

# Predição para conjunto de teste
# Extrair valores ambientais para pontos de teste usando terra
test_points <- vect(test_coords, geom = c("longitude", "latitude"), crs = crs(env_vars_final))
maxent_pred_test_values <- terra::extract(maxent_pred, test_points)
maxent_pred_test <- maxent_pred_test_values[,2]

# Verificar e corrigir valores NA se existirem
if(any(is.na(maxent_pred_test))) {
  cat("AVISO: Removendo", sum(is.na(maxent_pred_test)), "valores NA das predições\n")
  maxent_pred_test[is.na(maxent_pred_test)] <- 0
}

# Visualizar
plot(maxent_pred, col = green_red(100), 
     main = "MaxEnt - Distribuição potencial de Terminalia aucuminata")
points(ocorrencias, pch = 19, cex = 0.5)

# 10.2 GLM (Generalized Linear Model) ------------------------------------------

# Preparar dados ambientais (removendo coordenadas)
train_env <- train_data_env[, !names(train_data_env) %in% "presence"]
test_env <- test_data_env[, !names(test_data_env) %in% "presence"]

# Preparar fórmula
formula_glm <- as.formula(paste("presence ~", 
                                paste(names(train_env), collapse = " + ")))

# Mostrar fórmula
cat("Fórmula GLM:", deparse(formula_glm), "\n")

# Ajustar modelo

glm_model <- glm(formula_glm, data = train_data_env, 
                 family = binomial(link = "logit"))

# Resumo do modelo
summary(glm_model)

# Predição no conjunto de teste
glm_pred_test <- predict(glm_model, newdata = test_data_env, type = "response")

# Predição espacial
glm_pred_raster <- predict(env_vars_final, glm_model, type = "response")

# Visualizar
plot(glm_pred_raster, col = green_red(100), 
     main = "GLM - Distribuição potencial de Terminalia aucuminata")
points(ocorrencias, pch = 19, cex = 0.5)

# 10.3 Random Forest ---------------------------------------------------------

if (!require("randomForest")) {
  install.packages("randomForest")
  library(randomForest)
}



# Converter presence para factor (fazer cópia para não afetar GLM)
train_data_rf <- train_data_env
test_data_rf <- test_data_env
train_data_rf$presence <- as.factor(train_data_rf$presence)
test_data_rf$presence <- as.factor(test_data_rf$presence)

# Ajustar modelo
rf_model <- randomForest(formula_glm, data = train_data_rf, 
                         ntree = 5000, importance = TRUE)

# Resumo do modelo
print(rf_model)

# Importância das variáveis
importance(rf_model)
varImpPlot(rf_model)

# Predição no conjunto de teste
rf_pred_test <- predict(rf_model, newdata = test_data_rf, type = "prob")[,2]

# Predição espacial
rf_pred_raster <- predict(env_vars_final, rf_model, type = "prob", index = 2)

# Visualizar
plot(rf_pred_raster, col = green_red(100), 
     main = "Random Forest - Distribuição potencial de Terminalia aucuminata")
points(ocorrencias, pch = 19, cex = 0.5)

# 10.4 GAM (Generalized Additive Model) --------------------------------------

if (!require("mgcv")) {
  install.packages("mgcv")
  library(mgcv)
}



# Criar fórmula para GAM com smooth terms
# Usar k=3 ou k=4 para evitar overfitting com poucos dados
formula_gam <- as.formula(paste("presence ~ ", 
                                paste("s(", names(train_env), ", k=3)", 
                                      collapse = " + ")))

# Mostrar fórmula
cat("Fórmula GAM:", deparse(formula_gam), "\n")

# Ajustar modelo GAM
gam_model <- gam(formula_gam, 
                 data = train_data_env, 
                 family = binomial(link = "logit"),
                 method = "REML")

# Resumo do modelo
summary(gam_model)

# Verificar diagnósticos
gam.check(gam_model)

# Plotar smooth terms
plot(gam_model, pages = 1, residuals = TRUE, pch = 1, cex = 0.5)

# Predição no conjunto de teste
gam_pred_test <- predict(gam_model, newdata = test_data_env, type = "response")

# Predição espacial
gam_pred_raster <- predict(env_vars_final, gam_model, type = "response")

# Visualizar
plot(gam_pred_raster, col = green_red(100), 
     main = "GAM - Distribuição potencial de Terminalia aucuminata")
points(ocorrencias, pch = 19, cex = 0.5)

# 11. AVALIAÇÃO COMPLETA DOS MODELOS -----------------------------------------

# Função para calcular métricas incluindo TSS
calculate_metrics <- function(observed, predicted, threshold = 0.5) {
  # Converter para numérico se necessário
  if(is.factor(observed)) observed <- as.numeric(as.character(observed))
  
  pred_binary <- ifelse(predicted >= threshold, 1, 0)
  
  # Matriz de confusão
  conf_matrix <- table(observed, pred_binary)
  
  # Garantir que a matriz seja 2x2
  if(nrow(conf_matrix) == 1 || ncol(conf_matrix) == 1) {
    # Adicionar linha ou coluna faltante
    if(nrow(conf_matrix) == 1) {
      if(rownames(conf_matrix)[1] == "0") {
        conf_matrix <- rbind(conf_matrix, c(0, 0))
        rownames(conf_matrix)[2] <- "1"
      } else {
        conf_matrix <- rbind(c(0, 0), conf_matrix)
        rownames(conf_matrix)[1] <- "0"
      }
    }
    if(ncol(conf_matrix) == 1) {
      if(colnames(conf_matrix)[1] == "0") {
        conf_matrix <- cbind(conf_matrix, c(0, 0))
        colnames(conf_matrix)[2] <- "1"
      } else {
        conf_matrix <- cbind(c(0, 0), conf_matrix)
        colnames(conf_matrix)[1] <- "0"
      }
    }
  }
  
  print("Matriz de Confusão:")
  print(conf_matrix)
  
  # Métricas
  TP <- conf_matrix[2,2]
  TN <- conf_matrix[1,1] 
  FP <- conf_matrix[1,2]
  FN <- conf_matrix[2,1]
  
  accuracy <- (TP + TN) / sum(conf_matrix)
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  
  # TSS (True Skill Statistic)
  tss <- sensitivity + specificity - 1
  
  # Kappa
  po <- accuracy
  pe <- ((TP + FN) * (TP + FP) + (TN + FP) * (TN + FN)) / sum(conf_matrix)^2
  kappa <- (po - pe) / (1 - pe)
  
  return(data.frame(
    Accuracy = round(accuracy, 3),
    Sensitivity = round(sensitivity, 3),
    Specificity = round(specificity, 3),
    TSS = round(tss, 3),
    Kappa = round(kappa, 3)
  ))
}

# Função para encontrar threshold 
find_optimal_threshold <- function(observed, predicted) {
  if(is.factor(observed)) observed <- as.numeric(as.character(observed))
  
  thresholds <- seq(0, 1, by = 0.01)
  tss_values <- numeric(length(thresholds))
  
  for(i in 1:length(thresholds)) {
    pred_binary <- ifelse(predicted >= thresholds[i], 1, 0)
    conf_matrix <- table(observed, pred_binary)
    
    if(nrow(conf_matrix) == 2 && ncol(conf_matrix) == 2) {
      TP <- conf_matrix[2,2]
      TN <- conf_matrix[1,1] 
      FP <- conf_matrix[1,2]
      FN <- conf_matrix[2,1]
      
      sensitivity <- TP / (TP + FN)
      specificity <- TN / (TN + FP)
      tss_values[i] <- sensitivity + specificity - 1
    } else {
      tss_values[i] <- NA
    }
  }
  
  optimal_idx <- which.max(tss_values)
  return(list(threshold = thresholds[optimal_idx], 
              tss = tss_values[optimal_idx]))
}

# Calcular AUC
if (!require("pROC")) {
  install.packages("pROC")
  library(pROC)
}

# Preparar dados de observação para teste
test_obs <- as.numeric(as.character(test_data_env$presence))

# Avaliação MaxEnt

maxent_optimal <- find_optimal_threshold(test_obs, maxent_pred_test)
cat("Threshold:", maxent_optimal$threshold, "\n")
maxent_auc <- roc(test_obs, maxent_pred_test)
cat("AUC MaxEnt:", round(maxent_auc$auc, 3), "\n")
maxent_metrics <- calculate_metrics(test_obs, maxent_pred_test, 
                                    threshold = maxent_optimal$threshold)
print(maxent_metrics)

# Avaliação GLM

glm_optimal <- find_optimal_threshold(test_obs, glm_pred_test)
cat("Threshold:", glm_optimal$threshold, "\n")
glm_auc <- roc(test_obs, glm_pred_test)
cat("AUC GLM:", round(glm_auc$auc, 3), "\n")
glm_metrics <- calculate_metrics(test_obs, glm_pred_test, 
                                 threshold = glm_optimal$threshold)
print(glm_metrics)

# Avaliação Random Forest

rf_optimal <- find_optimal_threshold(test_obs, rf_pred_test)
cat("Threshold:", rf_optimal$threshold, "\n")
rf_auc <- roc(test_obs, rf_pred_test)
cat("AUC Random Forest:", round(rf_auc$auc, 3), "\n")
rf_metrics <- calculate_metrics(test_obs, rf_pred_test, 
                                threshold = rf_optimal$threshold)
print(rf_metrics)

# Avaliação GAM

gam_optimal <- find_optimal_threshold(test_obs, gam_pred_test)
cat("Threshold:", gam_optimal$threshold, "\n")
gam_auc <- roc(test_obs, gam_pred_test)
cat("AUC GAM:", round(gam_auc$auc, 3), "\n")
gam_metrics <- calculate_metrics(test_obs, gam_pred_test, 
                                 threshold = gam_optimal$threshold)
print(gam_metrics)

# 12. ENSEMBLE MODEL (COM OS 4 MODELOS) --------------------------------------

# Média ponderada das predições (pode ajustar pesos conforme performance)
# Usar AUC como peso
total_auc <- as.numeric(maxent_auc$auc) + as.numeric(glm_auc$auc) + 
  as.numeric(rf_auc$auc) + as.numeric(gam_auc$auc)
w_maxent <- as.numeric(maxent_auc$auc) / total_auc
w_glm <- as.numeric(glm_auc$auc) / total_auc
w_rf <- as.numeric(rf_auc$auc) / total_auc
w_gam <- as.numeric(gam_auc$auc) / total_auc


cat("\nPesos do Ensemble:\n")
cat("MaxEnt:", round(w_maxent, 3), "\n")
cat("GLM:", round(w_glm, 3), "\n")
cat("Random Forest:", round(w_rf, 3), "\n")
cat("GAM:", round(w_gam, 3), "\n")

# Ensemble ponderado
ensemble_pred_test <- (w_maxent * maxent_pred_test + 
                         w_glm * glm_pred_test + 
                         w_rf * rf_pred_test +
                         w_gam * gam_pred_test)

ensemble_pred_raster <- (w_maxent * maxent_pred + 
                           w_glm * glm_pred_raster + 
                           w_rf * rf_pred_raster +
                           w_gam * gam_pred_raster)

# Avaliação do ensemble

ensemble_optimal <- find_optimal_threshold(test_obs, ensemble_pred_test)
cat("Threshold:", ensemble_optimal$threshold, "\n")
ensemble_auc <- roc(test_obs, ensemble_pred_test)
cat("AUC Ensemble:", round(ensemble_auc$auc, 3), "\n")
ensemble_metrics <- calculate_metrics(test_obs, ensemble_pred_test, 
                                      threshold = ensemble_optimal$threshold)
print(ensemble_metrics)

# Visualizar ensemble
plot(ensemble_pred_raster, col = green_red(100), 
     main = "Ensemble Model - Distribuição potencial de Terminalia aucuminata")
points(ocorrencias, pch = 19, cex = 0.5)

# 13. COMPARAÇÃO DOS MODELOS E IDENTIFICAÇÃO DO MELHOR ----------------------

# Criar tabela comparativa
model_comparison <- data.frame(
  Modelo = c("MaxEnt", "GLM", "Random Forest", "GAM", "Ensemble"),
  AUC = c(round(as.numeric(maxent_auc$auc), 3),
          round(as.numeric(glm_auc$auc), 3),
          round(as.numeric(rf_auc$auc), 3),
          round(as.numeric(gam_auc$auc), 3),
          round(as.numeric(ensemble_auc$auc), 3)),
  TSS = c(maxent_metrics$TSS,
          glm_metrics$TSS,
          rf_metrics$TSS,
          gam_metrics$TSS,
          ensemble_metrics$TSS),
  Kappa = c(maxent_metrics$Kappa,
            glm_metrics$Kappa,
            rf_metrics$Kappa,
            gam_metrics$Kappa,
            ensemble_metrics$Kappa),
  Accuracy = c(maxent_metrics$Accuracy,
               glm_metrics$Accuracy,
               rf_metrics$Accuracy,
               gam_metrics$Accuracy,
               ensemble_metrics$Accuracy),
  Threshold = c(round(maxent_optimal$threshold, 3),
                round(glm_optimal$threshold, 3),
                round(rf_optimal$threshold, 3),
                round(gam_optimal$threshold, 3),
                round(ensemble_optimal$threshold, 3))
)

# Ordenar por AUC
model_comparison <- model_comparison[order(model_comparison$AUC, decreasing = TRUE), ]

print(model_comparison)

# Identificar melhor modelo
best_model <- model_comparison$Modelo[1]
cat("\n*** MELHOR MODELO:", best_model, "***\n")
cat("AUC:", model_comparison$AUC[1], "\n")
cat("TSS:", model_comparison$TSS[1], "\n")

# 14. PLOTAR CURVAS ROC ------------------------------------------------------

# Criar diretório
dir.create("resultados", showWarnings = FALSE)

# Plotar e salvar curvas ROC
png("resultados/curvas_roc_completas.png", width = 1000, height = 800, res = 100)
par(mar = c(5, 4, 4, 8), xpd = TRUE)

# Plotar curvas
plot(maxent_auc, col = "red", lwd = 2, 
     main = "Curvas ROC - Comparação de Todos os Modelos")
lines(glm_auc, col = "blue", lwd = 2)
lines(rf_auc, col = "green", lwd = 2)
lines(gam_auc, col = "orange", lwd = 2)
lines(ensemble_auc, col = "purple", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "gray")

# Legenda
legend("bottomright", 
       legend = c(paste("MaxEnt (AUC =", round(as.numeric(maxent_auc$auc), 3), ")"),
                  paste("GLM (AUC =", round(as.numeric(glm_auc$auc), 3), ")"),
                  paste("Random Forest (AUC =", round(as.numeric(rf_auc$auc), 3), ")"),
                  paste("GAM (AUC =", round(as.numeric(gam_auc$auc), 3), ")"),
                  paste("Ensemble (AUC =", round(as.numeric(ensemble_auc$auc), 3), ")")),
       col = c("red", "blue", "green", "orange", "purple"),
       lwd = 2,
       cex = 0.9)



# 15. GERAR SHAPEFILE DO POLÍGONO DE DISTRIBUIÇÃO ---------------------------

# Selecionar o melhor modelo para gerar o shapefile
if(best_model == "MaxEnt") {
  best_pred_raster <- maxent_pred
  best_threshold <- maxent_optimal$threshold
} else if(best_model == "GLM") {
  best_pred_raster <- glm_pred_raster
  best_threshold <- glm_optimal$threshold
} else if(best_model == "Random Forest") {
  best_pred_raster <- rf_pred_raster
  best_threshold <- rf_optimal$threshold
} else if(best_model == "GAM") {
  best_pred_raster <- gam_pred_raster
  best_threshold <- gam_optimal$threshold
} else {
  best_pred_raster <- ensemble_pred_raster
  best_threshold <- ensemble_optimal$threshold
}

# Criar raster binário
binary_raster <- best_pred_raster >= best_threshold
binary_raster[binary_raster == 0] <- NA

# Converter para polígono
distribution_poly <- as.polygons(binary_raster, dissolve = TRUE)

# Simplificar geometria (opcional, reduz tamanho do arquivo)
if(require("terra")) {
  distribution_poly <- simplifyGeom(distribution_poly, tolerance = 0.01)
}

# Adicionar atributos
distribution_poly$Species <- "Terminalia aucuminata"
distribution_poly$Model <- best_model
distribution_poly$Threshold <- best_threshold
distribution_poly$AUC <- model_comparison$AUC[model_comparison$Modelo == best_model]
distribution_poly$TSS <- model_comparison$TSS[model_comparison$Modelo == best_model]

# Salvar shapefile
writeVector(distribution_poly, "resultados/distribuicao_T_aucuminata.shp", 
            overwrite = TRUE)

# Visualizar
plot(distribution_poly, main = paste("Distribuição de T. aucuminata -", best_model))
points(ocorrencias, pch = 19, cex = 0.5, col = "red")

# 16. SALVAR TODOS OS RESULTADOS ---------------------------------------------
# Criar diretório principal
dir.create("resultados/rasters", recursive = TRUE, showWarnings = FALSE)

# Salvar rasters de predição contínua
writeRaster(maxent_pred, "resultados/rasters/maxent_prediction.tif", overwrite = TRUE)
writeRaster(glm_pred_raster, "resultados/rasters/glm_prediction.tif", overwrite = TRUE)
writeRaster(rf_pred_raster, "resultados/rasters/rf_prediction.tif", overwrite = TRUE)
writeRaster(gam_pred_raster, "resultados/rasters/gam_prediction.tif", overwrite = TRUE)
writeRaster(ensemble_pred_raster, "resultados/rasters/ensemble_prediction.tif", overwrite = TRUE)

# Salvar rasters binários (presença/ausência)
maxent_binary <- maxent_pred >= maxent_optimal$threshold
glm_binary <- glm_pred_raster >= glm_optimal$threshold
rf_binary <- rf_pred_raster >= rf_optimal$threshold
gam_binary <- gam_pred_raster >= gam_optimal$threshold
ensemble_binary <- ensemble_pred_raster >= ensemble_optimal$threshold

writeRaster(maxent_binary, "resultados/rasters/maxent_binary.tif", overwrite = TRUE)
writeRaster(glm_binary, "resultados/rasters/glm_binary.tif", overwrite = TRUE)
writeRaster(rf_binary, "resultados/rasters/rf_binary.tif", overwrite = TRUE)
writeRaster(gam_binary, "resultados/rasters/gam_binary.tif", overwrite = TRUE)
writeRaster(ensemble_binary, "resultados/rasters/ensemble_binary.tif", overwrite = TRUE)

# Salvar modelos
saveRDS(maxent_model, "resultados/maxent_model.rds")
saveRDS(glm_model, "resultados/glm_model.rds")
saveRDS(rf_model, "resultados/rf_model.rds")
saveRDS(gam_model, "resultados/gam_model.rds")

# Salvar tabelas de métricas
write.csv(model_comparison, "resultados/comparacao_modelos.csv", row.names = FALSE)

# Salvar thresholds
thresholds_df <- data.frame(
  Modelo = c("MaxEnt", "GLM", "Random Forest", "GAM", "Ensemble"),
  Threshold = c(maxent_optimal$threshold, glm_optimal$threshold, 
                rf_optimal$threshold, gam_optimal$threshold,
                ensemble_optimal$threshold)
)
write.csv(thresholds_df, "resultados/thresholds_otimos.csv", row.names = FALSE)

# Limpar memória
gc()


