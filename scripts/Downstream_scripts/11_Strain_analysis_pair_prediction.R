################################################################################
##### LL-NEXT: Viral strain transmission - Mother-infant pair prediction
### Author(s): Asier Fernández-Pato
### Last updated: 8th December, 2025
################################################################################

#****************
# Load libraries
#****************
library(dplyr)
library(tidymodels)
library(themis)       
library(pROC)  
library(purrr)

#****************
# Define function
#****************
# Function to run a logistic regression model to predict the phenotype
logreg_prediction <- function(data, outcome_col, label = NULL) {
  
  # Rename outcome column to Phenotype 
  data <- data %>%
    rename(Phenotype = {{outcome_col}})
  
  data$Phenotype <- factor(data$Phenotype)
  
  # Recipe
  recipe_data <- recipe(Phenotype ~ ., data = data) %>%
    step_normalize(all_predictors()) %>%
    step_downsample(Phenotype) # gives very similar results step_smote(Phenotype) 
  
  # Cross-validation (5-fold cross-validation, repeated 10 times)
  cv_folds <- vfold_cv(data, v = 5, repeats = 10, strata = Phenotype)
  
  # Logistic regression (no penalty tunning as only 1 predictor is used)
  logreg_spec <- logistic_reg() %>%
    set_engine("glm") %>%
    set_mode("classification")
  
  # Workflow (combine the model and preprocessing steps into a single object)
  logreg_workflow <- workflow() %>%
    add_model(logreg_spec) %>%
    add_recipe(recipe_data)
  
  # Cross-validated predictions (Run the model on each CV fold, save predictions from the test sets)
  logreg_cv_preds <- fit_resamples(
    logreg_workflow,
    resamples = cv_folds,
    metrics = metric_set(roc_auc),
    control = control_resamples(save_pred = TRUE)
  )
  
  # Combine predictions & compute ROC
  all_preds <- collect_predictions(logreg_cv_preds)
  names(all_preds)
  roc_obj <- pROC::roc(
    response  = all_preds$Phenotype,
    predictor = all_preds$.pred_Pair,        
    levels    = c("Not pair", "Pair"),
    direction = "<"
  )
  
  list(roc = roc_obj, label = label, auc = auc(roc_obj))
}

# Function to extract ROC curve data frame
roc_to_df <- function(res) {
  data.frame(
    fpr = 1 - res$roc$specificities,
    tpr = res$roc$sensitivities,
    label = res$label,
    auc = round(res$auc, 3)
  )
}

# Set working directory 
setwd("~/Desktop/PhD/Projects/Virome LL-Next/Final_analysis/VIRUSES/")


##************************************************************************
# 1. Load inStrain results for the LL-NEXT samples 
#*************************************************************************

# Read processed inStrain results 
inStrain_results <- read.delim("10_STRAIN_TRANSMISSION/inStrain_results_processed.txt")

##************************************************************************
# 2. Prediction of maternal infant pairs
#*************************************************************************
inStrain_results_W2 <- inStrain_results[inStrain_results$Infant_timepoint == "W2", c("popANI", "Mother_Infant_pair")]
inStrain_results_M1 <- inStrain_results[inStrain_results$Infant_timepoint == "M1", c("popANI", "Mother_Infant_pair")]
inStrain_results_M2 <- inStrain_results[inStrain_results$Infant_timepoint == "M2", c("popANI", "Mother_Infant_pair")]
inStrain_results_M3 <- inStrain_results[inStrain_results$Infant_timepoint == "M3", c("popANI", "Mother_Infant_pair")]
inStrain_results_M6 <- inStrain_results[inStrain_results$Infant_timepoint == "M6", c("popANI", "Mother_Infant_pair")]
inStrain_results_M9 <- inStrain_results[inStrain_results$Infant_timepoint == "M9", c("popANI", "Mother_Infant_pair")]
inStrain_results_M12 <- inStrain_results[inStrain_results$Infant_timepoint == "M12", c("popANI", "Mother_Infant_pair")]

# Run prediction for each timepoint
res_W2  <- logreg_prediction(inStrain_results_M1,  outcome_col = Mother_Infant_pair, label = "W2")
res_M1  <- logreg_prediction(inStrain_results_M1,  outcome_col = Mother_Infant_pair, label = "M1")
res_M2  <- logreg_prediction(inStrain_results_M2,  outcome_col = Mother_Infant_pair, label = "M2")
res_M3  <- logreg_prediction(inStrain_results_M3,  outcome_col = Mother_Infant_pair, label = "M3")
res_M6  <- logreg_prediction(inStrain_results_M6,  outcome_col = Mother_Infant_pair, label = "M6")
res_M9  <- logreg_prediction(inStrain_results_M9,  outcome_col = Mother_Infant_pair, label = "M9")
res_M12 <- logreg_prediction(inStrain_results_M12, outcome_col = Mother_Infant_pair, label = "M12")


# Extract results for ROC curve
# Bind all ROC data
roc_df <- bind_rows(roc_to_df(res_W2), roc_to_df(res_M1), roc_to_df(res_M2), 
                    roc_to_df(res_M3), roc_to_df(res_M6), roc_to_df(res_M9), roc_to_df(res_M12))

# Set colors for the curves
timepoint_order <- c("W2", "M1", "M2", "M3", "M6", "M9", "M12")
timepoint_colors <- c("#4B4FC5", "#3B81B3", "#4CA6B1", "#88CFA4", "#C1D97F", "#E8C26D", "#F4A219")

# Add label with AUC
roc_df <- roc_df %>%
  group_by(label, auc) %>%
  mutate(label_auc = paste0(label, " (AUC=", round(auc, 2), ")")) %>%
  ungroup()

# Reorder factor
roc_df$label_auc <- factor(
  roc_df$label_auc,
  levels = paste0(timepoint_order, " (AUC=", round(roc_df$auc[match(timepoint_order, roc_df$label)], 2), ")")
)

# Generate AUC plot
pdf('10_STRAIN_TRANSMISSION/Plots/ROC_pairs.pdf', width = 4, height = 3.5)
ggplot(roc_df, aes(x = fpr, y = tpr, color = label_auc)) +
  geom_line(size = 1, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", size = 0.9) +
  scale_color_manual(values = timepoint_colors) +
  labs(x = "False Positive Rate",
       y = "True Positive Rate",
       color = "Timepoint") +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    legend.position  = c(0.95, 0.2),    
    legend.justification = c("right","bottom"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75)
  )
dev.off()
