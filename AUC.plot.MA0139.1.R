
library(pROC)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
df <- fread("MA0139.1.CTCF.COLO829.auc_motifctl_scored.bed") %>%
      mutate(motif = as.numeric(V4), meth = as.numeric(V6)) %>%
      mutate(label = case_when(V5 == "vacant" ~ 0, .default = 1)) %>%
      mutate(rule = as.integer(motif == 1 & meth < 31.005))

#decide best methylation cutoff: 31.005 for MA0139.1
df_sub <- subset(df, motif == 1 & !is.na(meth))
roc_meth_m1 <- roc(df_sub$label, df_sub$meth, quiet = TRUE)
best_coords <- coords(roc_meth_m1, x = "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"))
print(best_coords)

#motif roc
roc_motif <- roc(df$label, df$motif, ci = TRUE, quiet = TRUE)
#meth roc
roc_meth <- roc(df[!is.na(meth)]$label, df[!is.na(meth)]$meth, ci = TRUE, quiet = TRUE)
#glm roc
 df_model <- df[!is.na(meth)]
 #Fit logistic regression using both predictors
 model <- glm(label ~ motif + meth, data = df_model, family = binomial)
 #Predict probability
 df_model[, pred := predict(model, type = "response")]
 #ROC for combined model
 roc_combined <- roc(df_model$label, df_model$pred, ci = TRUE, quiet = TRUE)

#dummy rule roc
 df_rule <- df[!is.na(rule)]
 roc_rule <- roc(df_rule$label, df_rule$rule, ci = TRUE, quiet = TRUE)

#function extract info
roc_df <- function(roc_obj, name) {
  df <- data.frame(
    FPR = 1 - roc_obj$specificities,          # 1 - specificity
    TPR = roc_obj$sensitivities,
    Model = name
  )
  return(df)
}

# Combine all curves
df_motif  <- roc_df(roc_motif,  "Motif")
df_meth   <- roc_df(roc_meth,   "Meth")
df_glm    <- roc_df(roc_combined, "Combined (GLM)")
df_rule   <- roc_df(roc_rule,   "Combined (Rule)")

df_all <- rbind(df_motif, df_meth, df_glm, df_rule)

ggplot(df_all, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 1.2) +
  geom_abline(linetype = "dashed", color = "gray") +
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  labs( title = "ROC: MA0139.1", x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)", color = "Model") +
  coord_fixed(ratio = 1) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )
  
roc.test(roc_combined, roc_motif)
roc.test(roc_meth, roc_motif)
