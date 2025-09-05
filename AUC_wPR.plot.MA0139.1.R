# --- add these packages ---
library(PRROC)      # PR curves + AUPRC
library(precrec)    # convenient multi-model PR/ROC (optional, used for plotting)
library(data.table)
library(pROC)
library(dplyr)
library(ggplot2)

# --- your data prep (unchanged) ---
df <- fread("MA0139.1.CTCF.COLO829.auc_motifctl_scored.bed") %>%
  mutate(motif = as.numeric(V4), meth = as.numeric(V6)) %>%
  mutate(label = if_else(V5 == "vacant", 0L, 1L)) %>%
  mutate(rule  = as.integer(motif == 1 & meth < 31.005))

# --- models (as you had) ---
df_model <- df[!is.na(meth)]
model <- glm(label ~ motif + meth, data = df_model, family = binomial)
df_model$pred <- predict(model, type = "response")  # higher = more positive

# --- Build a unified table of labels + scores for 4 models ---
# Important: for PR, scores must be higher for positives.
dt_scores <- data.table(
  label = df$label
)

# Motif is 0/1 already
dt_scores$Motif <- df$motif

# Meth: lower meth means more positive => flip sign so higher = more positive
# (You could also use max(meth) - meth; sign flip is enough for ranking-based metrics)
dt_scores$Meth <- -df$meth
dt_scores$Meth[is.na(dt_scores$Meth)] <- NA

# GLM (prob) available only where meth is not NA
dt_scores$`Combined (GLM)` <- NA_real_
dt_scores$`Combined (GLM)`[!is.na(df$meth)] <- df_model$pred

# Rule is 0/1 already
dt_scores$`Combined (Rule)` <- df$rule

# --- Helper to compute ROC AUC (pROC) and PR AUC (PRROC) for a single score vector ---
compute_metrics <- function(labels, scores, model_name) {
  ok <- !is.na(scores) & !is.na(labels)
  y  <- labels[ok]
  s  <- scores[ok]

  # ROC AUC with pROC (direction automatically detected)
  roc_obj <- roc(y, s, quiet = TRUE)
  roc_auc <- as.numeric(auc(roc_obj))

  # PR AUC (average precision) with PRROC
  # PRROC expects separate score vectors for positives and negatives
  pos_scores <- s[y == 1]
  neg_scores <- s[y == 0]
  pr_obj <- pr.curve(scores.class0 = pos_scores, scores.class1 = neg_scores, curve = TRUE)
  pr_auc <- as.numeric(pr_obj$auc.integral)

  list(
    name = model_name,
    roc_obj = roc_obj,
    pr_obj  = pr_obj,
    roc_auc = roc_auc,
    pr_auc  = pr_auc,
    prevalence = mean(y == 1)
  )
}

models <- c("Motif", "Meth", "Combined (GLM)", "Combined (Rule)")
res <- lapply(models, function(m) compute_metrics(dt_scores$label, dt_scores[[m]], m))

# --- Tabulate metrics ---
metrics <- rbindlist(lapply(res, function(x) {
  data.table(
    Model = x$name,
    ROC_AUC = round(x$roc_auc, 3),
    PR_AUC  = round(x$pr_auc, 3),
    Positive_Prevalence = round(x$prevalence, 3)
  )
}))
print(metrics)

# --- ROC plot (your style) ---
roc_df <- rbindlist(lapply(res, function(x) {
  data.table(
    FPR = 1 - x$roc_obj$specificities,
    TPR = x$roc_obj$sensitivities,
    Model = x$name
  )
}))
p_roc <- ggplot(roc_df, aes(FPR, TPR, color = Model)) +
  geom_line(size = 1.2) +
  geom_abline(linetype = "dashed") +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  labs(title = "ROC Curves", x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)") +
  theme(legend.position = "bottom")
print(p_roc)

# --- PR plot ---
# Use PRROC curves we already computed
pr_df <- rbindlist(lapply(res, function(x) {
  # PRROC returns a matrix with recall (x), precision (y)
  data.table(Recall = x$pr_obj$curve[,1],
             Precision = x$pr_obj$curve[,2],
             Model = x$name,
             Baseline = x$prevalence)
}))
# Add a baseline (horizontal) at prevalence (will vary per model if NA filtering differs);
# for simplicity, use the prevalence from GLM (largest filtered set) or compute global:
base_prev <- round(mean(dt_scores$label == 1, na.rm = TRUE), 3)

p_pr <- ggplot(pr_df, aes(Recall, Precision, color = Model)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = base_prev, linetype = "dotted") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  labs(title = "Precision–Recall Curves",
       x = "Recall",
       y = paste0("Precision (baseline=", base_prev, ")")) +
  theme(legend.position = "bottom")
print(p_pr)

# --- Requirements ---
library(PRROC)
library(data.table)

# labels: integer/numeric 0/1 vector (length N)
# scores_dt: data.frame/data.table with one column per model (same length/order as labels)
#            e.g., scores_dt[, c("Motif","Meth","Combined (GLM)","Combined (Rule)")]
boot_pr_auc_multi <- function(labels, scores_dt, n_boot = 2000, conf = 0.95, stratified = TRUE, seed = 123) {
  stopifnot(length(labels) == nrow(scores_dt))
  set.seed(seed)
  
  # helper for one model
  one_model <- function(y, s) {
    ok <- !is.na(y) & !is.na(s)
    y  <- y[ok]; s <- s[ok]
    pos <- s[y == 1]; neg <- s[y == 0]
    # point estimate
    pe <- pr.curve(scores.class0 = pos, scores.class1 = neg)$auc.integral
    
    # bootstrap
    n_pos <- length(pos); n_neg <- length(neg)
    aucs <- numeric(n_boot)
    for (b in seq_len(n_boot)) {
      if (stratified) {
        pos_b <- sample(pos, n_pos, replace = TRUE)
        neg_b <- sample(neg, n_neg, replace = TRUE)
      } else {
        # non-stratified: resample indices and then split
        idx <- sample(seq_along(s), length(s), replace = TRUE)
        yb  <- y[idx]; sb <- s[idx]
        pos_b <- sb[yb == 1]; neg_b <- sb[yb == 0]
      }
      aucs[b] <- pr.curve(scores.class0 = pos_b, scores.class1 = neg_b)$auc.integral
    }
    
    ci <- as.numeric(quantile(aucs, probs = c((1 - conf)/2, 1 - (1 - conf)/2)))
    list(
      pr_auc = unname(pe),
      ci_lo  = ci[1],
      ci_hi  = ci[2],
      n_pos  = n_pos,
      n_neg  = n_neg,
      prev   = n_pos / (n_pos + n_neg)
    )
  }
  
  # run for all score columns
  out <- rbindlist(lapply(names(scores_dt), function(m) {
    res <- one_model(labels, scores_dt[[m]])
    data.table(
      Model = m,
      PR_AUC = round(res$pr_auc, 3),
      CI_L   = round(res$ci_lo, 3),
      CI_U   = round(res$ci_hi, 3),
      Positives = res$n_pos,
      Negatives = res$n_neg,
      Prevalence = round(res$prev, 3)
    )
  }))
  
  setorder(out, -PR_AUC)
  out
}

score_cols <- c("Motif", "Meth", "Combined (GLM)", "Combined (Rule)")
pr_ci_table <- boot_pr_auc_multi(labels = dt_scores$label,
                                 scores_dt = dt_scores[, ..score_cols],
                                 n_boot = 2000, conf = 0.95, stratified = TRUE, seed = 1)
print(pr_ci_table)




