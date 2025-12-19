# =========================================================
# Task 1: Preliminary Data Analysis (EDA) - Robust R Script (Mac)
# File: dataset_(MSC7)-80971f18-5843-4b66-8d97-d159a8118873 (1).csv
# Location: ~/Downloads
# =========================================================

library(tidyverse)
library(corrplot)
library(readr)
library(dplyr)
library(ggplot2)
# 1) Set working directory to Downloads (Mac)
setwd(file.path(Sys.getenv("HOME"), "Downloads"))

# 2) Read CSV (check.names=TRUE makes names R-safe automatically)
df <- read.csv("dataset_(MSC7)-80971f18-5843-4b66-8d97-d159a8118873 (1).csv",
               check.names = TRUE)

# 3) Show column names so you can verify
cat("\n--- Column names in dataset (R-safe) ---\n")
print(colnames(df))

# 4) Identify the target column (bg+1:00 becomes something like bg.1.00)
#    We search for a "bg" column that contains "1" and "00"
target_candidates <- names(df)[grepl("^bg", names(df), ignore.case = TRUE) &
                                 grepl("1", names(df)) &
                                 grepl("00", names(df))]

if (length(target_candidates) == 0) {
  stop("Target column not found. Please check printed column names above.")
}

# If multiple matches, take the first and warn
target_col <- target_candidates[1]
if (length(target_candidates) > 1) {
  warning(paste("Multiple target candidates found:",
                paste(target_candidates, collapse = ", "),
                "\nUsing:", target_col))
}

# 5) Rename target to a clean name for modeling
df <- df %>% rename(bg_1h = all_of(target_col))

# 6) Define predictors (must exist exactly)
predictors <- c("bg_mean", "insulin_sum", "carbs_sum", "hr_mean", "steps_sum", "cals_sum")

missing_preds <- setdiff(predictors, names(df))
if (length(missing_preds) > 0) {
  stop(paste("Missing predictor columns:", paste(missing_preds, collapse = ", "),
             "\nCheck printed column names above and adjust predictors list."))
}

# 7) Keep only required columns and drop NA
df2 <- df %>%
  dplyr::select(bg_1h, all_of(predictors)) %>%
  drop_na()

# 8) Create time index
df2$time_index <- 1:nrow(df2)

# =========================================================
# A) Time-series plots
# =========================================================

# Target time-series
p_ts_target <- ggplot(df2, aes(x = time_index, y = bg_1h)) +
  geom_line() +
  labs(title = "Time Series: Future Blood Glucose (bg+1:00)",
       x = "Time Index",
       y = "bg+1:00 (mmol/L)") +
  theme_minimal()
print(p_ts_target)

# Predictors time-series (facet)
df_long <- df2 %>%
  pivot_longer(cols = all_of(predictors),
               names_to = "variable",
               values_to = "value")

p_ts_pred <- ggplot(df_long, aes(x = time_index, y = value)) +
  geom_line() +
  facet_wrap(~variable, scales = "free_y", ncol = 2) +
  labs(title = "Time Series: Predictor Variables",
       x = "Time Index",
       y = "Value") +
  theme_minimal()
print(p_ts_pred)

# =========================================================
# B) Distribution plots
# =========================================================

# Target distribution
p_dist_target <- ggplot(df2, aes(x = bg_1h)) +
  geom_histogram(bins = 30) +
  labs(title = "Distribution: bg+1:00",
       x = "bg+1:00 (mmol/L)",
       y = "Frequency") +
  theme_minimal()
print(p_dist_target)

# Predictor distributions
p_dist_pred <- ggplot(df_long, aes(x = value)) +
  geom_histogram(bins = 30) +
  facet_wrap(~variable, scales = "free", ncol = 2) +
  labs(title = "Distributions: Predictor Variables",
       x = "Value",
       y = "Frequency") +
  theme_minimal()
print(p_dist_pred)

# =========================================================
# C) Correlation matrix + heatmap
# =========================================================

cor_mat <- cor(df2 %>% dplyr::select(bg_1h, all_of(predictors)))
cat("\n--- Correlation Matrix ---\n")
print(round(cor_mat, 3))

corrplot(cor_mat,
         method = "color",
         type = "upper",
         addCoef.col = "black",
         tl.col = "black",
         number.cex = 0.8,
         title = "Correlation Matrix (Target + Predictors)",
         mar = c(0, 0, 2, 0))

# =========================================================
# D) Scatter plots (predictors vs target)
# =========================================================

for (p in predictors) {
  g <- ggplot(df2, aes(x = .data[[p]], y = bg_1h)) +
    geom_point(alpha = 0.3) +
    labs(title = paste("Scatter Plot:", p, "vs bg+1:00"),
         x = p,
         y = "bg+1:00 (mmol/L)") +
    theme_minimal()
  print(g)
}

# =========================================================
# Optional: Save plots
# =========================================================
# ggsave("01_timeseries_target.png", p_ts_target, width = 8, height = 4)
# ggsave("02_timeseries_predictors.png", p_ts_pred, width = 10, height = 7)
# ggsave("03_dist_target.png", p_dist_target, width = 8, height = 4)
# ggsave("04_dist_predictors.png", p_dist_pred, width = 10, height = 7)
# ---- Load data ----
file_path <- "dataset_(MSC7)-80971f18-5843-4b66-8d97-d159a8118873 (1).csv"  # <-- change path if needed
df <- read_csv(file_path, show_col_types = FALSE)

# Keep only required variables + drop missing rows (needed for least squares)
data <- df %>%
  select(bg_mean, insulin_sum, carbs_sum, hr_mean, steps_sum, cals_sum, `bg+1:00`) %>%
  na.omit()

# Define y and x’s
y  <- data$`bg+1:00`
x1 <- data$bg_mean
x2 <- data$insulin_sum
x3 <- data$carbs_sum
x4 <- data$hr_mean
x5 <- data$steps_sum
x6 <- data$cals_sum

# ============================================================
# Task 2.1: Estimate parameters for every candidate model using Least Squares
# ============================================================

# Model 1: y = b1*x1^3 + b2*x2^2 + b3*x3^2 + b4*x4 + b5*x5 + b6*x6 + b0
fit1 <- lm(y ~ I(x1^3) + I(x2^2) + I(x3^2) + x4 + x5 + x6)

# Model 2: y = b1*x1^2 + b2*x2^2 + b3*x3^3 + b4*x4 + b5*x5 + b6*x6 + b0
fit2 <- lm(y ~ I(x1^2) + I(x2^2) + I(x3^3) + x4 + x5 + x6)

# Model 3: y = b1*x1 + b2*x2 + b3*x3 + b4*x4^2 + b5*x5 + b6*x6^2 + b0
fit3 <- lm(y ~ x1 + x2 + x3 + I(x4^2) + x5 + I(x6^2))

# Model 4: y = b1*x1^2 + b2*x2^2 + b3*x3^2 + b4*x4^2 + b5*x5^2 + b6*x6^2 + b0
fit4 <- lm(y ~ I(x1^2) + I(x2^2) + I(x3^2) + I(x4^2) + I(x5^2) + I(x6^2))

# Model 5: linear + interactions
# y = b1*x1 + b2*x2 + b3*x3 + b4*x4 + b5*x5 + b6*x6
#   + b7*x1*x2 + b8*x3*x4 + b9*x2*x6 + b0
fit5 <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + I(x1*x2) + I(x3*x4) + I(x2*x6))

models <- list(M1 = fit1, M2 = fit2, M3 = fit3, M4 = fit4, M5 = fit5)

# ============================================================
# Task 2.2: Compute RSS (sum of squared residuals)
# ============================================================
RSS <- function(fit) sum(residuals(fit)^2)

# ============================================================
# Task 2.3: Compute log-likelihood for every candidate model
# (Gaussian residual assumption; use logLik() from lm)
# ============================================================
# (Gaussian residual assumption)
# ============================================================

log_likelihoods <- sapply(models, function(fit) {
  as.numeric(logLik(fit))
})

# Convert to data frame for clarity
logLik_df <- data.frame(
  Model  = names(log_likelihoods),
  LogLik = log_likelihoods
)

print(logLik_df)
# ============================================================
# Task 2.4: Compute AIC and BIC for every candidate model
# ============================================================

results <- lapply(names(models), function(name){
  fit <- models[[name]]
  data.frame(
    Model  = name,
    RSS    = RSS(fit),
    LogLik = as.numeric(logLik(fit)),
    AIC    = AIC(fit),
    BIC    = BIC(fit),
    k      = length(coef(fit))   # number of parameters including intercept
  )
}) %>% bind_rows()

print(results)

# ============================================================
# Task 2.5: Residual distribution checks (Histogram + Q-Q plot)
# ============================================================

# Create an output folder for plots
out_dir <- "task2_residual_plots"
dir.create(out_dir, showWarnings = FALSE)

for(nm in names(models)){
  fit <- models[[nm]]
  r <- residuals(fit)
  
  # Histogram
  png(file.path(out_dir, paste0(nm, "_residual_hist.png")), width = 900, height = 650)
  hist(r, breaks = 40, main = paste(nm, "Residual Histogram"), xlab = "Residual")
  dev.off()
  
  # Q-Q plot
  png(file.path(out_dir, paste0(nm, "_residual_qq.png")), width = 900, height = 650)
  qqnorm(r, main = paste(nm, "Residual Q-Q Plot"))
  qqline(r)
  dev.off()
}

# ============================================================
# Task 2.6: Select the best model based on AIC, BIC, residual diagnostics
# ============================================================

best_AIC <- results$Model[which.min(results$AIC)]
best_BIC <- results$Model[which.min(results$BIC)]

cat("\nBest by AIC:", best_AIC, "\n")
cat("Best by BIC:", best_BIC, "\n")

# Note: You justify final choice using AIC/BIC + which residual plots look most Gaussian.
# ============================================================
# Task 2.7: Train/Test split + fit best model + predict on test
#          + 95% prediction intervals + plot with error bars
# ============================================================

set.seed(123)  # reproducibility

# ----- 70/30 split -----
n <- nrow(data)
train_idx <- sample(seq_len(n), size = floor(0.7 * n), replace = FALSE)

train_data <- data[train_idx, ]
test_data  <- data[-train_idx, ]

# ----- Choose the "best" model (pick one rule) -----
# Option A: Use best by BIC (often preferred for parsimony)
best_model_name <- best_BIC

# Option B (alternative): use best by AIC
# best_model_name <- best_AIC

cat("\nTask 2.7 using best model:", best_model_name, "\n")

# ----- Fit the selected model on TRAINING data -----
# IMPORTANT: rebuild the model formula so it works with train_data/test_data columns directly.
# (Avoid x1..x6 vectors here; use column names.)

if(best_model_name == "M1"){
  best_fit <- lm(`bg+1:00` ~ I(bg_mean^3) + I(insulin_sum^2) + I(carbs_sum^2) + hr_mean + steps_sum + cals_sum,
                 data = train_data)
} else if(best_model_name == "M2"){
  best_fit <- lm(`bg+1:00` ~ I(bg_mean^2) + I(insulin_sum^2) + I(carbs_sum^3) + hr_mean + steps_sum + cals_sum,
                 data = train_data)
} else if(best_model_name == "M3"){
  best_fit <- lm(`bg+1:00` ~ bg_mean + insulin_sum + carbs_sum + I(hr_mean^2) + steps_sum + I(cals_sum^2),
                 data = train_data)
} else if(best_model_name == "M4"){
  best_fit <- lm(`bg+1:00` ~ I(bg_mean^2) + I(insulin_sum^2) + I(carbs_sum^2) + I(hr_mean^2) + I(steps_sum^2) + I(cals_sum^2),
                 data = train_data)
} else if(best_model_name == "M5"){
  best_fit <- lm(`bg+1:00` ~ bg_mean + insulin_sum + carbs_sum + hr_mean + steps_sum + cals_sum +
                   I(bg_mean*insulin_sum) + I(carbs_sum*hr_mean) + I(insulin_sum*cals_sum),
                 data = train_data)
} else {
  stop("Unknown best model name. Expected one of: M1..M5")
}

# ----- Predict on TEST data + 95% prediction intervals -----
pred_pi <- predict(best_fit, newdata = test_data, interval = "prediction", level = 0.95)

test_results <- test_data %>%
  mutate(
    y_true = `bg+1:00`,
    y_pred = pred_pi[, "fit"],
    lo_95  = pred_pi[, "lwr"],
    hi_95  = pred_pi[, "upr"]
  )

# ----- Optional: compute test error metrics (useful in write-up) -----
rmse <- sqrt(mean((test_results$y_true - test_results$y_pred)^2))
mae  <- mean(abs(test_results$y_true - test_results$y_pred))

cat("Test RMSE:", rmse, "\n")
cat("Test MAE :", mae,  "\n")

# ----- Plot: test samples + predictions + 95% PI error bars -----
# For plotting, order points by prediction (or by bg_mean). Here we order by bg_mean for readability.
plot_df <- test_results %>%
  arrange(bg_mean) %>%
  mutate(idx = row_number())

p <- ggplot(plot_df, aes(x = idx)) +
  geom_point(aes(y = y_true), alpha = 0.8) +
  geom_point(aes(y = y_pred), alpha = 0.8) +
  geom_errorbar(aes(ymin = lo_95, ymax = hi_95), width = 0.2, alpha = 0.5) +
  labs(
    title = paste0("Task 2.7: Test Data vs Prediction (", best_model_name, ") with 95% Prediction Intervals"),
    x = "Test sample index (sorted by bg_mean)",
    y = "bg+1:00 (mmol/L)"
  )

print(p)

# Save plot
ggsave(filename = "Task2_7_Test_Prediction_with_95PI.png", plot = p, width = 10, height = 6, dpi = 300)
# ============================================================
# Task 3: Approximate Bayesian Computation (ABC) - Rejection ABC
# Selected model: M5
# Only 2 parameters: the 2 largest |beta| from LS estimation (Task 2.1)
# Fix all other parameters at LS estimates
# Prior: Uniform around LS values
# ============================================================

library(dplyr)
library(ggplot2)

# -----------------------------
# Use M5 least-squares fit
# -----------------------------
sel_fit <- fit5   # M5 from Task 2.1

# Build model frame + design matrix that exactly match M5 structure
abc_data <- model.frame(sel_fit)                 # data used in sel_fit
y_obs    <- model.response(abc_data)             # observed y
X        <- model.matrix(sel_fit, data = abc_data)
beta_hat <- coef(sel_fit)
sigma_hat <- summary(sel_fit)$sigma

# -----------------------------
# -----------------------------
# Step 1: pick top-2 |beta| (exclude intercept + remove NA coefficients)
# -----------------------------
beta_hat <- coef(sel_fit)

beta_rank <- beta_hat[names(beta_hat) != "(Intercept)"]
beta_rank <- beta_rank[!is.na(beta_rank)]   # IMPORTANT: drop aliased/NA coefficients

if(length(beta_rank) < 2) stop("Not enough non-NA coefficients to select top 2 parameters.")

top2_names <- names(sort(abs(beta_rank), decreasing = TRUE))[1:2]

cat("\nM5 - Top 2 parameters by |LS estimate| (non-NA):\n")
print(top2_names)
cat("\nLS estimates for top 2:\n")
print(beta_hat[top2_names])

# -----------------------------
# Step 2: define Uniform priors around LS estimates (make sure numeric)
# -----------------------------
prior_frac <- 0.30

prior_ranges <- lapply(top2_names, function(pn){
  est <- as.numeric(beta_hat[pn])
  if(is.na(est)) stop(paste("Selected parameter has NA estimate:", pn))
  w <- prior_frac * max(abs(est), 1e-6)
  c(lower = est - w, upper = est + w)
})
names(prior_ranges) <- top2_names

cat("\nUniform prior ranges:\n")
print(prior_ranges)

# Sample from priors (safe numeric extraction)
p1 <- top2_names[1]; p2 <- top2_names[2]

theta1 <- runif(N,
                min = as.numeric(prior_ranges[[p1]]["lower"]),
                max = as.numeric(prior_ranges[[p1]]["upper"]))

theta2 <- runif(N,
                min = as.numeric(prior_ranges[[p2]]["lower"]),
                max = as.numeric(prior_ranges[[p2]]["upper"]))


# -----------------------------
# Step 3: Rejection ABC
# - sample from priors
# - simulate y under M5 with only these 2 betas varying
# - accept the closest samples by distance (MSE)
# -----------------------------
set.seed(123)
N <- 50000          # number of prior draws (increase if you want smoother posterior)
accept_rate <- 0.01 # accept best 1% (rejection threshold)

p1 <- top2_names[1]
p2 <- top2_names[2]

j1 <- which(colnames(X) == p1)
j2 <- which(colnames(X) == p2)
if(length(j1) != 1 || length(j2) != 1) stop("Could not match parameter names to design matrix columns.")

# Fixed betas (we only replace p1 and p2 each draw)
beta_fixed <- beta_hat
base_mu <- as.vector(X %*% beta_fixed)

# Sample from uniform priors
theta1 <- runif(N, min = prior_ranges[[p1]]["lower"], max = prior_ranges[[p1]]["upper"])
theta2 <- runif(N, min = prior_ranges[[p2]]["lower"], max = prior_ranges[[p2]]["upper"])

dist_vec <- numeric(N)

for(i in 1:N){
  mu_i <- base_mu +
    (theta1[i] - beta_hat[p1]) * X[, j1] +
    (theta2[i] - beta_hat[p2]) * X[, j2]
  
  # Simulate response under Gaussian noise
  y_sim <- mu_i + rnorm(length(y_obs), mean = 0, sd = sigma_hat)
  
  # Distance metric (mean squared error)
  dist_vec[i] <- mean((y_sim - y_obs)^2)
}

# Accept the best fraction (smallest distances)
n_acc <- max(1, floor(accept_rate * N))
acc_idx <- order(dist_vec)[1:n_acc]

post <- data.frame(
  b1 = theta1[acc_idx],
  b2 = theta2[acc_idx],
  dist = dist_vec[acc_idx]
)
names(post)[1:2] <- c(p1, p2)

cat("\nAccepted samples:", nrow(post), "out of", N, "\n")
cat("Accepted distance summary:\n")
print(summary(post$dist))
print(summ)

# -----------------------------
# Step 4: Plot joint and marginal posterior distributions
# -----------------------------
# Joint posterior
p_joint <- ggplot(post, aes_string(x = p1, y = p2)) +
  geom_point(alpha = 0.25) +
  geom_density_2d() +
  labs(
    title = "Task 3: ABC Rejection Posterior (Joint) - M5",
    subtitle = paste0("Accepted ", n_acc, " / ", N, " (", accept_rate*100, "%)"),
    x = p1, y = p2
  )
print(p_joint)
ggsave("Task3_M5_ABC_joint_posterior.png", plot = p_joint, width = 8, height = 6, dpi = 300)

# Marginal posterior for p1
p_m1 <- ggplot(post, aes_string(x = p1)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.6) +
  geom_density() +
  labs(title = paste("Task 3: ABC Marginal Posterior -", p1, "(M5)"), x = p1, y = "Density")
print(p_m1)
ggsave(paste0("Task3_M5_ABC_marginal_", gsub("[^A-Za-z0-9_]+","_", p1), ".png"),
       plot = p_m1, width = 8, height = 5, dpi = 300)

# Marginal posterior for p2
p_m2 <- ggplot(post, aes_string(x = p2)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.6) +
  geom_density() +
  labs(title = paste("Task 3: ABC Marginal Posterior -", p2, "(M5)"), x = p2, y = "Density")
print(p_m2)
ggsave(paste0("Task3_M5_ABC_marginal_", gsub("[^A-Za-z0-9_]+","_", p2), ".png"),
       plot = p_m2, width = 8, height = 5, dpi = 300)

# -----------------------------
# Step 5: Numerical summary (for your explanation)
# -----------------------------
summ <- post %>%
  summarise(
    mean_p1 = mean(.data[[p1]]),
    sd_p1   = sd(.data[[p1]]),
    mean_p2 = mean(.data[[p2]]),
    sd_p2   = sd(.data[[p2]]),
    corr    = cor(.data[[p1]], .data[[p2]])
  )

cat("\nPosterior summary (accepted ABC samples):\n")
print(summ)

cat("\nExplanation hints:\n")
cat("- If posterior means are close to LS estimates, ABC agrees with LS under Gaussian noise.\n")
cat("- Narrower posterior => parameter better identified by the data.\n")
cat("- Wider posterior => more uncertainty / weaker identifiability.\n")
cat("- Strong correlation in joint posterior => trade-off between the two parameters.\n")
coef(fit5)[is.na(coef(fit5))]
print(top2_names)
print(coef(fit5))