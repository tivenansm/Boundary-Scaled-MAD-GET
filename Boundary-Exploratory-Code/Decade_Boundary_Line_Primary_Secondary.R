# ==============================================================================
# 1. LIBRARIES AND SETTINGS
# ==============================================================================
library(dplyr)
library(hetGP)
library(ggplot2)
library(MASS)
library(sf)
library(rnaturalearth)

# Set your working directory
setwd("C:/Users/steph/Downloads/Ninety_Sixty_Ninety_Eighty")

# ==============================================================================
# 2. FUNCTIONS
# ==============================================================================

# Data Loader
# --- Updated Data Loading with Debugging ---
get_data <- function(year) {
  file <- paste0('Two_Lines', year, '_dup.csv')
  
  if(!file.exists(file)) {
    message(paste("Warning: File not found ->", file))
    return(NULL)
  }
  
  df <- read.csv(file)
  df$year <- as.numeric(year)
  return(df)
}

# Load data and check if it's empty
data_list <- lapply(time, get_data)
all_data <- bind_rows(data_list)

if (nrow(all_data) == 0) {
  stop("CRITICAL ERROR: No data was loaded. Check your working directory and file names.")
}

# Now perform the join
all_data <- all_data %>%
  inner_join(basis_df, by = "year")

# Comparison Plotting Function (GET)
compare_itcz_decades <- function(model_A, model_B, baseline_A, baseline_B, 
                                 test_grid, line_label = "ITCZ Line",
                                 n_sims = 2500, alpha = 0.05) {
  
  x_scaled <- as.matrix(test_grid[, 1])
  x_real   <- test_grid[, 2]
  
  # Predictions
  pA <- predict(object = model_A, x = x_scaled)
  pB <- predict(object = model_B, x = x_scaled)
  
  total_mean_A <- baseline_A + pA$mean
  total_mean_B <- baseline_B + pB$mean
  
  # Covariance Extraction Helper
  get_total_sigma <- function(mod, x_new, p_obj) {
    Kx  <- cov_gen(mod$X0, mod$X0, theta = mod$theta) + diag(mod$g, nrow(mod$X0))
    Kxs <- cov_gen(mod$X0, x_new, theta = mod$theta)
    Kss <- cov_gen(x_new, x_new, theta = mod$theta)
    Ki  <- solve(Kx)
    S <- mod$nu_hat * (Kss - t(Kxs) %*% Ki %*% Kxs) + diag(as.vector(p_obj$nugs))
    return((S + t(S)) / 2) # Force symmetry
  }
  
  Sigma_A <- get_total_sigma(model_A, x_scaled, pA)
  Sigma_B <- get_total_sigma(model_B, x_scaled, pB)
  
  # Simulation under the Null (Mean diff = 0)
  Ynull_A <- MASS::mvrnorm(n = n_sims, mu = rep(0, length(total_mean_A)), Sigma = Sigma_A)
  Ynull_B <- MASS::mvrnorm(n = n_sims, mu = rep(0, length(total_mean_B)), Sigma = Sigma_B)
  
  resid_sim <- Ynull_A - Ynull_B
  T0   <- colMeans(resid_sim)
  Tvar <- apply(resid_sim, 2, var)
  
  # Standardized GET Stat
  R_i <- sort(apply(resid_sim, 1, function(row) max(abs((row - T0) / sqrt(Tvar)))), decreasing = TRUE)
  obs_diff <- total_mean_A - total_mean_B
  R_t <- max(abs((obs_diff - T0) / sqrt(Tvar)))
  
  p_val <- sum(R_i >= R_t) / n_sims
  crit  <- R_i[ceiling(alpha * n_sims)]
  
  # Plotting DataFrame
  df_plot <- data.frame(
    x = x_real, obs = obs_diff,
    upper = T0 + (crit * sqrt(Tvar)),
    lower = T0 - (crit * sqrt(Tvar)),
    null = T0
  )
  
  # Plot
  p <- ggplot(df_plot, aes(x = x)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "skyblue", alpha = 0.3) +
    geom_line(aes(y = null), color = "black", linewidth = 0.8) +
    geom_line(aes(y = lower), color = "blue", linetype = "dashed", linewidth = 0.5) +
    geom_line(aes(y = upper), color = "blue", linetype = "dashed", linewidth = 0.5) +
    geom_line(aes(y = obs), color = "red", linewidth = .9) +
    labs(
         x = "Longitude", y = "Difference (Observed - Null Mean)") +
    theme_minimal(base_size = 14) +
    coord_cartesian(ylim = c(-4, 4))
  
  return(p)
}

# ==============================================================================
# 3. DATA PREPARATION & BASIS FUNCTIONS
# ==============================================================================

# Change this to "Des_Line" or "Non_Line"

time <- 1960:1989
time_scaled <- scale(time, center = TRUE, scale = FALSE)
T_periods <- c(6, 18, 3, 9, 12, 15)

basis_list <- lapply(T_periods, function(tp) {
  cbind(sin(2 * pi * 1 * time_scaled / tp), cos(2 * pi * 1 * time_scaled / tp),
        sin(2 * pi * 2 * time_scaled / tp), cos(2 * pi * 2 * time_scaled / tp))
})
rbf_basis <- do.call(cbind, basis_list)
colnames(rbf_basis) <- paste0("t", 1:24)
basis_df <- data.frame(year = time, rbf_basis)

all_data <- lapply(time, get_data) %>% 
  bind_rows() %>%
  inner_join(basis_df, by = "year")

# Filtering Logic
if(TARGET_CLASS == "Non_Line"){
  working_data <- all_data %>% filter(CLASS == "Non_Line", lon < 30)
} else {
  working_data <- all_data %>% filter(CLASS == "Des_Line", (lon < 30) | (lon >= 30 & lat > 9))
}

# ==============================================================================
# 4. LINEAR MODEL & GP TRAINING
# ==============================================================================

# Get Residuals
formula_lm <- as.formula(paste("lat ~", paste0("t", 1:24, collapse = " + ")))
model_lm <- lm(formula_lm, data = working_data)

# Create Training Object
X01 <- data.frame(
  lon = working_data$lon,
  Scaled_X = (working_data$lon - min(working_data$lon)) / (max(working_data$lon) - min(working_data$lon)),
  year = working_data$year,
  resids = as.vector(model_lm$residuals)
)

# Subset Decades
train60 <- subset(X01, year >= 1960 & year <= 1969)
train70 <- subset(X01, year >= 1970 & year <= 1979)
train80 <- subset(X01, year >= 1980 & year <= 1989)

# Fit GP Models
gp_60 <- mleHetGP(X = as.matrix(train60$Scaled_X), Z = as.matrix(train60$resids), covtype = "Matern3_2")
gp_70 <- mleHetGP(X = as.matrix(train70$Scaled_X), Z = as.matrix(train70$resids), covtype = "Matern3_2")
gp_80 <- mleHetGP(X = as.matrix(train80$Scaled_X), Z = as.matrix(train80$resids), covtype = "Matern3_2")

# ==============================================================================
# 5. HARD-CODED MEANS & PREDICTION GRID
# ==============================================================================

# Assign means based on class
if(TARGET_CLASS == "Non_Line") {
  Y_MEANS <- c(10.97166, 10.41671, 10.15725) # Non_Line means
} else {
  Y_MEANS <- c(13.35778, 12.77919, 12.53291) # Des_Line means (example)
}

# Create Grid
test_grid <- cbind(
  seq(0, 1, length = 500), 
  seq(min(X01$lon), max(X01$lon), length = 500)
)

# ==============================================================================
# 6. GENERATE & SAVE PLOTS
# ==============================================================================

# 1960s vs 1970s
p_60_70 <- compare_itcz_decades(
  model_A = gp_60, model_B = gp_70, 
  baseline_A = Y_MEANS[1], baseline_B = Y_MEANS[2], 
  test_grid = test_grid, line_label = paste(TARGET_CLASS, "60s vs 70s")
)

# 1970s vs 1980s
p_70_80 <- compare_itcz_decades(
  model_A = gp_70, model_B = gp_80, 
  baseline_A = Y_MEANS[2], baseline_B = Y_MEANS[3], 
  test_grid = test_grid, line_label = paste(TARGET_CLASS, "70s vs 80s")
)

# 1960s vs 1980s
p_60_80 <- compare_itcz_decades(
  model_A = gp_60, model_B = gp_80, 
  baseline_A = Y_MEANS[1], baseline_B = Y_MEANS[3], 
  test_grid = test_grid, line_label = paste(TARGET_CLASS, "60s vs 80s")
)



# Save
plot_path <- 'C:/Users/steph/Downloads/Picture_Paper_5_22_2026'
if(!dir.exists(plot_path)) dir.create(plot_path)

ggsave(file.path(plot_path, paste0(TARGET_CLASS, "_60_70_GET.png")), plot = p_60_70, width = 6, height = 4, dpi = 300, bg = "white")
ggsave(file.path(plot_path, paste0(TARGET_CLASS, "_70_80_GET.png")), plot = p_70_80, width = 6, height = 4, dpi = 300, bg = "white")
ggsave(file.path(plot_path, paste0(TARGET_CLASS, "_60_80_GET.png")), plot = p_60_80, width = 6, height = 4, dpi = 300, bg = "white")

print("Analysis Complete. Plots saved.")



####################################################################
####################################################################
####################################################################
####################################################################

# ==============================================================================
# 7. 1960s ABSOLUTE GLOBAL ENVELOPE TEST (GET) & POINT COUNTING
# ==============================================================================

# Extract baseline means
mean_60 <- Y_MEANS[1]
mean_70 <- Y_MEANS[2]
mean_80 <- Y_MEANS[3]

x_scaled <- as.matrix(test_grid[, 1])
x_real   <- test_grid[, 2]

# --- 1. Build the 1960s GET ---

# Predict 1960s Mean
p60 <- predict(object = gp_60, x = x_scaled)
total_mean_60 <- mean_60 + p60$mean

# Define Covariance Extraction Helper
get_total_sigma <- function(mod, x_new, p_obj) {
  Kx  <- cov_gen(mod$X0, mod$X0, theta = mod$theta) + diag(mod$g, nrow(mod$X0))
  Kxs <- cov_gen(mod$X0, x_new, theta = mod$theta)
  Kss <- cov_gen(x_new, x_new, theta = mod$theta)
  Ki  <- solve(Kx)
  S <- mod$nu_hat * (Kss - t(Kxs) %*% Ki %*% Kxs) + diag(as.vector(p_obj$nugs))
  return((S + t(S)) / 2) 
}

# Extract Covariance and Simulate
Sigma_60 <- get_total_sigma(gp_60, x_scaled, p60)
n_sims <- 2500
alpha <- 0.05

# Simulate curves to build the envelope
Ysim_60 <- MASS::mvrnorm(n = n_sims, mu = total_mean_60, Sigma = Sigma_60)

T0_60   <- colMeans(Ysim_60)
Tvar_60 <- apply(Ysim_60, 2, var)

# Standardized GET Stat for the single 1960s decade
R_i_60 <- sort(apply(Ysim_60, 1, function(row) max(abs((row - T0_60) / sqrt(Tvar_60)))), decreasing = TRUE)
crit_60 <- R_i_60[ceiling(alpha * n_sims)]

# Final 1960s Envelope Boundaries
upper_60 <- T0_60 + (crit_60 * sqrt(Tvar_60))
lower_60 <- T0_60 - (crit_60 * sqrt(Tvar_60))


# --- 2. COUNTING POINTS INSIDE THE 1960s GET ---

# METHOD A: Checking the RAW Observed Data Points
# Reconstruct the true latitude values for 70s/80s (baseline + residuals)
raw_lat_70 <- train70$resids + mean_70
raw_lat_80 <- train80$resids + mean_80

# Use linear interpolation to find the exact envelope bounds at the specific longitudes of the raw data
env_lower_70 <- approx(x = x_real, y = lower_60, xout = train70$lon)$y
env_upper_70 <- approx(x = x_real, y = upper_60, xout = train70$lon)$y

env_lower_80 <- approx(x = x_real, y = lower_60, xout = train80$lon)$y
env_upper_80 <- approx(x = x_real, y = upper_60, xout = train80$lon)$y

# Count how many raw points fall between the bounds
count_raw_70 <- sum(raw_lat_70 >= env_lower_70 & raw_lat_70 <= env_upper_70, na.rm = TRUE)
count_raw_80 <- sum(raw_lat_80 >= env_lower_80 & raw_lat_80 <= env_upper_80, na.rm = TRUE)

message("\n--- RAW DATA POINTS IN 1960s ENVELOPE ---")
message("1970s: ", count_raw_70, " out of ", nrow(train70), " raw points fall within the 60s GET.")
message("1980s: ", count_raw_80, " out of ", nrow(train80), " raw points fall within the 60s GET.")


# METHOD B: Checking the 500 predicted GP Grid Points
p70 <- predict(gp_70, x_scaled)
p80 <- predict(gp_80, x_scaled)

total_mean_70 <- mean_70 + p70$mean
total_mean_80 <- mean_80 + p80$mean

count_grid_70 <- sum(total_mean_70 >= lower_60 & total_mean_70 <= upper_60)
count_grid_80 <- sum(total_mean_80 >= lower_60 & total_mean_80 <= upper_60)

message("\n--- PREDICTED GP GRID POINTS IN 1960s ENVELOPE ---")
message("1970s: ", count_grid_70, " out of 500 GP predicted points fall within the 60s GET.")
message("1980s: ", count_grid_80, " out of 500 GP predicted points fall within the 60s GET.")