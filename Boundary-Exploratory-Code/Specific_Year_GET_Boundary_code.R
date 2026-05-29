library(dplyr)
library(hetGP)
library(splines)
library(ggplot2)
library(MASS)
library(GET)
library(rnaturalearth)
library(sf)

# --- 1. DATA PREPARATION FUNCTION ---
# Run this once to load all data and attach basis functions
# --- 1. DATA PREPARATION FUNCTION (UPDATED) ---
# Run this once to load all data and attach basis functions
prep_data <- function(data_dir = ".", start_year = 1960, end_year = 1989) {
  a <- 1
  time <- start_year:end_year
  time_scaled <- scale(time, center = TRUE, scale = FALSE)
  T_periods <- c(18, 3, 6, 9, 12, 15)
  
  # Generate harmonic basis functions
  rbf_basis <- matrix(NA, nrow = length(time), ncol = 24)
  idx <- 1
  for(tp in T_periods) {
    for(har in 1:2) {
      rbf_basis[, idx]   <- a * sin(2 * pi * har * time_scaled / tp)
      rbf_basis[, idx+1] <- a * cos(2 * pi * har * time_scaled / tp)
      idx <- idx + 2
    }
  }
  col_names <- paste0("t", 1:24)
  colnames(rbf_basis) <- col_names
  basis_df <- data.frame(year = time, rbf_basis)
  
  # Load data
  get_data <- function(year) {
    file_name <- file.path(data_dir, paste0('Two_Lines', year, '_dup.csv'))
    if(!file.exists(file_name)) return(NULL)
    df <- read.csv(file_name)
    df$year <- as.numeric(year)
    return(df)
  }
  
  all_data <- lapply(start_year:end_year, get_data) %>% 
    bind_rows() %>%
    inner_join(basis_df, by = "year") %>%
    # --- NEW CONDITIONAL SPATIAL FILTER ---
    filter(
      (CLASS == "Des_Line" & ((lon < 30) | (lon >= 30 & lat > 9))) |
        (CLASS == "Non_Line" & (lon < 30))
    ) %>%
    arrange(year, lon)
  
  return(all_data)
}

# --- 2. EMPIRICAL GET TEST FUNCTION ---
run_get_test_empirical <- function(model, x_scaled, actual_lat, null_mean, n_sims = 2500, alpha = 0.05) {
  x_mat <- as.matrix(x_scaled)
  p_basic <- predict(object = model, x = x_mat)
  
  # Extract Latent Covariance
  Kss <- cov_gen(X1 = x_mat, X2 = x_mat, theta = model$theta, type = model$covtype)
  Kxs <- cov_gen(X1 = model$X0, X2 = x_mat, theta = model$theta, type = model$covtype)
  Sigma_latent <- model$nu_hat * (Kss - t(Kxs) %*% model$Ki %*% Kxs)
  
  # Total Covariance + Symmetry Fix
  Sigma_total <- Sigma_latent + diag(as.vector(p_basic$nugs))
  Sigma_total <- (Sigma_total + t(Sigma_total)) / 2 
  
  # Simulate under the null
  obs_diff <- as.vector(actual_lat) - as.vector(null_mean)
  Ynull1 <- MASS::mvrnorm(n = n_sims, mu = null_mean, Sigma = Sigma_total)
  Ynull2 <- MASS::mvrnorm(n = n_sims, mu = null_mean, Sigma = Sigma_total)
  
  # Compute empirical baseline and variance
  T0 <- colMeans(Ynull1 - Ynull2)
  Tvar <- apply(Ynull1 - Ynull2, 2, var)
  
  # Standardized envelope statistic
  R_i <- sort(apply(Ynull1 - Ynull2, 1, function(row) max(abs((row - T0) / sqrt(Tvar)))), decreasing = TRUE)
  R_t <- max(abs((obs_diff - T0) / sqrt(Tvar)))
  
  q_index <- ceiling(alpha * length(R_i))
  critical_value <- R_i[q_index]
  
  upper_env <- T0 + (critical_value * sqrt(Tvar))
  lower_env <- T0 - (critical_value * sqrt(Tvar))
  p_val <- sum(R_i >= R_t) / length(R_i)
  
  return(list(
    resid = obs_diff, mean = T0, ymin = lower_env, ymax = upper_env, p_val = p_val
  ))
}


rolling_window <- function(full_data, line_class = "Des_Line", test_year = 1983, window_size = 7, gap_years = 1) {
  
  target_data <- full_data %>% filter(CLASS == line_class)
  
  # Fit Harmonic Mean Model
  col_names <- paste0("t", 1:24)
  formula_lm <- as.formula(paste("lat ~", paste(col_names, collapse = " + ")))
  model_lm <- lm(formula_lm, data = target_data)
  
  # Prep residuals
  target_data$resids <- model_lm$residuals
  lon_min <- min(target_data$lon); lon_max <- max(target_data$lon)
  target_data$Scaled_X <- (target_data$lon - lon_min) / (lon_max - lon_min)
  
  # Subset Window
  train_end <- test_year - gap_years - 1
  train_start <- train_end - window_size + 1
  train_data <- target_data %>% filter(year >= train_start, year <= train_end)
  
  if(nrow(train_data) == 0) stop("No training data found.")
  
  # Fit GP
  model_gp <- mleHetGP(X = as.matrix(train_data$Scaled_X), 
                       Z = as.matrix(train_data$resids), 
                       covtype = "Matern3_2")
  
  # Predict Test Year (Silencing the 'Rank-Deficient' warning)
  test_data <- target_data %>% filter(year == test_year)
  lm_pred <- suppressWarnings(predict(model_lm, newdata = test_data))
  gp_pred <- predict(x = as.matrix(test_data$Scaled_X), object = model_gp)
  total_mean <- lm_pred + gp_pred$mean
  
  # Run GET
  res_emp <- run_get_test_empirical(model_gp, test_data$Scaled_X, test_data$lat, total_mean)
  
  # Prepare Plotting
  df_plot <- data.frame(
    lon = test_data$lon,
    lat_actual = test_data$lat,
    resid = res_emp$resid,
    env_mean = res_emp$mean,
    env_up = res_emp$ymax,
    env_low = res_emp$ymin
  )
  
  plot_resid <- ggplot(df_plot, aes(x = lon)) +
    geom_ribbon(aes(ymin = env_low, ymax = env_up), fill = "skyblue", alpha = 0.25) +
    geom_line(aes(y = env_mean), color = "black") +
    geom_line(aes(y = env_low), color = "blue", linetype = "dashed") +
    geom_line(aes(y = env_up), color = "blue", linetype = "dashed") +
    geom_line(aes(y = resid), color = "red", linewidth = 1) +
    labs(title = paste(line_class, test_year, "| p:", round(res_emp$p_val, 3)),
         x = "Longitude", y = "Residual Lat (Obs - Pred)") +
    theme_minimal()
  
  return(list(p_val = res_emp$p_val, plot_resid = plot_resid))
}



# 1. Load the data once (Set your specific directory path here)
setwd("C:/Users/steph/Downloads/Ninety_Sixty_Ninety_Eighty")
master_data <- prep_data(data_dir = ".") 

# 2. Test the Desert Line for 1983 (7 year window, 1 year gap)
results_des_83 <- rolling_window(
  full_data   = master_data, 
  line_class  = "Des_Line", 
  test_year   = 1983, 
  window_size = 7, 
  gap_years   = 1
)

# Print plots
print(results_des_83$plot_resid)

# 2. Test the Desert Line for 1983 (7 year window, 1 year gap)
results_non_83 <- rolling_window(
  full_data   = master_data, 
  line_class  = "Non_Line", 
  test_year   = 1983, 
  window_size = 7, 
  gap_years   = 1
)

# Print plots
print(results_non_83$plot_resid)



# 3. Test the Non-Desert Line for 1985 (10 year window, 2 year gap)
results_non_84 <- rolling_window(
  full_data   = master_data, 
  line_class  = "Non_Line", # Adjust this to match your exact CSV spelling 
  test_year   = 1984, 
  window_size = 7, 
  gap_years   = 1
)

print(results_non_84$plot_resid)


results_des_84 <- rolling_window(
  full_data   = master_data, 
  line_class  = "Des_Line", # Adjust this to match your exact CSV spelling 
  test_year   = 1984, 
  window_size = 7, 
  gap_years   = 1
)

print(results_des_84$plot_resid)


plot_path <- 'C:/Users/steph/Downloads/Picture_Paper_5_22_2026'
if(!dir.exists(plot_path)) dir.create(plot_path)

ggsave(file.path(plot_path, "primary_83_GET.png"), plot = results_des_83$plot_resid, width = 6, height = 4, dpi = 300, bg = "white")
ggsave(file.path(plot_path,"secondary_83_GET.png"), plot = results_non_83$plot_resid, width = 6, height = 4, dpi = 300, bg = "white")
ggsave(file.path(plot_path, "primary_84_GET.png"), plot = results_des_84$plot_resid, width = 6, height = 4, dpi = 300, bg = "white")
ggsave(file.path(plot_path, "secondary_84_GET.png"), plot = results_non_84$plot_resid, width = 6, height = 4, dpi = 300, bg = "white")

print("Analysis Complete. Plots saved.")
