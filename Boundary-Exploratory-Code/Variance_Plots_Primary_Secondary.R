library(dplyr)
library(hetGP)
library(splines)
library(ggplot2)
library(MASS)
library(GET)
library(rnaturalearth)
library(sf)

# --- 1. SETUP & BASIS FUNCTIONS ---
a <- 1
time <- 1960:1989
time_scaled <- scale(time, center = TRUE, scale = FALSE)
T_periods <- c(18, 3, 6, 9, 12, 15)

# Generate harmonic basis functions (24 total)
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

# --- 2. DATA LOADING & CLEANING ---
# setwd("C:/Users/steph/Downloads/Ninety_Sixty_Ninety_Eighty") # Set your path here

get_data <- function(year) {
  file_name <- paste0('Two_Lines', year, '_dup.csv')
  if(!file.exists(file_name)) return(NULL)
  df <- read.csv(file_name)
  df$year <- as.numeric(year)
  return(df)
}

# all_data <- lapply(1960:1989, get_data) %>% 
#   bind_rows() %>%
#   inner_join(basis_df, by = "year") %>%
#   filter(CLASS == "Des_Line") %>%
#   filter((lon < 30) | (lon >= 30 & lat > 9)) %>%
#   arrange(year, lon)


all_data <- lapply(1960:1989, get_data) %>% 
  bind_rows() %>%
  inner_join(basis_df, by = "year") %>%
  filter(CLASS == "Non_Line") %>%
  filter((lon < 30)) %>%
  arrange(year, lon)

# --- 3. HARMONIC MEAN MODEL ---
formula_lm <- as.formula(paste("lat ~", paste(col_names, collapse = " + ")))
model_lm <- lm(formula_lm, data = all_data)

# Add residuals to data
all_data$resids <- model_lm$residuals
lon_min <- min(all_data$lon); lon_max <- max(all_data$lon)
all_data$Scaled_X <- (all_data$lon - lon_min) / (lon_max - lon_min)


# Example usage (assuming 'model_gp' and 'data_83' are loaded):
# var_plot_83 <- plot_gp_variance(model = model_gp, test_data = data_83)
# print(var_plot_83)




# --- 4. HETEROSKEDASTIC GP TRAINING ---
# Training on 1976-1982 to predict 1983
train_data <- all_data %>% filter(year >= 1960, year <= 1969)
model_gp1 <- mleHetGP(X = as.matrix(train_data$Scaled_X), 
                     Z = as.matrix(train_data$resids), 
                     covtype = "Matern3_2")

train_data <- all_data %>% filter(year >= 1970, year <= 1979)
model_gp2 <- mleHetGP(X = as.matrix(train_data$Scaled_X), 
                      Z = as.matrix(train_data$resids), 
                      covtype = "Matern3_2")

train_data <- all_data %>% filter(year >= 1980, year <= 1989)
model_gp3 <- mleHetGP(X = as.matrix(train_data$Scaled_X), 
                      Z = as.matrix(train_data$resids), 
                      covtype = "Matern3_2")




train_data <- all_data %>% filter(year >= 1975, year <= 1981)
model_gp4 <- mleHetGP(X = as.matrix(train_data$Scaled_X), 
                      Z = as.matrix(train_data$resids), 
                      covtype = "Matern3_2")


train_data <- all_data %>% filter(year >= 1976, year <= 1982)
model_gp5 <- mleHetGP(X = as.matrix(train_data$Scaled_X), 
                      Z = as.matrix(train_data$resids), 
                      covtype = "Matern3_2")

plot_gp_simulations <- function(model, test_data, n_sims = 1000) {
  # 1. Scale inputs
  lon_min <- min(test_data$lon)
  lon_max <- max(test_data$lon)
  x_scaled <- (test_data$lon - lon_min) / (lon_max - lon_min)
  x_mat <- as.matrix(x_scaled)

  # 2. Get Predicted Mean and Nuggets
  p_basic <- predict(object = model, x = x_mat)
  mu <- as.numeric(p_basic$mean)

  # 3. Build Full Covariance Matrix (The "By Hand" Logic)
  Kss <- cov_gen(X1 = x_mat, X2 = x_mat, theta = model$theta, type = model$covtype)
  Kxs <- cov_gen(X1 = model$X0, X2 = x_mat, theta = model$theta, type = model$covtype)

  # Latent uncertainty
  Sigma_latent <- model$nu_hat * (Kss - t(Kxs) %*% model$Ki %*% Kxs)

  # Total uncertainty (Latent + Nugget)
  Sigma_total <- Sigma_latent + diag(as.numeric(p_basic$nugs))
  Sigma_total <- (Sigma_total + t(Sigma_total)) / 2 # Symmetry fix

  # 4. Simulate Realizations (Null Distribution)
  # We use Cholesky decomposition to transform white noise into GP noise
  L <- chol(Sigma_total + diag(1e-8, nrow(Sigma_total))) # Small jitter for stability

  # Generate simulations: Mean + L * random_normals
  sims <- matrix(nrow = length(mu), ncol = n_sims)
  for(i in 1:n_sims) {
    sims[,i] <- mu + t(L) %*% rnorm(length(mu))
  }

  # 5. Extract Quantiles for the Ribbon (The "Correct" Plot Scale)
  df_plot <- data.frame(
    lon  = test_data$lon,
    mean = mu,
    ymin = apply(sims, 1, quantile, probs = 0.025),
    ymax = apply(sims, 1, quantile, probs = 0.975)
  )

  # 6. Plotting
  ggplot(df_plot, aes(x = lon)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "skyblue", alpha = 0.3) +
    geom_line(aes(y = ymax), color = "blue", linetype = "dashed") +
    geom_line(aes(y = ymin), color = "blue", linetype = "dashed") +
    geom_line(aes(y = mean), color = "red", linewidth = 1) +
    labs(title = "Simulated Total Uncertainty Envelope",
         x = "Longitude", y = "Latitudinal Residual") +
    theme_minimal()
}
#This should now run without the "non-numeric" error
test_data_grid <- seq(min(train_data$lon), max(train_data$lon), length.out = 1000)
plot_total_variance_proper(model_gp1, test_data_grid, all_data)

plot_total_variance_analytic <- function(model, lon_vector, training_data) {
  # 1. Scaling
  lon_min <- min(training_data$lon)
  lon_max <- max(training_data$lon)
  x_scaled <- (lon_vector - lon_min) / (lon_max - lon_min)
  x_mat <- as.matrix(x_scaled)
  
  # 2. Get Predictions
  p <- predict(object = model, x = x_mat)
  
  # 3. Calculate Components
  # v_latent: The uncertainty of the GP mean (the 'bloom' in empty spaces)
  # v_nugget: The estimated local noise
  v_latent <- as.numeric(p$sd2)
  v_nugget <- as.numeric(p$nugs)
  
  # The Red Line (Total Variance)
  v_total <- v_latent + v_nugget
  
  # 4. Analytic Quantiles (The "No-Simulation" Logic)
  # Instead of simulating, we use the latent SD as our Standard Error (SE).
  # We use the 97.5th percentile of a Normal distribution (1.96).
  se <- sqrt(v_latent)
  z_score <- qnorm(0.975) # This equals approximately 1.96
  
  df_plot <- data.frame(
    lon   = lon_vector,
    v_tot = v_total,
    # Symmetric Analytic Quantiles
    vmax  = v_total + (z_score * se),
    # Symmetric Subtraction with the hard zero floor
    vmin  = pmax(0, v_total - (z_score * se))
  )
  
  # 5. Visualization
  ggplot(df_plot, aes(x = lon)) +
    geom_ribbon(aes(ymin = vmin, ymax = vmax), 
                fill = "skyblue", alpha = 0.3) +
    geom_line(aes(y = vmax), color = "blue", linetype = "dashed", linewidth = 0.8) +
    geom_line(aes(y = vmin), color = "blue", linetype = "dashed", linewidth = 0.8) +
    geom_line(aes(y = v_tot), color = "red", linewidth = 1.2) +
    scale_y_continuous(limits = c(0, 2) ) +
    labs(
      x = "Longitude",
      y=NULL
      #y = "Variance", # (Degrees²)
    ) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank())
}


test_data_grid <- seq(min(train_data$lon), max(train_data$lon), length.out = 1000)
ninety_60<-plot_total_variance_analytic(model_gp1, test_data_grid, all_data)


ninety_70<-plot_total_variance_analytic(model_gp2, test_data_grid, all_data)


ninety_80<-plot_total_variance_analytic(model_gp3, test_data_grid, all_data)


ninety_83<-plot_total_variance_analytic(model_gp4, test_data_grid, all_data)


ninety_84<-plot_total_variance_analytic(model_gp5, test_data_grid, all_data)



plot_path <- 'C:/Users/steph/Downloads/Picture_Paper_5_22_2026'
if(!dir.exists(plot_path)) dir.create(plot_path)

ggsave(file.path(plot_path, "1960_sec_variance_plot.png"), plot = ninety_60, width = 6, height = 4, dpi = 300, bg = "white")
ggsave(file.path(plot_path,"1970_sec_variance_plot.png"), plot = ninety_70, width = 6, height = 4, dpi = 300, bg = "white")
ggsave(file.path(plot_path, "1980_sec_variance_plot.png"), plot = ninety_80, width = 6, height = 4, dpi = 300, bg = "white")
ggsave(file.path(plot_path, "1983_sec_variance_plot.png"), plot = ninety_83, width = 6, height = 4, dpi = 300, bg = "white")
ggsave(file.path(plot_path, "1984_sec_variance_plot.png"), plot = ninety_84, width = 6, height = 4, dpi = 300, bg = "white")

print("Analysis Complete. Plots saved.")










