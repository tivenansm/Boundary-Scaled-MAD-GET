library(shiny)
library(MASS)
library(splines2)
library(plotly)
library(dplyr)
library(rsconnect)

# Deployment Credentials
rsconnect::setAccountInfo(name='tivenansm',
                          token='DF4E0D538EB3D450A5469EBD0F360526',
                          secret='ms+KzO8guSqhn1k07oXQIN7JtwfMs8cjBORzRHrb')

# ==============================================================================
# 2. UI (USER INTERFACE) 
# ==============================================================================
ui <- fluidPage(
  titlePanel("Global Envelope Test: Interactive Simulation (Pooled Variance)"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Underlying Truth"),
      sliderInput("coef", "Spline Coefficient 13 (Effect Size):", 
                  min = 1, max = 5, value = 1, step = 0.5),
      
      hr(),
      h4("Simulation Parameters"),
      sliderInput("n_reps", "Number of Monte Carlo Repetitions:", 
                  min = 100, max = 2000, value = 1000, step = 100),
      sliderInput("n_sims", "Null Distribution Size (per rep):", 
                  min = 100, max = 500, value = 200, step = 50),
      sliderInput("n_loc", "Number of Grid Points:", 
                  min = 50, max = 200, value = 100, step = 10),
      
      hr(),
      h4("Noise Settings"),
      sliderInput("noise_var1", "Baseline GP Variance (Period A):", 
                  min = 0.05, max = 1.0, value = 0.1, step = 0.05),
      sliderInput("noise_var2", "Modified GP Variance (Period B):", 
                  min = 0.05, max = 1.0, value = 0.1, step = 0.05),
      
      br(),
      actionButton("run_sim", "Update Dashboard", 
                   class = "btn-primary", style = "width: 100%; font-size: 16px; font-weight: bold;")
    ),
    
    mainPanel(
      plotlyOutput("dashboard", height = "800px")
    )
  )
)

# ==============================================================================
# 3. SERVER - The Brains and Math
# ==============================================================================
server <- function(input, output, session) {
  
  # Helper to simulate GP noise
  simulate_gp_residuals <- function(x, n_sim = 1, lengthscale = 5, noise_var = 0.1) {
    outer_diff <- outer(x, x, "-")
    K <- exp(- (outer_diff)^2 / (2 * lengthscale^2))
    Sigma <- noise_var * K
    MASS::mvrnorm(n = n_sim, mu = rep(0, length(x)), Sigma = Sigma)
  }
  
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = 'Running Rank Envelope simulations...', value = 0, {
      
      n_loc  <- input$n_loc
      n_reps <- input$n_reps
      n_sims <- input$n_sims
      COEF   <- input$coef
      # Grab individual variances
      NOISE1 <- input$noise_var1
      NOISE2 <- input$noise_var2
      
      x <- seq(0, 15, length.out = n_loc)
      b_spline <- splines2::bSpline(x, degree = 3, knots = c(2:13), Boundary.knots = c(0, 15))
      
      base_coef <- c(0, -1, -1, -1, -2, -2, -2, -2.5, -3, -2, -3, -3, 1, -1, -3)
      alt_coef  <- c(0, -1, -1, -1, -2, -2, -2, -2.5, -3, -2, -3, -3, COEF, -1, -3)
      
      # Determine if we are testing the Null Hypothesis
      use_null <- (COEF == 1)
      
      if (use_null) {
        coef1 <- base_coef
        coef2 <- base_coef
      } else {
        coef1 <- base_coef
        coef2 <- alt_coef
      }
      
      f1 <- as.numeric(b_spline %*% coef1 + 15)
      f2 <- as.numeric(b_spline %*% coef2 + 15)
      
      pvals_global <- numeric(n_reps)
      df_envelope <- NULL
      df_splines  <- NULL
      
      # The variance of the difference (A-B) is the sum of the variances
      pooled_noise <- NOISE1 + NOISE2
      
      for (i in 1:n_reps) {
        incProgress(1/n_reps)
        
        # --- Observed Curves with Different Noises ---
        y1 <- f1 + simulate_gp_residuals(x, 1, noise_var = NOISE1)
        y2 <- f2 + simulate_gp_residuals(x, 1, noise_var = NOISE2)
        obs_diff <- as.numeric(y1 - y2)
        
        # --- Null simulations using Pooled Variance ---
        # Representing Var(A - B) = Var(A) + Var(B)
        sim_diff <- simulate_gp_residuals(x, n_sims, noise_var = pooled_noise)
        
        # --- Rank Calculation ---
        all_curves <- rbind(obs_diff, sim_diff)
        T_all <- apply(abs(all_curves), 1, max)
        
        ranks <- rank(-T_all, ties.method = "random")
        
        r_obs <- ranks[1]
        r_sim <- ranks[-1]
        
        # --- p-value ---
        pvals_global[i] <- (sum(r_sim <= r_obs) + 1) / (n_sims + 1)
        
        # --- Capture state for plotting (first iteration only) ---
        if (i == 1) {
          df_splines <- rbind(
            data.frame(x = x, y = y1, Group = "Baseline (Period A)"),
            data.frame(x = x, y = y2, Group = "Modified (Period B)")
          )
          
          # 95% Rank envelope bounds
          k_alpha <- floor(0.05 * (n_sims + 1))
          lower <- apply(sim_diff, 2, function(col) sort(col)[k_alpha])
          upper <- apply(sim_diff, 2, function(col) sort(col)[n_sims - k_alpha])
          
          df_envelope <- data.frame(
            x = x,
            Obs_Diff = obs_diff,
            Lower_Bound = lower,
            Upper_Bound = upper
          )
        }
      }
      
      list(
        splines = df_splines,
        envelope = df_envelope,
        pvals = pvals_global
      )
    })
  }, ignoreNULL = FALSE)
  
  output$dashboard <- renderPlotly({
    d <- sim_data()
    
    # Plot 1: The Raw Simulated Boundaries
    p1 <- plotly::plot_ly(d$splines, x = ~x, y = ~y, color = ~Group,
                          colors = c("black", "red"),
                          type = 'scatter', mode = 'lines') %>%
      plotly::layout(yaxis = list(title = "Lat"))
    
    # Plot 2: The Difference and the 95% Global Envelope
    p2 <- plotly::plot_ly(d$envelope, x = ~x) %>%
      plotly::add_ribbons(ymin = ~Lower_Bound, ymax = ~Upper_Bound,
                          fillcolor = "rgba(200,200,200,0.5)",
                          line = list(color = "transparent"),
                          name = "95% Rank Envelope") %>%
      plotly::add_lines(y = ~Obs_Diff, line = list(color = "blue"),
                        name = "Observed Difference") %>%
      plotly::layout(yaxis = list(title = "Difference"))
    
    # Plot 3: Distribution of P-values over MC Repetitions
    p3 <- plotly::plot_ly(x = ~d$pvals, type = "histogram", 
                          histnorm = "probability density", 
                          xbins = list(start = 0, end = 1, size = 0.05),
                          marker = list(color = "lightblue", line = list(color = "white", width = 1)),
                          showlegend = FALSE) %>%
      plotly::layout(
        yaxis = list(title = "Density"), 
        xaxis = list(title = "P-value", range = c(0, 1)),
        shapes = list(
          list(type = "line", x0 = 0, x1 = 1, y0 = 1, y1 = 1,
               xref = "paper", line = list(dash = "dot"))
        )
      )
    
    subplot(p1, p2, p3, nrows = 3, margin = 0.07, titleY = TRUE)
  })
}

shinyApp(ui = ui, server = server)