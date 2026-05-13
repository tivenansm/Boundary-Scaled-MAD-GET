library(shiny)
library(MASS)
library(splines2)
library(plotly)
library(dplyr)

ui <- fluidPage(
  titlePanel("Interactive Simulation for Global Scaled MAD Envelope Test"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Coeficient Slider"),
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

server <- function(input, output, session) {
  
  simulate_gp_residuals <- function(x, n_sim = 1, lengthscale = 5, noise_var = 0.1) {
    outer_diff <- outer(x, x, "-")
    K <- exp(- (outer_diff)^2 / (2 * lengthscale^2))
    Sigma <- noise_var * K
    MASS::mvrnorm(n = n_sim, mu = rep(0, length(x)), Sigma = Sigma)
  }
  
  sim_data <- eventReactive(input$run_sim, {
    # Renamed progress message to MAD-GET
    withProgress(message = 'Calculating Scaled MAD-GET...', value = 0, {
      
      n_loc  <- input$n_loc
      n_reps <- input$n_reps
      n_sims <- input$n_sims
      COEF   <- input$coef
      NOISE1 <- input$noise_var1
      NOISE2 <- input$noise_var2
      
      x <- seq(0, 15, length.out = n_loc)
      b_spline <- splines2::bSpline(x, degree = 3, knots = c(2:13), Boundary.knots = c(0, 15))
      
      base_coef <- c(0, -1, -1, -1, -2, -2, -2, -2.5, -3, -2, -3, -3, 1, -1, -3)
      alt_coef  <- c(0, -1, -1, -1, -2, -2, -2, -2.5, -3, -2, -3, -3, COEF, -1, -3)
      
      coef1 <- base_coef
      coef2 <- if(COEF == 1) base_coef else alt_coef
      
      f1 <- as.numeric(b_spline %*% coef1 + 15)
      f2 <- as.numeric(b_spline %*% coef2 + 15)
      
      pvals_global <- numeric(n_reps)
      df_envelope <- NULL
      df_splines  <- NULL
      df_resids   <- NULL # For the dark blue lines
      pooled_noise <- NOISE1 + NOISE2
      
      for (i in 1:n_reps) {
        incProgress(1/n_reps)
        y1 <- f1 + simulate_gp_residuals(x, 1, noise_var = NOISE1)
        y2 <- f2 + simulate_gp_residuals(x, 1, noise_var = NOISE2)
        obs_diff <- as.numeric(y1 - y2)
        sim_diff <- simulate_gp_residuals(x, n_sims, noise_var = pooled_noise)
        
        all_curves <- rbind(obs_diff, sim_diff)
        T_all <- apply(abs(all_curves), 1, max)
        ranks <- rank(-T_all, ties.method = "random")
        pvals_global[i] <- (sum(ranks[-1] <= ranks[1]) + 1) / (n_sims + 1)
        
        if (i == 1) {
          df_splines <- rbind(
            data.frame(x = x, y = y1, Group = "Baseline (Period A)"),
            data.frame(x = x, y = y2, Group = "Modified (Period B)")
          )
          
          # Capture first 10 null simulations as "Dark Blue Residuals"
          resids_to_plot <- sim_diff[1:min(10, n_sims), ]
          df_resids <- data.frame(
            x = rep(x, each = nrow(resids_to_plot)),
            y = as.numeric(t(resids_to_plot)),
            id = rep(1:nrow(resids_to_plot), times = length(x))
          )
          
          k_alpha <- floor(0.05 * (n_sims + 1))
          lower <- apply(sim_diff, 2, function(col) sort(col)[k_alpha])
          upper <- apply(sim_diff, 2, function(col) sort(col)[n_sims - k_alpha])
          df_envelope <- data.frame(x = x, Obs_Diff = obs_diff, Lower_Bound = lower, Upper_Bound = upper)
        }
      }
      list(splines = df_splines, envelope = df_envelope, pvals = pvals_global, resids = df_resids)
    })
  }, ignoreNULL = FALSE)
  
  output$dashboard <- renderPlotly({
    d <- sim_data()
    
    # Plot 1: Standard comparison
    p1 <- plotly::plot_ly(d$splines, x = ~x, y = ~y, color = ~Group, colors = c("black", "green"), type = 'scatter', mode = 'lines') %>%
      plotly::layout(yaxis = list(title = "Values"))
    
    # Plot 2: MAD-GET Difference
    p2 <- plotly::plot_ly(data = d$envelope, x = ~x) %>%
      
      # 1. The Ribbon (Skyblue with 0.25 alpha)
      plotly::add_ribbons(ymin = ~Lower_Bound, ymax = ~Upper_Bound, 
                          fillcolor = "rgba(135, 206, 235, 0.25)", 
                          line = list(color = "transparent"), 
                          name = "Simulations",            
                          showlegend = FALSE) %>%        
      
      # 2. Upper Edge (Blue Dashed)
      plotly::add_lines(y = ~Upper_Bound, 
                        line = list(color = "blue", dash = "dash", width = 1.5), 
                        showlegend = TRUE, name = "95% GET") %>%
      
      # 3. Lower Edge (Blue Dashed) - Set back to FALSE
      plotly::add_lines(y = ~Lower_Bound, 
                        line = list(color = "blue", dash = "dash", width = 1.5), 
                        showlegend = FALSE, name = "95% GET") %>%
      
      # 4. The Main Observed Difference (Red)
      plotly::add_lines(y = ~Obs_Diff, 
                        line = list(color = "red", width = 2), 
                        name = "Observed Diff",       
                        showlegend = TRUE) %>%
      
      # 5. Layout
      plotly::layout(
        showlegend = TRUE,            
        plot_bgcolor = "white",       
        paper_bgcolor = "white",
        xaxis = list(
          title = "",                 
          showgrid = TRUE,
          gridcolor = "#E5E5E5",      
          zeroline = FALSE
        ),
        yaxis = list(
          title = "Deviants",                 
          showgrid = TRUE,
          gridcolor = "#E5E5E5",
          zeroline = FALSE
        )
      )
    
    # Plot 3: Histogram (Silent Legend)
    p3 <- plotly::plot_ly(x = ~d$pvals, type = "histogram", histnorm = "probability density", 
                          xbins = list(start = 0, end = 1, size = 0.05),
                          marker = list(color = "lightblue", line = list(color = "white", width = 1)),
                          name = "", showlegend = FALSE) %>%
      plotly::layout(yaxis = list(title = "Probability Density"), xaxis = list(title = "P-value", range = c(0, 1)),
                     shapes = list(list(type = "line", x0 = 0, x1 = 1, y0 = 1, y1 = 1, xref = "paper", line = list(dash = "dot"))))
    
    # Combine subplots
    # Removed the plotly::style() overrides that were hiding the legends!
    subplot(p1, p2, p3, nrows = 3, margin = 0.07, titleY = TRUE)     
  })
}

shinyApp(ui = ui, server = server)