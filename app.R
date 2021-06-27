library(tidyverse)
library(shinydashboard)
library(shinydashboardPlus)
library(tseries)

df_sum <- readRDS("estimation_sum.rds")
ggplot2::update_geom_defaults("line", list(size = 1.3))

# DGP functions -------------------------------------------------------------------------

simulate_bivariate_0ci <- function(t = 250, innov_sd = .5) {
  # Simulated bivariate non-cointegrated system
  tibble(
    x = cumsum(rnorm(t, 0, innov_sd)), 
    y = cumsum(rnorm(t, 0, innov_sd))
  )
}

simulate_bivariate_1ci <- function(t = 250, p = .75, innov_sd = .5) {
  # Simulated bivariate cointegrated system
  y <- cumsum(rnorm(t, 0, innov_sd))
  x <- y + arima.sim(list(ar = p), innov = rnorm(t, 0, innov_sd), n = t)
  tibble(x, y)
}

simulate_trivariate_0ci <- function(t = 250, innov_sd = .5) {
  # Simulated bivariate non-cointegrated system
  tibble(
    x = cumsum(rnorm(t, 0, innov_sd)), 
    y = cumsum(rnorm(t, 0, innov_sd)), 
    z = cumsum(rnorm(t, 0, innov_sd))
  )
}

simulate_trivariate_1ci <- function(t = 250, p = .75, innov_sd = .5) {
  # Simulated trivariate cointegrated system with 1 cointegrating vector
  y <- cumsum(rnorm(t, 0, innov_sd))
  z <- cumsum(rnorm(t, 0, innov_sd))
  x <- .5*y + .5*z + arima.sim(list(ar= p), innov = rnorm(t, 0, innov_sd), n = t)
  tibble(x, y, z)
}

simulate_trivariate_2ci <- function(t = 250, p = .75, innov_sd = .5) {
  # Simulated trivariate cointegrated system with 2 cointegrating vectors
  z <- cumsum(rnorm(t, 0, innov_sd))
  x <- z + arima.sim(list(ar= p), innov = rnorm(t, 0, innov_sd), n = t)
  y <- z + arima.sim(list(ar= p), innov = rnorm(t, 0, innov_sd), n = t)
  tibble(x, y, z)
}

# UI -----------------------------------------------------------------------------


ui <- dashboardPage(skin = "yellow",
                    dashboardHeader(title = "Johansen test simulations"),
                    dashboardSidebar(
                      sidebarMenu(
                        h3("DGP"),
                        sliderInput("n", "Number of series:", min = 2, max = 3, step = 1, value = 3),
                        uiOutput("ci"),
                        hr(),
                        br(), br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
                        br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
                        br(),br(),br(),br(),br(),br(),br(),br(),br(),
                        h6("Simulation results: 10 000 trajectories"),
                        h6("with each parameter combination"),
                        uiOutput("control"),
                        uiOutput("arg1"),
                        uiOutput("arg2")
                      )
                    ),
                    dashboardBody(
                      h1("Introduction of the data generating process to the simulation"),
                      box(width = 5,
                          title = "The data generating process",
                          withMathJax(),
                          shiny::uiOutput("dgp_formula")
                      ),
                      box(width = 5,
                          title = "DGP as R code",
                          shiny::verbatimTextOutput("dgp_rcode")
                      ),
                      box(width = 2,
                          title = "Public GitHub repository",
                          attachmentBlock(
                            image = "https://newsignature.com/wp-content/uploads/2020/08/github-logo-crop.png",
                            title = "PairsTrading",
                            href = "https://github.com/MarcellGranat/PairsTrading",
                            "MarcellGranat"
                          )
                      ),
                      box(width = 12,
                          title = "One unique trajectory",
                          shiny::inputPanel(
                            shinyWidgets::sliderTextInput("dgp_arg1", "Length of series:", as.character(1:20 * 50, animate = T), selected = as.character(250)),
                            shinyWidgets::sliderTextInput("dgp_arg2", "Standard deviation of innovation:",
                                                          as.character(c(.2, .5, 1, 1.5, 2)),
                                                          selected = as.character(0.5),
                                                          width = "500px"
                            ),
                            shiny::uiOutput("dgp_arg3")
                          ),
                          shiny::actionButton(inputId = "new_seed", "Generate new random values!"),
                          shiny::plotOutput("dgp_plot")
                      ),
                      shiny::fluidRow(
                        title = h1("Results of the simulation"),
                        box(height = "600px",
                            background = "maroon",
                            shiny::checkboxInput("bar_revert", label = "Revert scales"),
                            shiny::uiOutput("bar_arg"),
                            shiny::plotOutput("bar_chart")
                        ),
                        box(heigth = "700px",
                            title = "Heatmap", background = "blue", solidHeader = TRUE,
                            selectInput("estimated_ci", label = "Estimated rank?", 0:3, 1),
                            shiny::checkboxInput("hm_revert", label = "Revert scales"),
                            shiny::plotOutput("heatmap")
                        )
                      )
                    )
)

server <- function(input, output) {
  
  output$dgp_formula <- renderUI({
    
    if (input$n == 2 & input$ci == 0) {
      out <- withMathJax(helpText('$$x_{t} = x_{t-1} + u_{x,t}$$  \n
                         $$y_{t} = y_{t-1} + u_{y,t},$$ \n
                         where $$u \\sim norm(0, \\sigma)$$'))
    }
    
    if (input$n == 2 & input$ci == 1) {
      out <- withMathJax(helpText('$$x_{t} = \\beta y_t + u_{x,t}$$  \n
                         $$y_{t} = y_{t-1} + u_{y,t},$$ \n
                         where $$u \\sim norm(0, \\sigma)$$'))
    }
    
    if (input$n == 3 & input$ci == 0) {
      out <-  withMathJax(helpText('$$x_{t} = x_{t-1} + u_{x,t}$$  \n
                         $$y_{t} = y_{t-1} + u_{y,t}$$ \n
                         $$z_{t} = z_{t-1} + u_{z,t},$$ \n
                         where $$u \\sim norm(0, \\sigma)$$'))
    }
    
    if (input$n == 3 & input$ci == 1) {
      out <-  withMathJax(helpText('$$x_{t} = \\beta y_{t-1} + \\beta z_{t-1} + u_{x,t}$$  \n
                         $$y_{t} = y_{t-1} + u_{y,t}$$ \n
                         $$z_{t} = z_{t-1} + u_{z,t},$$ \n
                         where $$u \\sim norm(0, \\sigma)$$'))
    }
    
    if (input$n == 3 & input$ci == 2) {
      out <-  withMathJax(helpText('$$x_{t} =  \\beta z_{t-1} + u_{x,t}$$  \n
                         $$y_{t} = \\beta z_{t-1} + u_{y,t}$$ \n
                         $$z_{t} = z_{t-1} + u_{z,t},$$ \n
                         where $$u \\sim norm(0, \\sigma)$$'))
    }
    
    out
    
  })
  
  output$dgp_rcode <- shiny::renderText({
    dgp <- ifelse(input$n == 2, "bivariate", "trivariate") %>% 
      {str_c("simulate_", ., "_", input$ci, "ci")}
    
    f_text <- capture.output(get(dgp))
    
    out <- f_text %>% 
      head(-1) %>% 
      str_c("\n")
    
    out
  })
  
  output$ci <- shiny::renderUI(
    sliderInput("ci", "Number of cointegrating vectors:", min = 0, max = input$n - 1, step = 1, value = 0),
  )
  
  output$control <- shiny::renderUI(
    if (input$ci != 0) {
      selectInput("control", "Control for:", 
                  c("AR term - Length of series",
                    "AR term - Standard deviation of innovation",
                    "Length of series - Standard deviation of innovation"
                  )
      )
    }
  )
  
  setup <- reactive({
    df_sum <- df_sum %>% 
      ungroup() %>% 
      filter(n == input$n & ci == input$ci)
    
    if (input$ci == 0) {
      df <- select(df_sum, t, innov_sd, p, estimated_rank, nn, rate, reject_rate) 
    } else {
      if(input$control == "AR term - Length of series") {
        df <- select(df_sum, p, t, innov_sd, estimated_rank, nn, rate, reject_rate) %>% 
          filter(innov_sd == .5)
      }
      if(input$control == "AR term - Standard deviation of innovation") {
        df <-  select(df_sum, p, innov_sd, t, estimated_rank, nn, rate, reject_rate) %>% 
          filter(t == 250)
      }
      if(input$control == "Length of series - Standard deviation of innovation") {
        df <-  select(df_sum, t, innov_sd, p, estimated_rank, nn, rate, reject_rate) %>% 
          filter(p == .75)
      }
    }
    
    names(df) <-  c("x", "y", "left_out", "estimated_rank", "nn", "rate", "reject_rate")
    
    df
  }
  )
  
  output$dgp_arg3 <- shiny::renderUI({
    if (input$ci != 0) {
      shinyWidgets::sliderTextInput("dgp_arg3", "AR term:", as.character(c(1:9 / 10, .92, .94, .96, .98), selected = as.character(0.75)))
    }
  })
  
  seed_number <- shiny::reactive({
    input$new_seed
    floor(runif(1, max = 1000))
  })
  
  output$dgp_plot <- shiny::renderPlot({
    
    seed_number <- seed_number()
    set.seed(seed_number)
    dgp <- ifelse(input$n == 2, "bivariate", "trivariate") %>% 
      {str_c("simulate_", ., "_", input$ci, "ci")}
    
    if (input$ci == 0) {
      invoke(dgp, .x = list(t = input$dgp_arg1, innov_sd = input$dgp_arg2)) %>% 
        ts() %>% 
        forecast::autoplot() +
        labs(color = NULL) +
        theme(text = element_text(size = 20))
    } else {
      invoke(dgp, .x = list(t = input$dgp_arg1, innov_sd = input$dgp_arg2, p = input$dgp_arg3)) %>% 
        ts() %>% 
        forecast::autoplot() +
        labs(color = NULL) +
        theme(text = element_text(size = 20))
    }
    
    # ggplot(cars, aes(cars$speed, cars$dist)) + geom_point() 
  })
  
  output$heatmap <- shiny::renderPlot({
    
    df <- setup() %>%
      filter(estimated_rank == input$estimated_ci)
    
    if (input$ci == 0) {
      labels <- "Length of series - Standard deviation of innovation"
    } else {
      labels <- input$control
    }
    
    labels <- labels %>% 
      str_split(" - ") %>% 
      .[[1]]
    
    if(input$hm_revert) labels <- rev(labels)
    
    if (input$ci == 0) {
      left_out_text <- ""
    } else {
      if(input$control == "AR term - Length of series") left_out_text <- "Standard deviation of innovation = 0.5"
      if(input$control == "AR term - Standard deviation of innovation") left_out_text <- "Length of series = 250"
      if(input$control == "Length of series - Standard deviation of innovation") left_out_text <- "AR term = 0.75"
    }
    
    
    
    df <- crossing(df$x, df$y) %>%
      set_names("x", "y") %>%
      left_join(df) %>%
      replace_na(list(rate = 0))
    
    p <- ggplot(data = df)
    
    if(!input$hm_revert)  {
      p <- p + aes(x, y, fill = rate) 
    }
    
    if(input$hm_revert)  {
      p <- p + aes(y, x, fill = rate)
    }
    if (input$ci != 0) {
      p <- p + ggtitle(left_out_text)
    }
    
    p +
      labs(x = labels[1], y = labels[2]) +
      geom_tile(color = "black") +
      scale_fill_viridis_c(limits = c(0, 1), option = "E",
                           guide = guide_colorsteps())
  })
  
  
  output$bar_arg <- shiny::renderUI({
    
    if (input$ci == 0) {
      control <- "Length of series - Standard deviation of innovation"
    } else {
      control <- input$control
    }
    control <- str_split(control, " - ") %>% 
      .[[1]]
    if (input$bar_revert) {
      control <- control[2]
    } else {
      control <- control[1]
    }
    
    if (control == "Length of series") {
      out <- shinyWidgets::sliderTextInput("bar_arg", "Length of series:", as.character(1:20 * 50, animate = T))
    }
    if (control == "AR term") {
      out <- shinyWidgets::sliderTextInput("bar_arg", "AR term:", as.character(c(1:9 / 10, .92, .94, .96, .98), animate = T))
    }
    if (control == "Standard deviation of innovation") {
      out <- shinyWidgets::sliderTextInput("bar_arg", "Standard deviation of innovation:", as.character(c(.2, .5, 1, 1.5, 2), animate = T))
    }
    out
  }
  )
  
  
  
  
  output$bar_chart <- shiny::renderPlot({
    
    if (input$ci == 0) {
      control <- "Length of series - Standard deviation of innovation"
    } else {
      control <- input$control
    }
    control <- str_split(control, " - ") %>% 
      .[[1]]
    if (input$bar_revert) {
      x_lab <- control[1]
    } else {
      x_lab <- control[2]
    }
    
    df <- setup()
    if (input$bar_revert) {
      names(df)[1] <- "y" 
      names(df)[2] <- "x" 
    }
    df  %>%
      filter(x == as.numeric(input$bar_arg)) %>% 
      ggplot() +
      aes(x = y, y = rate, fill = factor(estimated_rank, levels = input$n:0)) +
      geom_col(color = "black") +
      scale_fill_viridis_d(option = "E") +
      labs(fill = "Estimated rank", x = x_lab, y = "Rate of estimation") + 
      theme(legend.position = "bottom")
    
  })
}

shinyApp(ui, server)