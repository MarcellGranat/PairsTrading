library(tidyverse)
library(shinydashboard)
library(shinydashboardPlus)
library(tseries)

theme_set(theme(text = element_text(size = 20)))

df_sum <- readRDS("estimation_sum.rds")
df_omitted <- readRDS("estimation_omitted_sum.rds")
df_finance <- readRDS("estimation_financial.rds")
load("data.RData")
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
                    dashboardHeader(title = "ÚNKP"),
                    dashboardSidebar(
                      sidebarMenu(
                        h3("DGP"),
                        sliderInput("n", "Number of series:", min = 2, max = 3, step = 1, value = 3),
                        uiOutput("ci"),
                        menuItem("Introduction", tabName = "intro", icon = icon("play-circle")),
                        menuItem("Simulation 1.", icon = icon("th"), tabName = "sim1"),
                        menuItem("Simulation 2.", tabName = "sim2", icon = icon("project-diagram")),
                        menuItem("Real-life performance", tabName = "finance", icon = icon("hand-holding-usd")),
                        uiOutput("control"),
                        uiOutput("arg1"),
                        uiOutput("arg2")
                      )
                    ),
                    dashboardBody(
                      tabItems(
                        tabItem(tabName = "intro",
                                h1("Introduction"),
                                h2("The data generating process to the simulation"),
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
                                box(width = 12, height = "1000px",
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
                                    shiny::plotOutput("dgp_plot", height = "700px")
                                )
                        )
                        ,
                        
                        tabItem(tabName = "sim1",
                                h1("Results of the simulation"),
                                h2("Simulation results: 10 000 trajectories with each parameter combination"),
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
                        ),
                        tabItem(tabName = "sim2", width = 12,
                                h1("Results of the simulation 2."),
                                h2("Comparison of tests performed on 2 or 3 time-series"),
                                box(width = 12, heigth = "800px",
                                  shiny::inputPanel(
                                    shinyWidgets::sliderTextInput("sim2_t", "Length of series:", as.character(1:20 * 50, animate = T), selected = as.character(250)),
                                    shinyWidgets::sliderTextInput("sim2_p", "AR term:", as.character(c(1:9 / 10, .92, .94, .96, .98)), selected = as.character(0.7))
                                  ),
                                  shiny::plotOutput("sim2_plot")
                                )
                        ),
                        tabItem(tabName = "finance", width = 12,
                                h1("Testing on real data"),
                                box(width = 12, shiny::plotOutput("bankdata_plot")),
                                box(
                                  shiny::inputPanel(
                                    shinyWidgets::sliderTextInput("finance_window_length", "Length of time-window:", as.character(1:20 * 50), selected = as.character(250))
                                  ),
                                  shiny::plotOutput("finance_tsplot")
                                ),
                                box(
                                  title = "Rate of rank remaining by window-to-window",
                                  width = 3, 
                                  shinydashboard::valueBoxOutput("finance_persistency_3", width = "300px"),
                                  shinydashboard::valueBoxOutput("finance_persistency_2", width = "300px")
                                ),
                                
                                    shiny::uiOutput("img")
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
  
  output$sim2_plot <- shiny::renderPlot({
    df <- df_omitted %>% 
      filter(t == input$sim2_t, p == input$sim2_p) %>% 
      mutate(n = n/sum(n)) 
    
    crossing(df$rank_xyz, df$rank_xy_xz) %>% 
      set_names("rank_xyz", "rank_xy_xz") %>% 
      left_join(df) %>% 
      replace_na(list(n = 0)) %>% 
      ggplot() +
      aes(x = rank_xyz, y = rank_xy_xz, fill = n) +
      geom_tile(color = "black") + 
      scale_fill_viridis_c(limits = c(0, 1), option = "E",
                           guide = guide_colorsteps()) +
      labs(x = "Estimated rank of trivariate test (x-y-z)", y = "Sum of the estimated ranks from\nthe bivariate tests (x-y, x-z)")
  })
  
  output$bankdata_plot <- shiny::renderPlot({
    Bankdata %>% 
      pivot_longer(-1) %>% 
      ggplot() +
      aes(Date, value) +
      geom_line() +
      facet_wrap(vars(name), scales = "free_y") +
      labs("The investigated stocks")
  })
  
  output$finance_tsplot <- shiny::renderPlot({
    df_finance %>% 
      filter(window_length == input$finance_window_length) %>% 
      count(date, rank_xyz) %>% 
      mutate(
        rank_xyz = factor(rank_xyz, levels = as.character(3:0), ordered = TRUE),
      ) %>% 
      group_by(date) %>%
      group_modify(~ mutate(.x, n = n / sum(n))) %>%
      ungroup() %>%
      ggplot(aes(x = date, y = n, fill = rank_xyz)) +
      geom_col(color = "black") +
      scale_fill_viridis_d(option = "E", direction = -1) + 
      scale_y_continuous(labels = scales::percent) +
      labs(y = "Rate of estimated rank \nfrom trivariate tests", fill = NULL)
  })
  
  finance_lags <- reactive(
    df_finance %>% 
      filter(window_length == input$finance_window_length) %>% 
      group_by(stock1, stock2, stock3) %>% 
      group_modify(
        ~ mutate(.x, 
                 rank_xyz_l = lag(rank_xyz),
                 rank_xy_l = lag(rank_xy),
                 rank_xz_l = lag(rank_xz),
                 rank_yz_l = lag(rank_yz),
        )
      ) %>% 
      ungroup()
    
  )
  
  output$finance_persistency_3 <- shinydashboard::renderValueBox({
    finance_lags() %>% 
      summarise(
        rank_xyz = mean(rank_xyz == rank_xyz_l, na.rm = TRUE)
      ) %>% 
      pull(rank_xyz) %>% 
      scales::percent(accuracy = .01) %>% 
      first() %>% 
      {
        shinydashboard::valueBox("Trivariate tests", ., icon = icon("thums-up"))
      }
  })
  
  output$finance_persistency_2 <- shinydashboard::renderValueBox({
    finance_lags() %>% 
      summarise(
        rank_xy = mean(rank_xy == rank_xy_l, na.rm = TRUE),
        rank_xz = mean(rank_xz == rank_xz_l, na.rm = TRUE),
        rank_yz = mean(rank_yz == rank_yz_l, na.rm = TRUE)
      ) %>% 
      mutate(m = mean(rank_xy, rank_xz, rank_yz)) %>% 
      pull(m) %>% 
      scales::percent(accuracy = .01) %>% 
      first() %>% 
      {
        shinydashboard::valueBox("Bivariate tests", ., icon = icon("thums-up"))
      }
  })
  
  output$img <- renderUI({
    tags$img(src = "https://i.ytimg.com/vi/KYGwd0Y3-zk/hqdefault.jpg")
  })
  
  df_finance %>% 
    mutate_at(7:9, ~ . == 1) %>% 
    mutate(rank_sum = rank_xy + rank_xz + rank_yz) %>% 
    count(rank_xyz, rank_sum)
  
  df_finance %>% 
    filter(window_length == 250) %>% 
    group_by(stock1, stock2, stock3) %>% 
    group_modify(
      ~ mutate(.x, 
               rank_xyz_l = lag(rank_xyz),
               rank_xy_l = lag(rank_xy),
               rank_xz_l = lag(rank_xz),
               rank_yz_l = lag(rank_yz),
      )
    ) %>% 
    ungroup() %>% 
    {
      bind_rows(
        set_names(select(., rank_xyz_l, rank_xy, rank_xy_l), "rank_xyz", "rank", "rank_l"),
        set_names(select(., rank_xyz_l, rank_xz, rank_xz_l), "rank_xyz", "rank", "rank_l"),
        set_names(select(., rank_xyz_l, rank_yz, rank_yz_l), "rank_xyz", "rank", "rank_l")
      )
    } %>% 
    na.omit() %>% 
    count(rank_xyz, rank, rank_l) %>% 
    filter(rank == 1) %>% 
    ggplot() +
    aes(rank_xyz, rank_l, fill = n) %>% 
      geom_tile()
}

shinyApp(ui, server)