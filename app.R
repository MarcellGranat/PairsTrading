library(tidyverse)
library(shinydashboard)

df <- readRDS("C:/rprojects/PairsTrading/estimation.rds")
df_sum <- df %>% 
  count(n, ci, innov_sd, t, p, estimated_rank) %>% 
  arrange(desc(estimated_rank)) %>% 
  group_by(n, ci, innov_sd, t, p) %>% 
  group_modify(~ mutate(.x, rate = nn / sum(nn), reject_rate = cumsum(rate) - rate))

ui <- dashboardPage(skin = "yellow",
                    dashboardHeader(title = "Johansen test simulations"),
                    dashboardSidebar(
                      sidebarMenu(
                        h3("DGP"),
                        sliderInput("n", "Number of series:", min = 2, max = 3, step = 1, value = 3),
                        uiOutput("ci"),
                        uiOutput("control"),
                        uiOutput("arg1"),
                        uiOutput("arg2")
                      )
                    ),
                    dashboardBody(
                      box(
                        title = "Heatmap", background = "blue", solidHeader = TRUE,
                          selectInput("estimated_ci", label = "Estimated rank?", 0:3, 1),
                          shiny::checkboxInput("hm_revert", label = "Revert scales"),
                        shiny::plotOutput("heatmap")
                      )
                    )
)

server <- function(input, output) {
  output$ci <- shiny::renderUI(
    sliderInput("ci", "Number of ci:", min = 0, max = input$n - 1, step = 1, value = 0),
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
  
  
  # output$arg1 <- shiny::renderUI(
  #   if (input$ci != 0) {
  #     if (str_starts(input$control, "Length of series")) {
  #       shinyWidgets::sliderTextInput("arg1", "Length of series:", as.character(1:20 * 50, animate = T))
  #     } else {
  #       shinyWidgets::sliderTextInput("arg1", "AR term:", as.character(c(1:9 / 10, .92, .94, .96, .98), animate = T))
  #     }
  #   } else {
  #     shinyWidgets::sliderTextInput("arg1", "Length of series:", as.character(1:20 * 50, animate = T))
  #   }
  # )
  # 
  # output$arg2 <- shiny::renderUI(
  #   if (input$ci != 0) {
  #     if (str_ends(input$control, "Standard deviation of innovation")) {
  #       shinyWidgets::sliderTextInput("arg2", "Standard deviation of innovation:", as.character(c(.2, .5, 1, 1.5, 2), animate = T))
  #     } else {
  #       shinyWidgets::sliderTextInput("arg2", "Length of series:", as.character(1:20 * 50, animate = T))
  #     }
  #   } else {
  #     shinyWidgets::sliderTextInput("arg2", "Standard deviation of innovation:", as.character(c(.2, .5, 1, 1.5, 2), animate = T))
  #   }
  #   
  # )
  
  
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
      scale_fill_gradient(limits = c(0, 1))
  })
}

shinyApp(ui, server)