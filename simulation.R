# Setup -----------------------------------------------------------------------------

library(tidyverse)
library(tsDyn)
set.seed(2021)

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
  y2 <- cumsum(rnorm(t, 0, innov_sd))
  y1 <- y2 + arima.sim(list(ar= p), innov = rnorm(t, 0, innov_sd), n = t)
  tibble(x = y1, y = y2)
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
  y2 <- cumsum(rnorm(t, 0, innov_sd))
  y3 <- cumsum(rnorm(t, 0, innov_sd))
  y1 <- .5*y2 + .5*y3 + arima.sim(list(ar= p), innov = rnorm(t, 0, innov_sd), n = t)
  tibble(x = y1, y = y2, z = y3)
}

simulate_trivariate_2ci <- function(t = 250, p = .75, innov_sd = .5) {
  # Simulated trivariate cointegrated system with 2 cointegrating vectors
  y3 <- cumsum(rnorm(t, 0, innov_sd))
  y1 <- y3 + arima.sim(list(ar= p), innov = rnorm(t, 0, innov_sd), n = t)
  y2 <- y3 + arima.sim(list(ar= p), innov = rnorm(t, 0, innov_sd), n = t)
  tibble(x = y1, y = y2, z = y3)
}



df <- expand.grid(
  n = 2:3,
  ci = 0:2,
  control_1 = c("t", "p", "innov_sd"),
  control_2 = c("t", "p", "innov_sd"), stringsAsFactors = FALSE
) %>% 
  filter(control_1 < control_2 & 
           ci < n &
           (ci != 0 | (control_1 != "p" & control_2 != "p"))
  ) %>% 
  arrange(n, ci) %>% 
  tibble()

control_values <- list(
  t = 1:20 * 50,
  p = c(1:9 / 10, .92, .94, .96, .98),
  innov_sd = c(.2, .5, 1, 1.5, 2)
)

df <- df %>% 
  mutate(
    fixed = map2_chr(control_1, control_2, ~ setdiff(c("t", "p", "innov_sd"), c(.x, .y))),
    control_value_1 = map(control_1, ~ control_values[[.]]),
    control_value_2 = map(control_2, ~ control_values[[.]]),
  ) %>% 
  unnest(control_value_1) %>% 
  unnest(control_value_2) 

df <- df %>% 
  mutate(
    f = str_c(ifelse(n == 2, "simulate_bivariate_", "simulate_trivariate_"), ci, "ci"),
    innov_sd = ifelse(control_1 == "innov_sd", control_value_1, .5),
    t = ifelse(control_2 == "t", control_value_2, 250),
    p = ifelse(control_1 == "p", control_value_1, .75),
    p = ifelse(control_2 == "p", control_value_2, p)
  ) %>% 
  select(n, ci, innov_sd, t, p, f) %>% 
  slice(rep(1:n(), each = 10000))

simulate_and_estimate_rank <- function(f, innov_sd, t, p) {
  if (f != "simulate_bivariate_0ci" & f != "simulate_trivariate_0ci") {
    invoke(f, list(innov_sd = innov_sd, t = t, p = p)) %>% 
      tsDyn::VECM(lag = 0, estim = "ML", include = "none") %>%  
      tsDyn::rank.test(type = 'trace', cval = 0.05) %>% 
      .$r
  } else {
    invoke(f, list(innov_sd = innov_sd, t = t)) %>% 
      tsDyn::VECM(lag = 0, estim = "ML", include = "none") %>%  
      tsDyn::rank.test(type = 'trace', cval = 0.05) %>% 
      .$r
  }
}




library(parallel)
library(progressr)
cl <- makeCluster(8)
clusterEvalQ(cl, library(tidyverse))
clusterExport(cl, list("simulate_bivariate_0ci", "simulate_bivariate_1ci", "simulate_trivariate_0ci", "simulate_trivariate_1ci", 
                       "simulate_trivariate_2ci", "df", "simulate_and_estimate_rank"), envir = environment())
evaluated_df <- tibble()

for (i in 1:1000) {
  
  current_df <- df %>% 
    filter(cut(row_number(), 1000, F) == i)
  
  estimated_rank <- parApply(cl, current_df, 1, FUN = function(x) {
    simulate_and_estimate_rank(f = x[6], as.numeric(x[3]), as.numeric(x[4]), as.numeric(x[5]))
  })
  
  current_df$estimated_rank <- estimated_rank
  evaluated_df <- bind_rows(evaluated_df, current_df)
  
  if (i %% 10 == 0) {
    cat(str_c(str_c(rep("=", i/1000*100 - 1), collapse = ""), ">", str_c(rep(" ", 100-i/1000*100), collapse = ""), "| ", i/1000*100, "%\n", collapse = ""))
  }
}

write_rds(evaluated_df, "C:/rprojects/PairsTrading/estimation.rds")

df_sum <- evaluated_df %>% 
  count(n, ci, innov_sd, t, p, estimated_rank) %>% 
  arrange(desc(estimated_rank)) %>% 
  group_by(n, ci, innov_sd, t, p) %>% 
  group_modify(~ mutate(.x, rate = nn / sum(nn), reject_rate = cumsum(rate) - rate)) %>% 
  ungroup()

write_rds(df_sum, file = "C:/rprojects/PairsTrading/estimation_sum.rds")

