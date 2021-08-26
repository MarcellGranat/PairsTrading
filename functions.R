library(tidyverse)
library(tsDyn)
library(tseries)

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
