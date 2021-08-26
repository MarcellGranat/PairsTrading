source("functions.R")

load("data.Rdata")

df <- Bankdata

v <- df %>% 
  names() %>% 
  .[-1]

evaluated_df <- tibble()

for (window_length in seq(from = 50, to = 1000, by = 50)) {
  
  current_df <- crossing(stock1 = v, stock2 = v, stock3 = v) %>% 
    filter(stock1 < stock2, stock2 < stock3) %>% 
    crossing(t = seq(nrow(df) %/% window_length)) %>% 
    mutate(stocks = pmap(list(stock1, stock2, stock3, t), .f = function(x, y, z, t) {
      select(df, x, y, z) %>% 
        head(-(nrow(df) %% window_length)) %>% 
        filter(cut(row_number(), nrow(df) %/% window_length, FALSE) == t)
    })
    ) %>% 
    mutate(
      rank_xyz = map_dbl(stocks, ~ {tsDyn::VECM(., lag = 0, estim = "ML", include = "none") %>%  
          tsDyn::rank.test(type = 'trace', cval = 0.05) %>% 
          .$r}
      ),
      rank_xy = map_dbl(stocks, ~ {
        .[, 1:2] %>% 
          tsDyn::VECM(lag = 0, estim = "ML", include = "none") %>%  
          tsDyn::rank.test(type = 'trace', cval = 0.05) %>% 
          .$r}
      ),
      rank_xz = map_dbl(stocks, ~ {
        .[, c(1, 3)] %>% 
          tsDyn::VECM(lag = 0, estim = "ML", include = "none") %>%  
          tsDyn::rank.test(type = 'trace', cval = 0.05) %>% 
          .$r}
      ),
      rank_yz = map_dbl(stocks, ~ {
        .[, 2:3] %>% 
          tsDyn::VECM(lag = 0, estim = "ML", include = "none") %>%  
          tsDyn::rank.test(type = 'trace', cval = 0.05) %>% 
          .$r}
      ),
    ) %>% 
    mutate(window_length = window_length) %>% 
    select(window_length, everything(), -stocks)
  
  evaluated_df <- bind_rows(evaluated_df, current_df)
  print(window_length)
}


evaluated_date <- tibble()
for (window_length in seq(from = 50, to = 1000, by = 50)) {
  
current_date <- df %>% 
  head(-(nrow(df) %% window_length)) %>% 
  mutate(t = cut(row_number(), nrow(df) %/% window_length, FALSE)) %>% 
  group_by(t) %>% 
  group_map(~ {head(.x, nrow(.x) / 2) %>% 
      pull(Date) %>% 
      last()
  }
  ) %>% 
  reduce(c) %>% 
  enframe() %>% 
  set_names("t", "date") %>% 
  mutate(
    window_length = window_length,
    date = lubridate::ymd(date)
    )
evaluated_date <- bind_rows(evaluated_date, current_date)
}

evaluated_df <- evaluated_df %>% 
  left_join(evaluated_date)

write_rds(evaluated_df, file = "C:/rprojects/PairsTrading/estimation_financial.rds")
