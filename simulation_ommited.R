source("functions.R")

df <- crossing(
  t = 1:20 * 50,
  p = c(1:9 / 10, .92, .94, .96, .98),
  trajectory = 1:1000
) 

evaluated_df <- tibble()

for (i in 1:100) {
  current_df <- df %>% 
    filter(cut(row_number(), 100, F) == i) %>% 
    mutate(
      dgp = pmap(.l = list(t = t, p = p), .f = simulate_trivariate_1ci),
      dgp_xy = map(dgp, ~ select(., x, y)),
      dgp_xz = map(dgp, ~ select(., x, z)),
      dgp_yz = map(dgp, ~ select(., y, z)),
    ) %>% 
    pivot_longer(cols = dgp:dgp_yz, names_to = "contained", values_to = "dgp") %>% 
    mutate(
      rank = map_dbl(dgp, ~ {tsDyn::VECM(., lag = 0, estim = "ML", include = "none") %>%  
          tsDyn::rank.test(type = 'trace', cval = 0.05) %>% 
          .$r}
      )
    ) %>% 
    select(-dgp) %>% 
    pivot_wider(names_from = contained, values_from = rank) %>% 
    rename(rank_xyz = dgp) %>% 
    rename_all(~ str_replace_all(., "dgp", "rank")) %>% 
    mutate(rank_xy_xz = rank_xy + rank_xz)
  evaluated_df <- bind_rows(evaluated_df, current_df)
  print(i)
}

df_sum <- evaluated_df %>%
  count(t, p, rank_xyz, rank_xy_xz)

write_rds(df_sum, file = "C:/rprojects/PairsTrading/estimation_omitted_sum.rds")