summer_anoxic_periods <- read_csv('output/all_data.csv') 


# Results stats ---------------------
# frequency of reoxygenation
summer_anoxic_periods |> 
  group_by(year(datetime)) |>
  summarise(n = sum(DO_5 > 1), tot = n()) |>
  mutate(n/tot)

for (year_use in c(2018, 2019)) {
  
  period_use <- summer_anoxic_periods |> 
    filter(year == year_use) 
  
  rle_use <- rle(period_use$DO_5 > 1)
  
  data.frame(length = rle_use$lengths,
             oxic = rle_use$values) |> 
    mutate(end = cumsum(length),
           start = end - length + 1,
           date_start = period_use$datetime[start],
           date_end = period_use$datetime[end]) |> 
    filter(oxic == T) |> 
    slice_max(length) |> 
    print()
}

         