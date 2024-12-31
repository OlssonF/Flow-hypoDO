#' @param DO a long format dataframe of DO profiles
#' @param threshold what is the threshold for anoxia? Default 1
#' @param depth_use what depth should be used to evaluate anoxia? Default 5
calc_anoxic_period <- function(DO, 
                               threshold = 1,
                               depth_use = 5) {
  
  # DO <- daily_DO
  
  DO <- DO |> 
    filter(depth == depth_use) |> 
    mutate(greater_than_1 = ifelse(DO_mgl_daily >= threshold, T, F),
           greater_than_1 = ifelse(is.na(greater_than_1),
                                   T, greater_than_1))
  
  rle_DO_5 <- rle(DO$greater_than_1 == T)
  
  rle_dates <- data.frame(exceed_1 = rle_DO_5$values,
                          length = rle_DO_5$lengths)  |>
    
    mutate(exceed_1 = ifelse(is.na(exceed_1),
                             T, exceed_1),
           end_row = cumsum(length),
           end_date = DO$datetime[end_row],
           # Get the start of each run
           start_row = end_row - length + 1,
           start_date = DO$datetime[start_row],
           year = year(end_date)) 
  
  #extract the start and end dates for each year periods
  
  rle_dates |> 
    group_by(year, exceed_1) |> 
    summarise(start = min(end_date),
              end = max(end_date), 
              .groups = 'drop') |> 
    pivot_longer(start:end, 
                 names_to = 'period',
                 values_to = 'datetime') |> 
    filter(exceed_1 == T & period == 'start' | # last day of oxic conditions
             exceed_1 == F & period == 'end') |>  # last day of anoxic conditions 
    select(year, datetime, period) |> 
    pivot_wider(values_from = datetime, names_from = period)

}

