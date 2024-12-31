
calc_strat_dates <- function(density_diff = 0.1,
                             use_depths = c(1,6),
                             wtr) {
  
  all_dates <- data.frame(datetime = seq.Date(min(wtr$datetime),
                                              max(wtr$datetime),
                                              by = 'day'))
  ## extract the depths that will be used to calculate the density difference (surface, bottom)
  wtr_use <- wtr |>
    full_join(all_dates, by = 'datetime') |> 
    pivot_longer(-datetime, 
                 values_to = 'temp',
                 names_to = 'depth', 
                 names_prefix = 'wtr_', names_transform = list(depth = as.numeric)) |> 
    dplyr::filter(depth %in% use_depths) |> 
    dplyr::mutate(location = ifelse(depth == min(use_depths), 'top', 'bottom'))
  
  
  
  density_obs <-
  wtr_use |> 
    mutate(density = rLakeAnalyzer::water.density(temp)) |> 
    # select(datetime, density, temp, depth) |> 
    pivot_wider(values_from = c(density, temp), 
                names_from = location, id_cols = c(datetime)) |> 
    mutate(dens_diff = density_bottom - density_top,
           strat = ifelse(abs(dens_diff > 0.1) & temp_top > temp_bottom, 1, 0),
           strat = imputeTS::na_locf(strat, maxgap = 3)) # in case of a gap the last obs is carried forward
  
  
  # extract the dates of the stratified periods
  #using a loop function to go through each year and do the rle function
  
  strat <- data.frame(year = unique(year(density_obs$datetime)), 
                      length = NA,
                      start = NA,
                      end = NA)
  
  for (i in 1:nrow(strat)) {
    year_use <- strat$year[i]
    
    temp.dens <- density_obs %>%
      filter(year(datetime) == year_use)
    
    if (nrow(temp.dens) >= 300) {
      #run length encoding according to the strat var
      temp.rle <- rle(temp.dens$strat)
      
      #what is the max length for which the value is "norm"
      strat$length[i] <- max(temp.rle$lengths[temp.rle$values==1], 
                             na.rm = T)
      
      #stratification dates
      rle.strat <- data.frame(strat = temp.rle$values, 
                              lengths = temp.rle$lengths)
      
      # Get the end of ech run
      rle.strat$end <- cumsum(rle.strat$lengths)
      # Get the start of each run
      rle.strat$start <- rle.strat$end - rle.strat$lengths + 1
      
      # Sort rows by whehter it is stratified or not
      rle.strat <- rle.strat[order(rle.strat$strat), ]
      
      start.row <- rle.strat$start[which(rle.strat$length == max(rle.strat$lengths)
                                         & rle.strat$strat == 1)] 
      #gets the row with the start date
      #of the run which has the max length and is 1
      
      end.row <- rle.strat$end[which(rle.strat$length == max(rle.strat$lengths)
                                     & rle.strat$strat == 1)] 
      #gets the row with the end date
      #of the run which has the max length and is TRuE
      
      strat$start[which(strat$year == year_use)] <- as.character(temp.dens$datetime[start.row])
      strat$end[which(strat$year == year_use)] <- as.character(temp.dens$datetime[end.row])
      
    }
    

  }
  strat <- strat |> 
    mutate(start = as_date(start),
           end = as_date(end))
  return(strat)
}
