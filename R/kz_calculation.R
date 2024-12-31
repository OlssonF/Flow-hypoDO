# Kz calculation - gradient flux method ----------
# this method of estimating Kz uses water temperature profiles and
# applies the gradient-flux method

# Kz = sum(ak * Vk)/ Gi * Ai

# top of the equation is the temporal increase in temperature * volume in the box
# for each box below the specified depth e.g. for Kz at 4m it is the sum 
# of the ak* volume for 4-5 and 5-6m... etc.

# ak = dT/dt in each box below k=1
# V = volume of box i

# Gi = vertical temperature gradient at the upper boundary of box i
# Gi = changeT(i-1) / changeD
# Ai is area at upper boundary of box i

#' @param wtr data frame of water temperatures at measured depths (1,2,3m), requires temperatures at the 'middle' of each box
#' data are formatted like rLakeAnalyzer (datetime, wtr_*)
#' @param bathy data frame of depths and areas/volumes, requires bathymetry at the top/bottom of boxes (0.5, 1.5, 2.5 etc.)
#' interpolated from the 'observations'
#' @return dataframe of kz values at different depths
#' 
ts.kz <- function(wtr, bathy) {
  # message('calculating Kz via gradient heat flux method')
  # Top of the equation -------------
  # calculate dT/dt for the water temperatures (numeric)
  dT_dt <- wtr %>%
    mutate_at(vars(matches('wtr')), dT.dt)
  
  # then multiply the dT_dt by the volume of the box (Vi)
  # dT/dt * Vk
  dT_dt_Vk <- dT_dt %>%
    mutate(across(matches('wtr'), ~dT.dt.V(column = .x, 
                                           depth = as.numeric(gsub('wtr_', '', cur_column())),
                                           bathy = bathy)))
  
  
  # for each depth you add up that layer plus the ones below
  sum_dT_dt_Vk <- 
    dT_dt_Vk %>%
    select(datetime, wtr_4, wtr_5, wtr_6) %>%
    # writing over the columns but should be okay because don't then use
    # that column in the next mutate
    mutate(wtr_4 = wtr_4 + wtr_5 + wtr_6, # box = 4 m
           wtr_5 = wtr_5 + wtr_6,         # box = 5 m
           wtr_6 = wtr_6)                # box = 6 m
  
  # message('top of calc done')
  
  # bottom of the equation --------------------
  # temperature gradient between the boxes
  # gradient between each layer divided by the height of the layer (1m)
  Gi <- data.frame(date = wtr$datetime,
                   
                   # wtr[,2:6] is wtr_1:wtr_5
                   # wtr[,3:7] is wtr_2:wtr_6
                   # so the calc is 1 - 2; 2 - 3; 3 - 4 etc..
                   (wtr[,2:(ncol(wtr) - 1)] - wtr[,3:ncol(wtr)])/1) # then divide by height of layer (ie 1m)
  
  # use the boundary values as column names
  # only 1.5 - 5.5m
  colnames(Gi) <- c("Date", paste0("wtr_", depths_boundaries[-c(1,7)])) 
  
  
  # areas at boundaries
  Gi_Ai <- Gi |> 
    as_tibble() |> 
    # mutate_if(is.numeric, Gi.Ai) |> 
    mutate(across(matches('wtr'), ~Gi.Ai(column = .x, 
                                         depth = as.numeric(gsub('wtr_', '', cur_column())),
                                         bathy = bathy))) |> 
    # only valid for the deeper layers
    select(Date, wtr_3.5, wtr_4.5, wtr_5.5)
  
  # message('bottom of calc done')
  
  # final calculation of Kz using the gradient heatflux method -------
  # kz = sum_dT_dt_V / Gi_Ai
  Kz_ghf <- data.frame(Date = Gi_Ai$Date,
                       #calculate Kz for each layer [-date]
                       sum_dT_dt_Vk[,-1] / Gi_Ai[,-1])
  
  colnames(Kz_ghf) <- c("Date", paste0("Kz_", gsub("wtr_", "", colnames(Gi_Ai[-1]))))
  
  
  # Kz values cannot be less than the rate of molecular diffusion
  rate_diffusion <- 1.4e-7
  
  Kz_ghf <- Kz_ghf %>%
    # if value is less than molecular diffusion set to rate of diffusion
    mutate_if(is.numeric, 
              ~ifelse(. < rate_diffusion, # value is < rate of diffusion
                      rate_diffusion,
                      .)) 
  
  return(as_tibble(Kz_ghf))
  
}


# ==== top of the equation ====

# calculate ak = dT/dt

# change from previous value
# for daily values divide by number of s in a day = 86400
dT.dt <- function(var) {
  abs_change <-  var - lag(var)
  # divide by number of seconds in a day
  dT_dt <- abs_change/86400
  
  return(dT_dt)
}

# function to look-up the depth volume for each column and multiply by this value
dT.dt.V <- function(column, depth, bathy) {
  # subtract 0.5 to get the right layer, substitute() does something to get the colname
  # e.g.for the temperatures measured at 1m this requires the volume layer between 0.5 and 1.5 m
  V <- bathy$volume[which(bathy$depths == as.numeric(depth - 0.5))]
  dT_dt_V <- column * V
  return(dT_dt_V)
}

# multply the temperature gradient at the boundary by the area of that boundary
Gi.Ai <- function(column, depth, bathy) {
  # get the area for the right boundary
  A <- bathy$areas[which(bathy$depths == depth)]
  Gi_Ai <- column * A
  return(Gi_Ai)
}