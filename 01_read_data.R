# Data manipulation ------------------------
# Read in the raw data
# Do data processing e.g. data aggregation
# Generate a large data frame for analysis

# ------------------------------------------#

# Packages
library(tidyverse)
library(rLakeAnalyzer)
# DO data -----------------
# read in the QC-ed data, raw data that has uneven timesteps, roughly every 15 minutes
DO_dat <- read_csv("data/DO_QC.csv") |> 
  mutate(DateTime = ymd_hms(DateTime)) 

# calculate the hourly average DO at each depth
hourly_DO <- DO_dat  |> 
  mutate(Hour = floor_date(DateTime, unit  = "hour")) |> # round the hour to get an average
  group_by(P, Hour) |>
  summarise(DO_mgl = mean(DO_man, na.rm = T),  # mean DO
            Temp_C = mean(Temp_C, na.rm = T), # mean temperature
            n_DO = sum(!(is.na(DO_man))), # count the non-NA values
            .groups = 'drop') |> 
  mutate(DO_mgl = ifelse(n_DO < 3, # if there are less than 3 values used to calculate the mean
                         NA, DO_mgl)) #|> # set this to NA

# get all time, depth combinations
P_hour <- expand.grid(Hour=seq(min(hourly_DO$Hour),
                               max(hourly_DO$Hour),
                               "hour"),
                      P=1:4)
# complete the df
hourly_DO <- full_join(P_hour, hourly_DO, by = join_by(Hour, P))

# interpolate small gaps
hourly_DO_int <- hourly_DO |>
  group_by(P) |> # for each sensor
  mutate(DO_mgl = na.approx(DO_mgl, maxgap = 3, # interpolate gaps < 4hrs
                            na.rm = F)) 

# Calculate daily average DO
daily_DO <- hourly_DO_int |> # use the interpolated hourly data
  mutate(datetime = ymd(format(Hour, "%y-%m-%d"))) |> # round the hour to get an average
  group_by(P, datetime) |>
  summarise(DO_mgl_daily = mean(DO_mgl, na.rm = T),  # mean DO
            Temp_C_daily = mean(Temp_C, na.rm = T), # mean temperature
            n_DO = sum(!(is.na(DO_mgl))) # count the non-NA values
            , .groups = 'drop') |> 
  mutate(DO_mgl_daily = ifelse(n_DO < 18, # 2/3 of timesteps required
                               # if there are less than 18 values (= 18 hours) used to calculate the mean
                               NA, DO_mgl_daily)) # set this to NA

#get all time, depth combinations
P_day <- expand.grid(datetime=seq(min(daily_DO$datetime),
                              max(daily_DO$datetime),
                              "day"), P=1:4)

daily_DO <- full_join(P_day, daily_DO, 
                      by = join_by(datetime, P)) |> 
  mutate(depth = factor(P, levels = 1:4, labels = c(0.5,3,5,6))) |> 
  select(datetime, depth, DO_mgl_daily, Temp_C_daily)

# Generate a wide formats needed for some calculations
wide_daily_DO <- daily_DO %>%
  pivot_wider(names_from = depth,
              names_prefix = "wtr_", # wtr prefix needed in the lakeanalyzer functions
              values_from = "DO_mgl_daily",
              id_cols = "datetime") 

wide_daily_T <- daily_DO %>%
  pivot_wider(names_from = depth,
              names_prefix = "wtr_", # wtr prefix needed in the lakeanalyzer functions
              values_from = "Temp_C_daily",
              id_cols = "datetime") 

## Identify the period used for analysis ------------
# truncate the timeseries so we just have the "summer" period
# for DO this is defined as when the DO_5 is no longer continuously > 1mg/L.
source('R/anoxic_period.R')
anoxic_dates <- calc_anoxic_period(DO = daily_DO)

DO_anoxic_periods <- daily_DO |> 
  filter(between(datetime, anoxic_dates$start[1], anoxic_dates$end[1])|
           between(datetime, anoxic_dates$start[2], anoxic_dates$end[2]))

# Hypsometry/bathymetry ------------------
bathy <- read.delim("data/elter_bathymetry.dat")

depths <-  data.frame(depths = seq(0.5, 6.5, 0.5))
depths_boundaries <- seq(0.5, 6.5, 1)
depths_measurements <- seq(1, 6, 1)

bathy_int <- bathy |> 
  full_join(depths, by = "depths") |> 
  arrange(depths) |> 
  #interpolate the areas
  mutate(areas = ifelse(depths == 6.5, 0, 
                        na.approx(areas, na.rm = F))) |> 
  # then extract the boxes needed and calculate the volume based on these areas
  filter(depths %in% depths_boundaries) |> 
  # calculate the volume of each of these trapeziums
  # (a+b)/2 * h
  mutate(volume = ((lead(areas) + areas)/2) * (lead(depths) - depths))


# ----------------------------------------#

# EIDC flow data -------------------
temp <- file.path('data', 'raw_data', 'flowdata.zip')
exdir <- file.path('data', 'raw_data', 'flowdata') 
dir.create(exdir, showWarnings = F)

download.file("https://data-package.ceh.ac.uk/data/2883aaf1-6148-49cb-904a-d271a028c716.zip", temp)
unzip(zipfile = temp, exdir = exdir)

flow_dat <- read_csv(file.path(exdir, list.files(exdir, recursive = T, 
                                                 pattern = 'elter_inflow_2012-2019.csv')))
unlink(temp)

# calculate daily inflow data
flow_dat_daily <- 
  flow_dat |> 
  mutate(datetime = as_date(DateTime)) |> 
  reframe(.by = datetime,
          across(inflow_Q:inflow_T, ~ mean(.x, na.rm = T)))

# EIDC phys/chem/bio data -------------
temp2 <- file.path('data', 'raw_data', 'insitudata.zip')
exdir2 <- file.path('data', 'raw_data', 'insitudata') 
dir.create(exdir2, showWarnings = F)

download.file("https://data-package.ceh.ac.uk/data/37f17f6c-66f6-454c-bd52-7c601ef20ca2.zip", temp2)
unzip(zipfile = temp2, exdir = exdir2)

thermistor_dat <- read_csv(file.path(exdir2, list.files(exdir2, recursive = T, 
                                                        pattern = 'elter_t_profiles_2018-2019.csv')))

unlink(temp2)

thermistor_cols <- c(wtr_1 = 'T_1m',
                     wtr_2 = 'T_2m',
                     wtr_3 = 'T_3m',
                     wtr_4 = 'T_4m',
                     wtr_5 = 'T_5m',
                     wtr_6 = 'T_6m')

# calculate daily average profiles
thermistor_dat_daily <- 
  thermistor_dat |> 
  mutate(datetime = as_date(DateTime)) |> 
  reframe(.by = datetime,
          across(T_1m:T_6m, ~ mean(.x, na.rm = T))) |> 
  rename(any_of(thermistor_cols)) # rename columns for rlakeanalyzer


## Thermistor data and derived physical metrics ----------
### Schmidt stability ------------------------------------
schmidt_stability_daily <- ts.schmidt.stability(thermistor_dat_daily, bathy = bathy) 

### Kz (vertical diffusivity)-----------------------------
source('R/kz_calculation.R')
kz_daily <- ts.kz(wtr = thermistor_dat_daily,
                  bathy = bathy_int) # requires the bathymetry at the interpolate depths

### Nominal intrusion depth --------------------------------
source('R/inflowD_calculation.R')
intrusion_depth <- thermistor_dat_daily |> 
  inner_join(flow_dat_daily, by = 'datetime') |> 
  select(-inflow_Q) |> 
  rename(wtr_inflow = inflow_T) |>
  pivot_longer(wtr_1:wtr_6, names_to = 'depth', 
               names_prefix = 'wtr_', names_transform = list(depth = as.numeric),
               values_to = 'temp') |> 
  reframe(.by = datetime,
          inflowD = intrusion.depth(wtr = temp,
                                    depths = depth, 
                                    inflowT = wtr_inflow))
# ---------------------------------------#

### Density differences--------------------
density_difference_daily <- thermistor_dat_daily |> 
  mutate(density_difference = water.density(wtr_1) - water.density(wtr_6)) |> 
  select(datetime, density_difference)
#------------------------------------------#

### Stratification distance ---------------
# Time since/before onset/overturn
source('R/strat_dates.R')
strat_dates <- calc_strat_dates(wtr = thermistor_dat_daily, 
                                use_depths = c(1,6), 
                                density_diff = 0.1)
#-----------------------------------------#


# Date metrics ----------------------------- 
# date_metrics <-
date_metrics <- daily_DO |> 
  distinct(datetime) |> 
  mutate(doy = yday(datetime),
         year = year(datetime)) |> 
  left_join(strat_dates, by = join_by(year)) |> 
  mutate(n_onset = as.numeric(datetime - start), # days since onset
         n_overturn = as.numeric(end - datetime)) |> # days until overturn
  select(any_of(c('datetime', 'doy', 'year', 'n_onset', 'n_overturn')))

# Combine data -----------------------------
# Combine all dataframes to be used in analysis
