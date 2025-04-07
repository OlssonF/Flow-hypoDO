# Data manipulation ------------------------
# Read in the raw data
# Do data processing e.g. data aggregation
# Generate a large data frame for analysis

# ------------------------------------------#

# Packages
library(tidyverse)
library(rLakeAnalyzer)
library(zoo)
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

all_dates <- bind_rows(data.frame(datetime = seq.Date(anoxic_dates$start[1],
                                                      anoxic_dates$end[1],
                                                      'day')),
                       data.frame(datetime = seq.Date(anoxic_dates$start[2],
                                                      anoxic_dates$end[2],
                                                      'day')))


all_dates_hourly <- bind_rows(data.frame(datetime = seq(as_datetime(anoxic_dates$start[1]),
                                                        as_datetime(anoxic_dates$end[1]),
                                                        'hour')),
                              data.frame(datetime = seq(as_datetime(anoxic_dates$start[2]),
                                                        as_datetime(anoxic_dates$end[2]),
                                                        'hour')))

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

# Inflow data (EIDC) -------------------
temp <- file.path('data', 'raw_data', 'flowdata.zip')
exdir <- file.path('data', 'raw_data', 'flowdata') 
dir.create(exdir, showWarnings = F)

if (!file.exists(exdir)) {
  download.file("https://data-package.ceh.ac.uk/data/2883aaf1-6148-49cb-904a-d271a028c716.zip", temp)
  unzip(zipfile = temp, exdir = exdir)  
}


inflow_dat <- read_csv(file.path(exdir, list.files(exdir, recursive = T, 
                                                   pattern = 'elter_inflow_2012-2019.csv')))
unlink(temp)

# calculate daily inflow data
inflow_dat_daily <- 
  inflow_dat |> 
  mutate(datetime = as_date(DateTime)) |> 
  reframe(.by = datetime,
          across(inflow_Q:inflow_T, ~ mean(.x, na.rm = T)))

## Inflow DO estimates -----------------
inflow_dat_daily_LM <- 
  inflow_dat_daily |> 
  select(datetime, inflow_T) |> 
  rename(wtr = inflow_T)

# using lake metabolizer to estimate the 100% saturation DO concentration based on water temperature
inflow_dat_daily <- inflow_dat_daily |> 
  mutate(inflow_DOconc_mgL = LakeMetabolizer::o2.at.sat(inflow_dat_daily_LM)$do.sat, # assume a 100% sat mg/L
         # use the concentration to estimate the flux (mg s-1)
         # convert mg L-1 to mg m-3 ( * 1000)
         # multiply by discharge m3 s-1 ==> flux (mg s-1)
         inflow_DOflux_mgs = (inflow_DOconc_mgL * 1000) * inflow_Q, 
         
         # convert the flux to kg day-1
         # multiply by s in a day (86400) and divide by 1e6 (mg --> kg)
         inflow_DOmass_kg = (inflow_DOflux_mgs * 86400) / 1e6) 



# Thermistor data (EIDC) -------------
temp2 <- file.path('data', 'raw_data', 'insitudata.zip')
exdir2 <- file.path('data', 'raw_data', 'insitudata') 
dir.create(exdir2, showWarnings = F)

if (!file.exists(exdir2)) {
  download.file("https://data-package.ceh.ac.uk/data/37f17f6c-66f6-454c-bd52-7c601ef20ca2.zip", temp2)
  unzip(zipfile = temp2, exdir = exdir2)  
}

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

# write hourly data for lhfa
thermistor_dat |> 
  right_join(all_dates_hourly, by = join_by(DateTime == datetime)) |>  
  mutate(dateTime = format(DateTime, "%Y-%m-%d %H:%M")) |> 
  rename(temp1 = T_1m) |> 
  select(any_of(c('dateTime', 'temp1'))) |> 
  mutate(temp1 = na.approx(temp1)) |> 
  write_delim(file = 'HeatFluxAnalyzer/data/Elterwater.wtr',
              delim = '\t',
              quote = "none")
## Derived physical metrics ---------------
### Schmidt stability ------------------------------------
schmidt_stability_daily <- ts.schmidt.stability(thermistor_dat_daily, bathy = bathy) 

### Kz (vertical diffusivity)-----------------------------
source('R/kz_calculation.R')
kz_daily <- ts.kz(wtr = thermistor_dat_daily,
                  bathy = bathy_int) |> # requires the bathymetry at the interpolate depths
  select(datetime, Kz_4.5)

### Nominal intrusion depth --------------------------------
source('R/inflowD_calculation.R')
intrusion_depth <- thermistor_dat_daily |> 
  inner_join(inflow_dat_daily, by = 'datetime') |> 
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
  mutate(density_difference = water.density(wtr_6) - water.density(wtr_1)) |> 
  select(datetime, density_difference)
#------------------------------------------#

### Stratification distance ---------------
# Time since/before onset/overturn
source('R/strat_dates.R')
strat_dates <- calc_strat_dates(wtr = thermistor_dat_daily, 
                                use_depths = c(1,6), 
                                density_diff = 0.1)
#-----------------------------------------#


# Meteorological (EIDC) -------------------
library(suntools)
temp3 <- file.path('data', 'raw_data', 'metdata1.zip')
temp4 <- file.path('data', 'raw_data', 'metdata2.zip')
temp5 <- file.path('data', 'raw_data', 'metdata3.zip')
exdir3 <- file.path('data', 'raw_data', 'metdata1') 
exdir4 <- file.path('data', 'raw_data', 'metdata2') 
exdir5 <- file.path('data', 'raw_data', 'metdata3') 
dir.create(exdir3, showWarnings = F)
dir.create(exdir4, showWarnings = F)
dir.create(exdir5, showWarnings = F)

if (!file.exists(exdir3)) {
  download.file("https://data-package.ceh.ac.uk/data/603629d9-618c-4a26-8a1d-235a4c8f4791.zip", temp3)
  unzip(zipfile = temp3, exdir = exdir3)  
}

if (!file.exists(exdir)) {
  download.file("https://data-package.ceh.ac.uk/data/467942bd-2cd5-4038-bc4f-5d77535d99f1.zip", temp4)
  unzip(zipfile = temp4, exdir = exdir4)  
}

if (!file.exists(exdir5)) {
  download.file("https://data-package.ceh.ac.uk/data/3df05e85-2c56-4bd9-9918-44b760e20b2e.zip", temp5)
  unzip(zipfile = temp5, exdir = exdir5)  
}

met_dat_2018 <- read_csv(file.path(exdir3, list.files(exdir3, recursive = T, 
                                                      pattern = 'blel-2016_2017_2018.csv')),
                         skip = 2, 
                         col_names = c('datetime','wtr_0.5', 'wtr_1','wtr_2','wtr_3','wtr_4','wtr_5','wtr_6',
                                       'wtr_7','wtr_8','wtr_9','wtr_10','wtr_12', 'airT', 'sw', 'wnd')) |> 
  mutate(datetime = dmy_hm(datetime)) |> 
  select(datetime, airT, sw, wnd) |> 
  filter(year(datetime) == 2018)

met_dat_2019 <- read_csv(file.path(exdir4, list.files(exdir4, recursive = T, 
                                                      pattern = 'blel-2019.csv')),
                         skip = 2, 
                         col_names = c('datetime','wtr_0.5', 'wtr_1','wtr_2','wtr_3','wtr_4','wtr_5','wtr_6',
                                       'wtr_7','wtr_8','wtr_9','wtr_10','wtr_12', 'airT', 'sw', 'wnd'))  |> 
  mutate(datetime = dmy_hm(datetime)) |> 
  select(datetime, airT, sw, wnd) |> 
  filter(year(datetime) == 2019)


rh_2018_2019 <- read_csv(file.path(exdir5, list.files(exdir5, recursive = T, 
                                                      pattern = 'blel-rh_2012_2019.csv')),
                         skip = 1,
                         col_names = c('date', 'hour', 'RH', 'N')) |> 
  mutate(date = dmy(date)) |> 
  filter(year(date) %in% c(2018, 2019)) |> 
  mutate(datetime = ymd_h(paste0(date, hour))) |>
  select(datetime, RH)

# combine the dataframes from 2018, 2019, and the rh dataset
met_dat <- bind_rows(met_dat_2018, met_dat_2019) |> 
  full_join(rh_2018_2019) |> 
  # simple QAQC
  mutate(sw = ifelse(sw < 0, 0, sw),
         wnd = ifelse(wnd < 0, 0, wnd))

unlink(temp3)
unlink(temp4)
unlink(temp5)


# Do some QA/QC
long_lat <- matrix(c( -3.0350, 54.4287), nrow = 1)
check_dates <- as.POSIXct(unique(format(met_dat$datetime, "%Y-%m-%d")), tz = 'GMT')
sunrise <- suntools::sunriset(crds = long_lat,
                              dateTime = check_dates,
                              direction = 'sunrise',
                              POSIXct.out = T)[2] |>
  rename(sunrise = time) |> 
  mutate(date = as_date(format(sunrise, "%Y-%m-%d")))

sunset <- suntools::sunriset(crds = long_lat,
                             dateTime = check_dates,
                             direction = 'sunset',
                             POSIXct.out = T)[2] |> 
  rename(sunset = time) |> 
  mutate(date = as_date(format(sunset, "%Y-%m-%d")))


# check that the sw is 0 between sunrise and sunset
met_dat <- met_dat |> 
  mutate(date = as_date(datetime)) |> 
  full_join(sunrise, by = join_by(date), relationship = 'many-to-many') |> 
  full_join(sunset, by = join_by(date), relationship = 'many-to-many') |> 
  # set sw to zero if it's nighttime
  mutate(sw = ifelse(between(datetime, sunrise, sunset), sw, 0)) |> 
  select(datetime, airT, sw, wnd, RH)


# Export the met data for use in MATLAB Lake Heat Flux Analyzer
vars <- c('airT', 'sw', 'wnd', 'RH')

for (var in vars) {
  met_dat |> 
    right_join(all_dates_hourly) |> 
    mutate(dateTime = format(datetime, "%Y-%m-%d %H:%M")) |> 
    select(any_of(c('dateTime', var))) |> 
    mutate(across(all_of(var), na.approx, .names = var)) |> 
    
    head(n = 20) |> 
    write_delim(file = paste0('HeatFluxAnalyzer/data/Elterwater.',var),
                delim = '\t',
                quote = "none")
}

# THIS NEXT STEP NEEDS TO BE RUN IN MATLAB! ---------------
#   THE SCRIPT IS AVAILABLE IN THE SUBDIRECTORY ("./HeatFluxAnalyzer/run_lhfa.m")
#   The directory contains the configuration files (.hfx) and the matlab code to run LHFA

# This generates a .txt file with the output that we now read in. 
lhfa_daily <- read_delim('HeatFluxAnalyzer/data/Elterwater_results.txt', delim = '\t') 

# simplify column names
colnames(lhfa_daily) <- str_split_i(colnames(lhfa_daily), pattern = ' ', i = 1)

lhfa_daily <- lhfa_daily |> 
  select(DateTime,
         u10, # wind speed at 10 m
         Qtot, # total surface heat flux
         C_D10) |>  # transfer coefficient of momentum at 10 m
  mutate(P10 = C_D10 * 1.2 * (u10^3),
         # calculate windpower
         # C_D (transfer coefficient) * 1.2 (density of air) * cube of the wind
         datetime = ymd(format(DateTime, "%Y-%m-%d"))) |> 
  select(datetime, P10, Qtot)

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

DO_anoxic_periods |> 
  pivot_wider(names_from = depth, names_prefix = 'DO_',
              values_from = DO_mgl_daily,
              id_cols = datetime) |> 
  select(-DO_0.5, -DO_6) |> 
  left_join(date_metrics, by = join_by(datetime)) |> 
  left_join(inflow_dat_daily, by = join_by(datetime)) |> 
  left_join(density_difference_daily, by = join_by(datetime)) |> 
  left_join(intrusion_depth, by = join_by(datetime)) |> 
  left_join(kz_daily, by = join_by(datetime)) |> 
  left_join(lhfa_daily, by = join_by(datetime)) |>
  left_join(schmidt_stability_daily, by = join_by(datetime)) |> 
  write_csv(file = 'output/all_data.csv')
