#merging all the raw data

library("plyr")
library("dplyr")
library("purrr")
library("tidyr")
library("lubridate")
library("ggplot2")
library("zoo")

file_list <- list.files(path = "data/raw_data/", pattern ="^P.*\\.csv$") #starts with P and ends in .csv

# Generates a single merged file with all time periods (T) and locations (P)
for (file in file_list){
  
  # if the merged o2_data doesn't exist, create it
  if (!exists("o2_data")){
    o2_data <- read.table(paste0("data/raw_data/",file), header=TRUE, sep=",")
  }
  
  # if the merged o2_data does exist, append to it
  if (exists("o2_data")){
    temp_dataset <-read.table(paste0("data/raw_data/",file), header=TRUE, sep=",")
    o2_data<-rbind(o2_data, temp_dataset)
    rm(temp_dataset)
  }
  
}

# Make sure the columns are the right format
o2_data$DateTime <- ymd_hms(o2_data$DateTime)
o2_data$P <- as.factor(o2_data$P)

glimpse(o2_data)

#remove duplicate values
o2_data <- distinct(o2_data,
                    .keep_all = T)
#=================================#



# ======= QA/QC =========
# remove rows where Q < 0.7, manufacturer quality value
o2_data_QA <- o2_data[which(o2_data$Q >= 0.7),]


##### check 1 Gross range tests on DO and Temp
# DO ranges for each sensor and each season
DO_max <- 16 # gross max
DO_min <- 0 # gross min

# more specific limits
DO_winter_min <- 5 # winter not going to drop below 5mgL fully mixed
surface_DO_min <- 0 # surface shouldn't drop below 0


# temperature ranges for each sensor and each season
temp_max <- 27
temp_min <- 0

# more specific limits
bottom_temp_max <- 15 # bottom water won't exceed 15 

# GR (gross range) flag will be true if either the DO concentration or the temperatures are 
# outside the ranges specified about
# after the first check the flag is not over-written unless T
o2_data_QA_GR <- o2_data_QA %>%
  mutate(Season = factor(format(as.yearqtr(as.yearmon(DateTime) + 1 / 12), "%q"),
                         levels = 1:4,
                         labels = c("Winter", "Spring", "Summer", "Autumn")),
         GR_flag = ifelse(DO_mgl >= DO_min & # more than min
                            DO_mgl <= DO_max,  # less than max
                          F, 
                          T), # if not within range --> NA
         GR_flag = ifelse(Season == "Winter" & DO_mgl < DO_winter_min, 
                          # DO can't be less than 5 mgL during winter
                          T, 
                          GR_flag),
         GR_flag = ifelse(P == 1 & DO_mgl < surface_DO_min, 
                          #surface DO should not drop below 5 mgL
                          T,
                          GR_flag),
         GR_flag = ifelse(Temp_C >= temp_min & Temp_C <= temp_max,
                          GR_flag, 
                          T),
         GR_flag = ifelse(P == 4 & Temp_C > bottom_temp_max, # max temp for bottom sensor
                          T, 
                          GR_flag))

#### checks 2    
# spike test (abs change betweeen values)
# flat line test (values remain the same )
# rate of change test (does the value change rapidly)
o2_data_QA_allflags <- o2_data_QA_GR %>%
  arrange(P, DateTime) %>%
  group_by(P) %>% # do the checks per sensor
  mutate(spike_ref =  (lag(DO_mgl) + lead(DO_mgl))/2,# average of the neighbouring values
         DO_spike = abs(DO_mgl - spike_ref),
         
         abs_change = abs(DO_mgl - lag(DO_mgl)),
         sd25 = rollapply(DO_mgl, width =25, FUN = sd, fill = NA, align = "r"), 
         spike_flag = ifelse(DO_spike > 5, 
                             T,
                             F),
         roc_flag = ifelse(abs_change > 4*sd25 & # rate of change not greater than 4 SD of previous 25 values
                             floor(DO_mgl) != 0, 
                           T, F),
         fl_flag = ifelse(floor(DO_mgl) != 0 &
                            DO_mgl == lag(DO_mgl, 1) & 
                            DO_mgl == lag(DO_mgl, 2), # DO value equal to previous 2
                          T, F), 
         abschange_flag = ifelse(abs_change > 4, 
                                 T, F)) #abs change of  > 8 is suspected erro

# if any of the data has been flagged set to NA
test_QC <- o2_data_QA_allflags %>%
  mutate(DO_auto = ifelse(GR_flag == T | # is there a flag? make NA
                            roc_flag == T |
                            fl_flag == T| 
                            spike_flag == T|
                            abschange_flag == T,
                          NA, DO_mgl),
         # if there is a spike. absolute or rate of change flag
         #  set the next value to NA
         DO_test = ifelse(lag(abschange_flag) == T | lag(spike_flag) == T | lag(roc_flag) == T,
                          NA, DO_auto),
         # if the previous value NA and there is a spike. absolute or rate of change flag previously
         #  set the next value to NA
         DO_test = ifelse(is.na(lag(DO_test)) & 
                            lag(abschange_flag, 2) == T | lag(spike_flag, 2) == T | lag(roc_flag, 2) == T,
                          NA, DO_test),
         # if the next value NA and there is a spike. absolute or rate of change flag previously
         #  set the previous value to NA
         DO_test = ifelse(lead(abschange_flag) == T | lead(spike_flag) == T | lead(roc_flag) == T,
                          NA, DO_test))

ggplot(test_QC, aes(x=DateTime, y=DO_test, colour =P)) +
  geom_line() +
  #geom_point(aes(colour = flag)) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%d %b") +
  coord_cartesian(xlim = c(ymd_hms("2018-09-01 00:00:00"),
                           ymd_hms("2018-12-01 00:00:00"))) +
  facet_wrap(~P)


# seems to still be some errors, 
# check these manually by plotting month by month
# P1
# 5/7/19 22:00 - 6/7/19 03:00
# P2
# 22/03/19 07:00 - 09:15
# 11/7/19 03:00 - 04:30
# 11/7/19 18:45 - 20:00
test_QC1 <- test_QC %>%
  mutate(DO_man = ifelse(P == 1 & 
                           between(DateTime, 
                                   ymd_hms("2019-07-05 22:00:00"), 
                                   ymd_hms("2019-07-06 03:00:00")),
                         NA, 
                         DO_test), 
         DO_man = ifelse(P == 2 & 
                           between(DateTime, 
                                   ymd_hms("2019-03-22 07:00:00"), 
                                   ymd_hms("2019-03-22 09:15:00")) |
                           between(DateTime, 
                                   ymd_hms("2019-07-11 03:00:00"), 
                                   ymd_hms("2019-07-11 04:30:00")) |
                           between(DateTime, 
                                   ymd_hms("2019-07-11 18:45:00"), 
                                   ymd_hms("2019-07-11 20:00:00")),
                         NA, 
                         DO_man))

test_QCed <- test_QC1[,c("DateTime", "P", "DO_man", "Temp_C")] 

write.csv(test_QCed, file = "./data/DO_QC.csv", row.names = F)


# how many missing values at 3 and 5 m
missing_3_5 <- test_QC1 %>%
  select(DateTime, P, DO_man) %>%
  filter(P == 2 | P==3 ) # sensor 2 at 3m, sensor 3 at 5m

length(which(is.na(missing_3_5$DO_man)))/ length(missing_3_5$DO_man) 

