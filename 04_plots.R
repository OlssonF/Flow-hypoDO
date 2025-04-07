library(ggpubr)
library(tidyverse)
source("R/plotting_functions.R")
possible_predictors <- c("inflow_Q", "inflow_T", 
                         "inflow_DOconc_mgL", 
                         "density_difference","schmidt.stability" ,  
                         "inflowD", "Kz_4.5",
                         "P10", "Qtot")


# Figure 2 - DO observations
ggarrange(make_obs_plot(variable = 'DO_5', 
                        all_data = summer_anoxic_periods,
                        ylab = expression(atop(paste("Dissolved oxygen"), 
                                               paste("concentration ", (mg~L^-1))))),
          make_obs_plot(variable = 'inflow_Q', 
                        all_data = summer_anoxic_periods,
                        ylab = expression(atop(paste("Discharge"), 
                                               paste((m^3~s^-1))))),
          make_obs_plot(variable = 'schmidt.stability', 
                        all_data = summer_anoxic_periods,
                        ylab = expression(atop(paste("Schmidt stability"),
                                               paste((W*m^-2))))),
          make_obs_plot(variable = 'inflowD',
                        all_data = summer_anoxic_periods,
                        ylab = expression("Inflow depth"~~(m))),
          make_obs_plot(variable = 'density_difference',
                        all_data = summer_anoxic_periods,
                        ylab = expression(atop(paste("Inflow - 5 m density"),
                                               paste("difference"~~(kg~m^-3))))),
          make_obs_plot(variable = 'Kz_4.5',
                        all_data = summer_anoxic_periods,
                        ylab = expression(atop(paste("log"[10]*"vertical diffusivity"~~(K[z])),
                                               paste("(",m^2*s^-1,")")))),
          make_obs_plot(variable = 'P10',
                        all_data = summer_anoxic_periods,
                        ylab = expression("Wind power"~~(mW~m^-2))),
          make_obs_plot(variable = 'Qtot',
                        all_data = summer_anoxic_periods,
                        ylab = expression("Advective heat flux"~~(W~m^-2))),
          
          ncol = 1, align = "hv")
