# Fitting Random Forest
library(caret)
library(h2o)
library(tidyverse)
library(rsample)
# 1. Read in data ---------------

all_dat <- read_csv('output/all_data.csv') |> 
  na.omit() |> 
  select(-any_of(c('DO_3', 'year', 'Kz_3.5', 'Kz_5.5')))


corrplot::corrplot(cor(all_dat %>% select_if(is.numeric), 
                       method = "spearman"))

# Filter for correlations above a certain threshold (e.g., 0.5)
cor_matrix <- cor(all_dat |> select_if(is.numeric), 
                  method = "spearman") |> 
  as.data.frame() |> 
  rownames_to_column(var = 'var1') |>  
  gather(var2, value, -var1) |> 
  filter(abs(value) > 0.85,
         var1 != var2,
         var1 != 'DO_5',
         var2 != 'DO_5')

# remove highly correlated values
all_dat <- all_dat |> select(-any_of(c('datetime', 'doy', 'n_overturn', 'schmidt.stability', 'inflow_T', 'inflow_DOflux_mgs', 'inflow_DOmass_kg')))
#### 2. Nested CV assessment of Random Forest model ####

# this section adapted from https://bradleyboehmke.github.io/HOML/random-forest.html
# other potentially useful links:
# https://machinelearningmastery.com/nested-cross-validation-for-machine-learning-with-python/#comment-578556
# https://stats.stackexchange.com/questions/65128/nested-cross-validation-for-model-selection
# https://machinelearningmastery.com/nested-cross-validation-for-machine-learning-with-python/

# Train & cross-validate a RF model
# random forest is useful when one or two predictor variables may swamp out the signal
# by randomly subsetting the possible predictors at each node, other predictor variables
# can be more "seen" - wind direction and speed was coming out as overly important in
# regression tree analysis, so RF might help assess a wider range of predictors

# main hyperparameters when tuning random forests: 1) number of trees, 2) number of features 
# to consider at any given split (mtry), 3) complexity of each tree, 4) sampling scheme
# and 5) splitting rule to use during tree construction.

# to make reproducible
set.seed(1)

# split dataset
n_folds <- 7
folds <- createFolds(all_dat$DO_5, k = n_folds)

# Please download and install the latest version of h2o from http://h2o.ai/download/
# make sure h2o connection is fresh
h2o.shutdown(prompt = FALSE)
h2o.init() 
Sys.sleep(10)

for(i in 1:n_folds){
  h2o.shutdown(prompt = FALSE)
  Sys.sleep(10)
  h2o.init() 
  
  data_in <- all_dat[-folds[[i]],]
  data_h2o <- as.h2o(data_in)
  
  response <-"DO_5"
  predictors <- setdiff(names(data_h2o),response)
  
  # number of features
  n_features <- length(setdiff(names(data_h2o), "DO_5"))
  
  hyper_grid <- list(mtries = floor(n_features * c(.25, .4, .5, .6, .75)),
                     min_rows = c(1, 3, 5, 10),
                     max_depth = c(10, 20, 30),
                     sample_rate = c(1))
  
  # random grid search strategy
  search_criteria <- list(
    strategy = "RandomDiscrete",
    stopping_metric = "RMSE",
    stopping_tolerance = 0.001,   # stop if improvement is < 0.1%
    stopping_rounds = 10,         # over the last 10 models
    max_runtime_secs = 60*5      # or stop search after 5 min.
  )
  
  # perform grid search 
  random_grid <- h2o.grid(
    algorithm = "randomForest",
    grid_id = "rf_random_grid",
    x = predictors,
    nfolds=3,
    y = response, 
    training_frame = data_h2o,
    hyper_params = hyper_grid,
    ntrees = n_features * 10,
    seed = 123,
    search_criteria = search_criteria
  )
  
  # collect the results and sort by our model performance metric 
  # of choice
  random_grid_perf <- h2o.getGrid(grid_id = "rf_random_grid", 
                                  sort_by = "RMSE", 
                                  decreasing = FALSE)
  # random_grid_perf
  # save the top model, by mse
  best_1 <- h2o.getModel(random_grid_perf@model_ids[[1]])
  
  if(i==1){
    save_1 <- (best_1@model[["model_summary"]])
    save_1$min_rows<-as.numeric(random_grid_perf@summary_table[1,2])
    save_1$mtries<-as.numeric(random_grid_perf@summary_table[1,3])
    save_1$inner_rmse<-as.numeric(random_grid_perf@summary_table[1,6])
  }else{
    save_1[i,]<-(best_1@model[["model_summary"]])
    save_1$min_rows[i]<-as.numeric(random_grid_perf@summary_table[1,2])
    save_1$mtries[i]<-as.numeric(random_grid_perf@summary_table[1,3])
    save_1$inner_rmse[i]<-as.numeric(random_grid_perf@summary_table[1,6])
  }
  
  #then evaluate model on the test set
  data_test <- all_dat[folds[[i]],]
  test_h2o <- as.h2o(data_test)
  best_perf <- h2o.performance(model=best_1,newdata = test_h2o)
  rmse_in <- best_perf@metrics$RMSE
  if(i==1){
    save_1$outer_rmse <- as.numeric(rmse_in)
  }else{
    save_1$outer_rmse[i] <- as.numeric(rmse_in)
  }
  
}


#### 2a. Fit Final Random Forest model #### 
# With the above assessment of the model fiting routine, we now run the inner loop on the whole dataset 
# to establish final model parameters
h2o.shutdown(prompt = FALSE)
h2o.init()

data_in <- all_dat
data_h2o <- as.h2o(data_in)

response <- "DO_5"
predictors <- setdiff(names(data_h2o),response)

# number of features
n_features <- length(setdiff(names(data_h2o), response))

hyper_grid <- list(mtries = floor(n_features * c(.25, .4, .5, .6, .75)),
                   min_rows = c(1, 3, 5, 10),
                   max_depth = c(10, 20, 30),
                   sample_rate = c(0.6,0.75,0.9,1))

# random grid search strategy
search_criteria <- list(strategy = "RandomDiscrete",
                        stopping_metric = "RMSE",
                        stopping_tolerance = 0.001,   # stop if improvement is < 0.1%
                        stopping_rounds = 10,         # over the last 10 models
                        max_runtime_secs = 60*5)      # or stop search after 5 min.

# perform grid search 
random_grid <- h2o.grid(algorithm = "randomForest",
                        grid_id = "rf_random_grid",
                        x = predictors,
                        nfolds = 3,
                        y = response, 
                        training_frame = data_h2o,
                        hyper_params = hyper_grid,
                        ntrees = n_features * 10,
                        seed = 123,
                        search_criteria = search_criteria)

# collect the results and sort by our model performance metric 
# of choice
random_grid_perf <- h2o.getGrid(grid_id = "rf_random_grid", 
                                sort_by = "RMSE", 
                                decreasing = FALSE)
random_grid_perf

# save the top model, by rmse
best_1 <- h2o.getModel(random_grid_perf@model_ids[[1]])


#Now, create the final model using parameters from above grid search
h2o.init()

#split into training and test data
set.seed(2) 
data_split <- initial_split(all_dat, strata = "DO_5")
DO_train <- training(data_split)
DO_test <- testing(data_split)
train_h2o <- as.h2o(DO_train)
test_h2o <- as.h2o(DO_test)

#assign response variable
Y <- "DO_5" 

#assign predictor variables
X <- setdiff(names(train_h2o), Y)

# use parameters from top fitting model:
best_rf <- h2o.randomForest(x = X, y = Y, 
                            training_frame = train_h2o,
                            ntrees = best_1@params$actual$ntrees, 
                            mtries = best_1@params$actual$mtries,
                            validation_frame = test_h2o,
                            sample_rate = 1,
                            max_depth = best_1@params$actual$max_depth, 
                            min_rows = best_1@params$actual$min_rows,
                            nfolds = best_1@params$actual$nfolds,
                            fold_assignment = "Modulo",
                            keep_cross_validation_predictions = best_1@params$actual$keep_cross_validation_predictions,
                            seed = 123, 
                            stopping_rounds = 50, 
                            stopping_metric = "RMSE", score_each_iteration = T,
                            stopping_tolerance = 0)

best_rf

#pull final model RMSE w/hold-out test data
h2o.performance(best_rf, newdata = test_h2o)@metrics$RMSE

h2o.varimp(best_rf)

h2o.varimp_plot(best_rf)

h2o.pd_plot(best_rf, test_h2o, column = 'inflow_Q', show_rug = T)

# requires plot3Drgl - generate a heatmap
library(plot3Drgl)
test <- h2o.partialPlot(best_rf, test_h2o, col_pairs_2dpdp = list(c("inflow_Q", "Kz_4.5")), 
                        plot = FALSE)

ggplot(test, aes(x = inflow_Q, y=Kz_4.5, fill = mean_response)) +
  geom_tile()+
  scale_fill_viridis_c() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))
