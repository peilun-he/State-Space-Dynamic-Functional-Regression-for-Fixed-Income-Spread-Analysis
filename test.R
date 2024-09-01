library(reticulate) # import Python packages
library(ggplot2)
library(dplyr)
library(expm)
library(tidyr)
library(YieldCurve)
library(zoo) 
library(pracma)
library(plotly)
library(ftsa)
library(glmnet)
library(lmridge)
library(lmtest)
library(statespacer)
library(scales)

source("Functions/NS_loading.R")
source("Functions/match_maturity.R")
source("Functions/KF_NS.R")
source("Functions/kpca.R")

US <- read.csv("Data/US_yields.csv")
UK <- read.csv("Data/UK_yields.csv")

US <- US[, -1]
UK <- UK[, -1]

use_python("/Users/HPL/anaconda3/bin/python3") # select python version
pd <- import("pandas")
np <- import("numpy")
sk_dec <- import("sklearn.decomposition")
sk_met <- import("sklearn.metrics")
sk_ms <- import("sklearn.model_selection")

# Estimate hyper-parameters
score <- function(estimator, X, Y = NULL) {
  X_reduced <- estimator$fit_transform(X)
  X_preimage <- estimator$inverse_transform(X_reduced)
  return(-sk_met$mean_squared_error(X, X_preimage))
}

param_grid <- list(gamma = seq(from = 0.001, to = 1, by = 0.001))

Q <- 3 # number of factors. (max = 12)
kpca_machine <- sk_dec$KernelPCA(kernel = "rbf",  
                                 n_components = as.integer(Q), 
                                 fit_inverse_transform = TRUE)

set.seed(1234)
grid_search <- sk_ms$GridSearchCV(kpca_machine, 
                                  param_grid, 
                                  cv = as.integer(Q), scoring = score)
hp_fit <- grid_search$fit(t(as.matrix(US))) # hyper-parameter 

hp_fit$best_params_

# Extract factors
Q <- 3
U <- kpca(data = US, kernel = "rbf", gamma = hp_fit$best_params_$gamma, Q = Q)$U

write.csv(U, "Data/US_factors.csv", row.names = FALSE)
