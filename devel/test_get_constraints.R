rm(list = ls())
source("spotfwdratedef.R")
source("estim_nss_zeroyields.R")

myconstraints <- get_constraints("sv", c(0.2, 3, 0.3, 0.5))
