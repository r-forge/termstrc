## Package dependencies
library("rgl")
library("urca")
library("sandwich")

## Estimation kernel
source("estim_cs.R")
source("estim_nss_couponbonds.R")
source("estim_nss_dyncouponbonds.R")
source("estim_nss_zeroyields.R")
source("bondpricing.R")

## Parametric forms of spot rate functions
source("spotfwdratedef.R")
source("cubicfunc.R")

## Data handling
source("couponbonds_data.R")
source("create_cf_m.R")

## Methods for estimation results
source("methods_curves.R")
source("methods_dyntermstrc_nss.R")
source("methods_dyntermstrc_param.R")
source("methods_dyntermstrc_yields.R")
source("methods_termstrc_cs.R")
source("methods_termstrc_nss.R")
source("methods_zeroyields.R")
source("gof.R")


