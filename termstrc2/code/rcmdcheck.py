#! /usr/bin/env python

import os
import shutil

## Estimation kernel
os.system("cp estim_cs.R ~/WORK/termstrc/pkg/R")
os.system("cp estim_nss_couponbonds.R ~/WORK/termstrc/pkg/R")
os.system("cp estim_nss_dyncouponbonds.R ~/WORK/termstrc/pkg/R")
os.system("cp estim_nss_zeroyields.R ~/WORK/termstrc/pkg/R")
os.system("cp bondpricing.R ~/WORK/termstrc/pkg/R")

## Parametric forms of spot rate functions
os.system("cp spotfwdratedef.R ~/WORK/termstrc/pkg/R")
os.system("cp cubicfunc.R ~/WORK/termstrc/pkg/R")

## Data handling
os.system("cp couponbonds_data.R ~/WORK/termstrc/pkg/R")
os.system("cp create_cf_m.R ~/WORK/termstrc/pkg/R")
os.system("cp methods_couponbonds.R ~/WORK/termstrc/pkg/R")

## Methods for estimation results
os.system("cp methods_curves.R ~/WORK/termstrc/pkg/R")
os.system("cp methods_dyntermstrc_nss.R ~/WORK/termstrc/pkg/R")
os.system("cp methods_dyntermstrc_param.R ~/WORK/termstrc/pkg/R")
os.system("cp methods_dyntermstrc_yields.R ~/WORK/termstrc/pkg/R")
os.system("cp methods_termstrc_cs.R ~/WORK/termstrc/pkg/R")
os.system("cp methods_termstrc_nss.R ~/WORK/termstrc/pkg/R")
os.system("cp methods_zeroyields.R ~/WORK/termstrc/pkg/R")
os.system("cp gof.R ~/WORK/termstrc/pkg/R")

## Copy Rd files
os.system("cp -R  ~/WORK/termstrc/termstrc2/helpfiles/*.Rd  ~/WORK/termstrc/pkg/man")

## Perform R CMD check
os.system("R CMD check ~/WORK/termstrc/pkg")

## Print Log file
logfile = os.path.join(os.getcwd(),"pkg.Rcheck/00install.out")

inp = open(logfile,"r")
for line in inp.readlines():
    print(line)
    
inp.close()

os.system("R CMD build ~/WORK/termstrc/pkg")
os.system("R CMD INSTALL termstrc_1.2.tar.gz")

