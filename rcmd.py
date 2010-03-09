#! /usr/bin/env python

import os
import shutil

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
