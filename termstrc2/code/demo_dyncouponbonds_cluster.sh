#!/bin/sh
#$ -N dynbonds

R-g --vanilla < demo_dyncouponbonds_cluster.R

#$ -m bae
#$ -M robert.ferstl@wiwi.uni-regensburg.de
# -q node.q
