#$ -N termstrc -q bignode.q

R-g --vanilla < test_cluster_estim_ns.R

#$ -M rferstl@wu.ac.at

echo "#### finished ####"

