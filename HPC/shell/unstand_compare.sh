#!/bin/bash
#$ -cwd -S /bin/bash
#$ -o /dev/hull
#$ -e /dev/null
#$ -l mem=16G
#$ -l time=:10000:
#$ -M eag2186@cumc.columbia.edu

$MODULESHOME/init/bash
module load R/3.6.0

clear

R CMD BATCH --no-save everything_on_cluster_unstand.R


