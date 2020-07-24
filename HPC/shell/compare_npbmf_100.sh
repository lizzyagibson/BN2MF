#!/bin/bash
#$ -cwd -S /bin/bash
#$ -l mem=8G
#$ -l time=:2000:
#$ -M eag2186@cumc.columbia.edu

$MODULESHOME/init/bash
module load R/3.6.0

clear

R CMD BATCH everything_on_cluster_100.R


