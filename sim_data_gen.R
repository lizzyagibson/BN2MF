# Simulate data-generating process
# 6/21/2019
options(scipen = 999)
library(MNdata)
library(reshape2)
library(compositions)
library(MASS)
library(tidyverse)
library(gridExtra)

# Simulate dataset with 1000 individuals and 35 chemicals to approx Mothers and Newborns cohort data.

#####################
## MN Data #########
#####################

# 8 phenols, 9 phthalates, 10 pbdes, 8 pcbs
# phenols and phthalates are specific gravity adjusted
# all values < LOD are replaced by participant-specific LOD/sqrt(2)

mn_data <- mn_edc %>% dplyr::select(1:52) %>%
  select_if(~n_distinct(.) > 1) %>% #drop column if all values are the same
  dplyr::select(1:18, pcb105, pcb74, pcb99, pcb118,
                pcb138_158, pcb153, pcb187, pcb180,
                grep("BDE", colnames(.))) %>%
  dplyr::select(-sid)

