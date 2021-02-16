# Combine .RDA files from "combine_dim.R"
# into single dataframe
# This runs 1 time

# Packages
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

## Get functions
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/compare_functions.R")
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/factor_correspondence.R")

# Read in metrics for
# PCA, FA, NMF
n = 2700 # (3*3*3*100)

metrics = tibble()
for (i in 1:n) {
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/combo_dim/dim_out_", i, ".RDA"))
  output_all = output_all %>% mutate(id = i)
  metrics = bind_rows(metrics, output_all)
}

# Read in BN2MF metrics
bn2mf_m = tibble()
for (i in 1:n) {
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/combo_bnmf/bnmf_dim_out_", i, ".RDA"))
  output_all = output_all
  bn2mf_m = bind_rows(bn2mf_m, output_all)
}

dim_all = left_join(metrics, bn2mf_m) %>% 
  mutate(patterns = paste0(patterns, " Patterns"),
         patterns = ifelse(patterns == "1 Patterns", "1 Pattern", patterns),
         patterns = fct_relevel(patterns, "4 Patterns", after = 1),
         participants = paste0(participants, " Participants"),
         chemicals = paste0(chemicals, " Chemicals"),
         chemicals = fct_inorder(chemicals))

save(dim_all, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/sup_dim.RDA")
