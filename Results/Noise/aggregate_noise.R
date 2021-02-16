# Combine .RDA files from "combine_noise.R"
# into single dataframe
# This runs 1 time

# Packages
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

## Get functions
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/compare_functions.R")
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/factor_correspondence.R")

# Read in metrics for
# PCA, FA, NMF
n = 300

noise = tibble()
for (i in 1:n) {
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/combo_noise/noise_out_", i, ".RDA"))
  output_all = output_all %>% mutate(id = i)
  noise = bind_rows(noise, output_all)
}

# Read in BN2MF metrics
bn2mf_n = tibble()
for (i in 1:n) {
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/combo_bnmf_noise/bnmf_noise_out_", i, ".RDA"))
  output_all = output_all
  bn2mf_n = bind_rows(bn2mf_n, noise)
}

noise_all = left_join(noise, bn2mf_n) %>% 
  mutate(patterns = paste0(patterns, " Patterns"),
         patterns = ifelse(patterns == "1 Patterns", "1 Pattern", patterns),
         patterns = fct_relevel(patterns, "4 Patterns", after = 1),
         participants = paste0(participants, " Participants"),
         chemicals = paste0(chemicals, " Chemicals"),
         chemicals = fct_inorder(chemicals))

save(noise_all, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/sup_noise.RDA")
