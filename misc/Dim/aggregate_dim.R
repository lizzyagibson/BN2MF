# Combine .RDA files from "combine_dim.R"
# into single dataframe
# This runs 1 time

# Packages
#library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

## Get functions
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/compare_functions.R")
#source("/ifs/scratch/msph/ehs/eag2186/npbnmf/factor_correspondence.R")

# Read in metrics for
# PCA, FA, NMF
n = 13500

dim_r_metrics = tibble()
for (i in 1:n) {
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/dim/combo_r_dim/dim_out_", i, ".RDA"))
  # output_all = output_all %>% mutate(id = i)
  dim_r_metrics = bind_rows(dim_r_metrics, output_all)
  if (i %% 100 == 0) {print(i)}
}
save(dim_r_metrics, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/dim/dim_r_metrics.RDA")

# Read in BN2MF metrics
bn2mf_rank = tibble()
bn2mf_metrics = tibble()
bn2mf_prop = tibble()

for (j in 1:n) {
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/dim/dim_prop_m/dim_prop_m_", j, ".RDA"))
  bn2mf_prop = bind_rows(bn2mf_prop, vci_metrics)
  
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/dim/dim_rank_m/dim_rank_m_", j, ".RDA"))
  bn2mf_rank = bind_rows(bn2mf_rank, dim_rank)
  
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/dim/dim_metrics_m/dim_metrics_m_", j, ".RDA"))
  bn2mf_metrics = bind_rows(bn2mf_metrics, dim_metrics)
  
  if (j %% 100 == 0) {print(j)}
}

save(bn2mf_metrics, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/dim/dim_m_metrics.RDA")
save(bn2mf_prop, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/dim/dim_m_prop.RDA")
save(bn2mf_rank, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/dim/dim_m_rank.RDA")
