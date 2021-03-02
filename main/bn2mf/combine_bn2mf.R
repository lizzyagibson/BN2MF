## Combine BN2MF with VCI, separate sim results in MATLAB
## Combine other methods, main results in R

# Run one time!

## Packages
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/compare_functions.R")

# Combine R tibbles
# Model output from PCA, FA, and NMF
m_rank = tibble()
m_metrics = tibble()
m_prop = tibble()

for (j in 1:12100) {
  if (file.exists(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_prop_m/sep_prop_m_", j, ".RDA"))) {
    load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_prop_m/sep_prop_m_", j, ".RDA"))
    m_prop = bind_rows(m_prop, vci_metrics)
  } else {m_prop = bind_rows(m_prop, vci_metrics = tibble(NA))}
  
  if (j %% 100 == 0) {print(j)}
}

save(m_prop, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/m_prop.RDA")
rm(m_prop)
rm(vci_metrics)

for (i in 1:12100) {
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_rank_m/sep_rank_m_", i, ".RDA"))
  m_rank = bind_rows(m_rank, sep_rank)
  if (i %% 100 == 0) {print(i)}
}

save(m_rank, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/m_rank.RDA")
rm(sep_rank)
rm(m_rank)

for (j in 1:12100) {
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_metrics_m/sep_metrics_m_", j, ".RDA"))
  m_metrics = bind_rows(m_metrics, sep_metrics)
  if (j %% 100 == 0) {print(j)}
}

save(m_metrics, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/m_metrics.RDA")



