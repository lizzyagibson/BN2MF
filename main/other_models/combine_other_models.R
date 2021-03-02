## Combine BN2MF with VCI, separate sim results in MATLAB
## Combine other methods, main results in R

# Run one time!

## Packages
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

# Combine R tibbles
# Model output from PCA, FA, and NMF
other_rank = tibble()
other_metrics = tibble()

# Read these from `other_model_metrics.R`
for (i in 1:12100) {
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_rank/sep_rank_", i, ".RDA"))
  other_rank = bind_rows(other_rank, sep_rank)
  if (i %% 100 == 0) {print(i)}
}

save(other_rank, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/other_rank.RDA")
rm(sep_rank)
rm(other_rank)

for (j in 1:12100) {
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_metrics/sep_metrics_", j, ".RDA"))
  other_metrics = bind_rows(other_metrics, sep_metrics)
  if (j %% 100 == 0) {print(j)}
}

save(other_metrics, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/other_metrics.RDA")



