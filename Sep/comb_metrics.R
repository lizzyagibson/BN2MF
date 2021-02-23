## Combine BN2MF with VCI, separate sim results in MATLAB
## Combine other methods, main results in R

# Run one time!

## Packages
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

# Combine R tibbles
# Model output from PCA, FA, and NMF
all_rank = tibble()
all_metrics = tibble()

for (i in 1:12100) {
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_rank/sep_rank_", i, ".RDA"))
  all_rank = bind_rows(all_rank, sep_rank)
  if (i %% 100 == 0) {print(i)}
}

save(all_rank, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/all_rank.RDA")
rm(sep_rank)
rm(all_rank)

for (j in 1:12100) {
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_metrics/sep_metrics_", j, ".RDA"))
  all_metrics = bind_rows(all_metrics, sep_metrics)
  if (i %% 100 == 0) {print(j)}
}

save(all_metrics, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/all_metrics.RDA")



