## Combine BN2MF with VCI, separate sim results in MATLAB
## Combine other methods, main results in R

## Take metrics for each method
# Symmetric subspace distance
# Cosine distance
# Relative error
# on sims
# on loadings
# on scores

# Run one time!

## Packages
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

## Combine R tibbles from `sep_other_models.R`
## Model output from PCA, FA, and NMF
# sep_r = tibble()
# for (i in 1:1210) {
#   load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_models_out/sep_out_", i, ".RDA"))
#   sep_r = bind_rows(sep_r, sep_out)
#   print(i)
# }

save(sep_r, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/combined_models.RDA")




