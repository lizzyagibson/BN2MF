# Run one time!

## Packages
source("./functions/compare_functions.R")

# Combine R tibbles from BN2MF models
m_rank = tibble()
m_metrics = tibble()
m_prop = tibble()

job_num = 1
# originally ran this on 12100 jobs

for (j in 1:length(job_num)) {
  if (file.exists(paste0("./main/bn2mf/output/sep_prop_m_", j, ".RDA"))) {
    load(paste0("./main/bn2mf/output/sep_prop_m_", j, ".RDA"))
    m_prop = bind_rows(m_prop, vci_metrics)
  } else {m_prop = bind_rows(m_prop, vci_metrics = tibble(NA))}
  
  if (j %% 1000 == 0) {print(j)}
}

# This file has all 12100 rows
# save(m_prop, file = "./main/bn2mf/output/m_prop.RDA")

for (j in 1:length(job_num)) {
  load(paste0("./main/bn2mf/output/sep_rank_m_", j, ".RDA"))
  m_rank = bind_rows(m_rank, sep_rank)
  if (i %% 1000 == 0) {print(i)}
}

# This file has all 12100 rows
# save(m_rank, file = "./main/bn2mf/output/m_rank.RDA")

for (j in 1:length(job_num)) {
  load(paste0("./main/bn2mf/output/sep_metrics_m_", j, ".RDA"))
  m_metrics = bind_rows(m_metrics, sep_metrics)
  if (j %% 1000 == 0) {print(j)}
}

# save(m_metrics, file = "./main/bn2mf/output/m_metrics.RDA")
# This file has all 12100 rows

#### R models ####

## Combine other methods, main results in R

# Combine R tibbles
# Model output from PCA, FA, and NMF
other_rank = tibble()
other_metrics = tibble()

# Read these from `other_model_metrics.R`
for (j in 1:length(job_num)) {
  load(paste0("./main/other_models/output/other_rank_", j, ".RDA"))
  other_rank = bind_rows(other_rank, sep_rank)
  if (i %% 1000 == 0) {print(i)}
}

# This file has all 12100 rows
# save(other_rank, file = "./main/other_models/output/other_rank.RDA")

for (j in 1:length(job_num)) {
  load(paste0("./main/other_models/output/other_metrics_", j, ".RDA"))
  other_metrics = bind_rows(other_metrics, sep_metrics)
  if (j %% 1000 == 0) {print(j)}
}

# This file has all 12100 rows
# save(other_metrics, file = "./main/other_models/output/other_metrics.RDA")


