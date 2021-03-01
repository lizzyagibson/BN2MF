# Get bootstrap coverage

####  Load packages ####
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(R.matlab, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

#### For 600 sims ####

# 600 simulated data sets
datasets = read_csv("/ifs/scratch/msph/ehs/eag2186/Data/bs_ids.csv")
bs_prop = tibble()

for (i in 1:nrow(datasets)) {

  set = datasets$value[i]
  
  # BS 95% CI
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_lower_wa_", set, ".RDA"))
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_upper_wa_", set, ".RDA"))
  
  bs_lower_wa = bs_lower_wa[,2:5]
  bs_upper_wa = bs_upper_wa[,2:5]
  
  # True Scores
  true_scores = read_csv(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_scores/scores_sep_", set, ".csv"))
  # True patterns
  true_patterns = read_csv(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_patterns/patterns_sep_", set, ".csv"))
  
  #### Normalize truth ####
  denom = apply(true_patterns,1, sum)
  patterns_scaled = true_patterns/denom
  scores_scaled = as.matrix(true_scores) %*% diag(denom)
  
  prop = sum(scores_scaled >= bs_lower_wa & scores_scaled <= bs_upper_wa)/(nrow(scores_scaled)*ncol(scores_scaled))

  prop_row = tibble(set, prop)
  bs_prop = bind_rows(bs_prop, prop_row)
  print(paste0("Sim Number: ", i))
}

# Save prop
save(bs_prop, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bs_prop.RDA")

