library(tidyverse)

## read job number from system environment
## This only works if run on cluster!!
job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
job_num

load(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sim_dim/sim_dim", job_num, ".RDA"))
to_save

sim_t = to_save %>% dplyr::select(sim)

sim_m = sim_t[[1]][[1]]
sim_d = as_tibble(sim_m)

write_csv(sim_m, path = paste0("/ifs/scratch/msph/ehs/eag2186/Data/sim_csv/sim_dim", job_num, ".csv"))
