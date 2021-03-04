# Bootstrap NMF

## Run script 1:90000 (150 bootstraps * 600 datasets)

## Get job number
job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
job_num

## Get bootstrap ids
bs_ids = read.csv("/ifs/scratch/msph/ehs/eag2186/Data/bs_ids.csv")
dim(bs_id)

bootstraps = 1:150

grid = expand.grid(bootstraps, bs_id[,1])
head(grid)

boot  = grid[job_num,1]
boot

sim_num  = grid[job_num,2]
sim_num

# If output exists, quit R
path = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf/"
if (file.exists(paste0(path, "nmf_bs_sim_", sim_num, "_bs_", boot,  ".RDA"))) {quit(save="no")}

## Get functions
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/compare_functions.R")

## Get sims

#chem <- read_csv(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_chem/chem_sep_", sim_num, ".csv"))

truescores <- read_csv(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_scores/scores_sep_", sim_num, ".csv")) %>% as.matrix()

sim <- read_csv(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_sim/sim_sep_", sim_num, ".csv")) %>% as.matrix()

truepatterns <- read_csv(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_patterns/patterns_sep_", sim_num, ".csv")) %>% as.matrix() 

## Normalize truth
patterns_denom  = apply(truepatterns, 1, sum)
truepatterns_scaled = truepatterns / patterns_denom
truescores_scaled   = truescores %*% diag(patterns_denom)

## Take a bootstrapped sample
n = nrow(sim)
p = ncol(sim)
sam = sample(1:n, n, replace = TRUE)
sim_sample = sim[sam,]
 
## Run nmf
nmf_sample <- nmf(sim_sample, 4, nrun = 100, method = "brunet")

scores <- basis(nmf_sample)
loadings <- coef(nmf_sample)
pred <- scores %*% loadings

## Normalize loadings to L1 norm across chemicals
denom = apply(loadings, 1, sum)
loadings_scaled = loadings / denom

## Scale scores by corresponding normalization constant
scores_scaled = scores %*% diag(denom)

## Rearrange solution to match truth
fc_out = factor_correspondence(t(truepatterns_scaled), t(loadings_scaled))
loadings_perm = t(fc_out$rearranged)
scores_perm = scores_scaled %*% fc_out$permutation_matrix

loadings_final = as_tibble(loadings_perm) %>% nest(data = everything()) %>% rename(loadings_final = data)
scores_sam   = as_tibble(cbind(sam, scores_perm)) %>% nest(data = everything()) %>% rename(scores_final = data)
pred_sam = as_tibble(cbind(sam, pred)) %>% nest(data = everything()) %>% rename(pred = data)

nmf_boot = bind_cols(id = sim_num, bootstrap = boot, loadings_final, scores_sam, pred_sam)

path = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf/"
save(nmf_boot, file = paste0(path, "nmf_bs_sim_", sim_num, "_bs_", boot,  ".RDA"))





