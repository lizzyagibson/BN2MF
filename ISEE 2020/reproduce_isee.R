load(file = "./Sims/sim_over_old_1012.RDA")
load(file = "./Sims/sim_dist_old_1012.RDA")
load(file = "./Sims/sim_cor_old_1012.RDA")

sim_dist_1012 <- sim_dist
sim_over_1012 <- sim_over
sim_cor_1012 <- sim_cor

load(file = "./Sims/sim_over_old.RDA")
load(file = "./Sims/sim_dist_old.RDA")
load(file = "./Sims/sim_cor_old.RDA")

all.equal(sim_dist, sim_dist_1012)
all.equal(sim_over, sim_over_1012)
all.equal(sim_cor, sim_cor_1012)

