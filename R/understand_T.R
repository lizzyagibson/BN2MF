library(tidyverse)
library(R.matlab)

grid_out1 = readMat("./MATLAB/grid_out.m")[[1]] %>% as_tibble()
grid_out2 = readMat("./MATLAB/grid_out2.m")[[1]] %>% as_tibble()

grid_out = rbind(grid_out1, grid_out2)

colnames(grid_out) = c("t0", "base", "power", "rank", "score", "var", "eh", "var_eh", "iter")

grid_out %>% 
  #filter(t0 == 0 & power ==1) %>% 
  filter(rank == 4 & base == 0.75 & power == 1) %>% 
  arrange(desc(score))

grid_out %>% 
  mutate(rank4 = ifelse(rank == 4, 1, 0)) %>% 
  group_by(rank4) %>% 
  summarize(n = n(),
            t0_max = max(t0))

df_out %>% 
  filter(basepow == 0.99) %>% 
  ggplot(aes(x = iter, y = powC)) +
  geom_tile(aes(fill = temp)) +
  scale_fill_viridis()

df_out %>% filter(temp <= 1.01) %>% 
  ggplot(aes(x = basepow, y = powC)) +
  geom_tile(aes(fill = iter)) +
  scale_fill_viridis()

#####
# Mean diff in variance
mean_out1 = readMat("./MATLAB/grid_out3.m")[[1]] %>% as_tibble() %>% mutate(run = 1)
mean_out2 = readMat("./MATLAB/grid_out4.m")[[1]] %>% as_tibble() %>% mutate(run = 2)
mean_out = rbind(mean_out1, mean_out2)
colnames(mean_out) = c("t0", "rank", "score", "iter", "mean", "run")

mean_out %>% 
  filter(rank != 4) %>% 
  arrange(desc(rank))

mean_o = readMat("./MATLAB/grid_out4.mat")[[1]] %>% as_tibble()
colnames(mean_o) = c("t0", "rank", "score", "iter", "mean_diff", "mean_T", "mean_reg")

mean_o %>% 
  filter(rank == 4) %>% 
  arrange((t0))

# Here we ran BN2MF with annealing 10x and chose the best
mean_10 = readMat("./MATLAB/grid_out5.mat")[[1]] %>% as_tibble()
colnames(mean_10) = c("t0", "rank", "score", "iter", "mean_diff", "mean_T", "mean_reg")

mean_10 %>% 
  filter(rank == 4) %>% 
  arrange(desc(t0))

