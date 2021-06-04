pat1=c(0,0,5,10)
pat2=c(0,5,10,0)
pat3=c(5,10,0,0)
pat4=c(10,0,0,5)

library(textshape)
library(LearnBayes)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggridges)

source("./functions/compare_functions_2.R")

load(file = "./sims/sim_sep.RDA")
sim_sep

reps = 5
# `reps` is how overlapping the patterns are
# reps = 0 -- completely overlapping
# reps = 10 completely distinct
start = c()
if (reps != 0) { 
  for (i in 1:reps) {
    start = cbind(start,diag(1,4))
  }}

patterns = cbind(start, t(rdirichlet(10-reps, pat1)), t(rdirichlet(10-reps, pat2)),
                 t(rdirichlet(10-reps, pat3)), t(rdirichlet(10-reps, pat4)))
(patterns)

p5 = as_tibble(cluster_matrix(patterns))
colnames(p5) = 1:ncol(p5)

# pdf("./figures/loadings_plot5.pdf", width = 10, height = 10)
p5 %>%
  mutate(Pattern = 1:nrow(.),
         Pattern = str_c("Pattern ", Pattern)) %>%
  pivot_longer(1:40) %>%
  mutate(name = fct_inorder(str_remove(name, "V"))) %>%
  ggplot(aes(x = name, y = value)) +
  geom_col(fill = "orange", color = "black") +
  facet_wrap(.~Pattern) +
  labs(y = "Simulated Loadings", x = "Simulated Chemicals") +
  theme_light(base_size = 45) +
  theme(panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size = 45, colour = 'black'))
# dev.off()

# Ridges ####
# Viz whole distributions
# VCI and bootstrap
# For 1 example sim

# Read in distributions ####

# Truth ####
chem = read_csv("./misc/ci_ex/true_chem.csv")
patterns = read_csv("./misc/ci_ex/true_patterns.csv")
scores = read_csv("./misc/ci_ex/true_scores.csv")

denom = apply(patterns, 1, sum)
patterns_scaled = patterns / denom
scores_scaled = as.matrix(scores) %*% diag(denom)

# Variational distributions ####

vci_pred      = readMat("./misc/vci_ex/q_pred_9562.mat")[[1]]
vci_eh_dist   = readMat("./misc/vci_ex/sep_distEH_9562.mat")[[1]]
vci_ewa_dist  = readMat("./misc/vci_ex/sep_distWA_9562.mat")[[1]]
vci_eh        = readMat("./misc/vci_ex/sep_eh_scaled9562.mat")[[1]]
vci_ewa       = readMat("./misc/vci_ex/sep_ewa_scaled9562.mat")[[1]]
vci_eh_lower  = readMat("./misc/vci_ex/sep_lowerH_9562.mat")[[1]]
vci_ewa_lower = readMat("./misc/vci_ex/sep_lowerWA_9562.mat")[[1]]
vci_eh_upper  = readMat("./misc/vci_ex/sep_upperH_9562.mat")[[1]]
vci_ewa_upper = readMat("./misc/vci_ex/sep_upperWA_9562.mat")[[1]]
dim(vci_ewa_dist)

# Plot VCI ####
theme_set(theme_bw(base_size = 20) + 
            theme(strip.background = element_rect(fill="white"),
                  legend.direction = "horizontal",
                  legend.title = element_blank()) # top, then right, bottom and left. 
)

# Read output
dim(vci_ewa_dist)

dim(vci_ewa_lower)
dim(vci_ewa_upper)
dim(vci_ewa)

lower = vci_ewa_lower %>% 
  as_tibble() %>% 
  mutate(id = 1:nrow(.)) %>% 
  pivot_longer(V1:V4,
               names_to = "pattern",
               values_to = "lower")

ewa = vci_ewa %>% 
  as_tibble() %>% 
  mutate(id = 1:nrow(.)) %>% 
  pivot_longer(V1:V4,
               names_to = "pattern",
               values_to = "ewa")

upper = vci_ewa_upper %>% 
  as_tibble() %>% 
  mutate(id = 1:nrow(.)) %>% 
  pivot_longer(V1:V4,
               names_to = "pattern",
               values_to = "upper")

vci = full_join(lower, ewa) %>% full_join(., upper) %>% 
  mutate(pattern = str_replace(pattern, "V", "P"))

vci %>% group_by(pattern) %>% 
  summarise(max = max(ewa),
            min = min(ewa))
vci %>% 
  arrange((ewa))

# Plot distribution ####
dist_1 =  vci_ewa_dist[,1,] %>% 
  as_tibble() %>% 
  mutate(id = as_factor(1:nrow(.))) %>% 
  pivot_longer(V1:V1000) %>% 
  mutate(pattern = "Pattern 1")

dist_2 =  vci_ewa_dist[,2,] %>% 
  as_tibble() %>% 
  mutate(id = as_factor(1:nrow(.))) %>% 
  pivot_longer(V1:V1000) %>% 
  mutate(pattern = "Pattern 2")

dist_3 =  vci_ewa_dist[,3,] %>% 
  as_tibble() %>% 
  mutate(id = as_factor(1:nrow(.))) %>% 
  pivot_longer(V1:V1000) %>% 
  mutate(pattern = "Pattern 3")

dist_4 =  vci_ewa_dist[,4,] %>% 
  as_tibble() %>% 
  mutate(id = as_factor(1:nrow(.))) %>% 
  pivot_longer(V1:V1000) %>% 
  mutate(pattern = "Pattern 4")

vci_dist = bind_rows(dist_1, dist_2, dist_3, dist_4)

vci_dist %>% group_by(pattern) %>% 
  summarise(max = max(value),
            min = min(value))

sum(scores_scaled <= vci_ewa_upper & scores_scaled >= vci_ewa_lower)/4000

as_tibble(scores_scaled <= vci_ewa_upper & scores_scaled >= vci_ewa_lower) %>% 
  mutate(row = 1:nrow(.)) %>% 
  filter(!V1 & !V3)

vci_ewa_lower[183,]
scores_scaled[183,]
vci_ewa_upper[183,]

scores_scaled[c(107, 170, 183, 336, 340, 412),] %>% as_tibble()

trueall = scores_scaled[c(333, 712, 56, 183),] %>% as_tibble()

truths = bind_cols(tibble(idx= as_factor(LETTERS[1:4])), trueall) %>% 
  pivot_longer(V1:V4,
               names_to = "pattern",
               values_to = "true") %>% 
  mutate(pattern = str_replace(as_factor(pattern), "V", "Pattern "))

ridges = vci_dist %>% 
  filter(id %in% c(333, 712, 56, 183)) %>%
  arrange(pattern, id) %>% 
  mutate(idx = as_factor(case_when(id == 333 ~ "A",
                                   id == 183 ~ "D",
                                   id == 712 ~ "B",
                                   id == 56 ~ "C"))) %>% 
  left_join(., truths) %>% 
  mutate(idx = fct_relevel(idx, "D", "C", "B", "A")) %>% 
  ggplot(aes(x = value, y = idx)) +
  geom_density_ridges(aes(fill=pattern), bandwidth=2, scale = 1, alpha = 0.6,
                      quantile_lines=TRUE,
                      #quantiles = c(0.025, 0.5, 0.975),
                      quantile_fun=function(x,...) 
                      {return(c(quantile(x, 0.025), mean(x), quantile(x, 0.975)))}
  ) + 
  geom_segment(aes(x = true, xend = true, y = as.numeric(idx), color = pattern,
                   yend = as.numeric(idx) + 1), linetype = "dashed") +
  theme(axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title = element_text(size=30),
        legend.text = element_text(size = 18),
        # legend.position = "bottom",
        legend.direction = "vertical",
        legend.position = c(0.85, 0.1275), # c(1,0) right bottom, c(1,1) right top.
        ) +
  guides(fill = guide_legend(override.aes = list(linetype = 0))) +
  labs(y ="Simulated individuals",
       x = "Pattern scores") +
  scale_fill_nejm()

ggsave("./figures/ridges.pdf", plot=ridges,  width=7, height=7)

# ERROR ####
# Defense figures

# Example patterns half overlapping
half_sep = sim_sep %>% 
  filter(sep_num %in% c(0, 5, 10))

half_sep %>% slice(2003)

p1 = half_sep$true_patterns[[1]]
p_re1 = cluster_matrix(p1, dim="col")
p1 = as_tibble(p_re1)
colnames(p1) = 1:ncol(p1)
p1long = p1 %>%
  mutate(Pattern = 1:nrow(.),
         Pattern = str_c("Pattern ", Pattern)) %>%
  pivot_longer(1:40) %>%
  mutate(name = fct_inorder(name),
         which = "Overlapping")

p2 = half_sep$true_patterns[[1002]]
p_re2 = cluster_matrix(p2)
p2 = as_tibble(p_re2)
colnames(p2) = 1:ncol(p2)
p2long = p2 %>%
  mutate(Pattern = 1:nrow(.),
         Pattern = str_c("Pattern ", Pattern)) %>%
  pivot_longer(1:40) %>%
  mutate(name = fct_inorder(name),
         which = "Hybrid")

p3 = half_sep$true_patterns[[2003]]
p_re3 = cluster_matrix(p3)
p3 = as_tibble(p_re3)
colnames(p3) = 1:ncol(p3)
p3long = p3 %>%
  mutate(Pattern = 1:nrow(.),
         Pattern = str_c("Pattern ", Pattern)) %>%
  pivot_longer(1:40) %>%
  mutate(name = fct_inorder(name),
         which = "Distinct")

#pdf("./figures/loadings_plot_bnmf.pdf", width=10, height=10)
bind_rows(p1long, p2long, p3long) %>% 
  filter(Pattern == "Pattern 1" & which!="Hybrid") %>% 
  filter(value != 0) %>% 
  ggplot(aes(x = name, y = value)) +
  geom_col(fill = "orange", color = "black") +
  facet_wrap(~which) +
  labs(y = "Simulated Loadings", x = "Simulated Chemicals") +
  theme_light(base_size = 45) +
  theme(panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size = 45, colour = 'black'))
#dev.off()

# example correlation matrix, no noise, patterns half overlapping
brewer.pal(8, "Spectral")
display.brewer.pal(8, "Spectral")

get_lower_tri <-function(x){
  x[upper.tri(x)] <- NA
  return(x)
}

corr = cor(half_sep$sim[[1]], method = "spearman")

corr_re = cluster_matrix(corr)

cormat <- as.data.frame(get_lower_tri(corr_re)) %>% 
  rownames_to_column(var = "Chem") %>% 
  as_tibble() %>% 
  mutate(Chem = fct_inorder(Chem),
         Chem = as_factor(1:nrow(.))) %>% 
  pivot_longer(V1:V40) %>% 
  mutate(name = fct_inorder(name),
         name = as_factor(rep(1:40,40)))

corplot = cormat %>%
  mutate(outline = ifelse(is.na(value), FALSE, TRUE))
corplot$outline[!corplot$outline] <- NA

# pdf("./figures/sim_corr_bn2mf.pdf")
corplot %>% 
  ggplot(aes(x = Chem, y = name, fill = value)) + 
  geom_tile() +
  geom_tile(data = corplot[!is.na(corplot$outline), ], aes(color = outline), size = 0.75) +
  scale_color_manual(guide = FALSE, values = c(`TRUE` = "black")) +
  labs(x = "Simulated Chemicals", y = "", fill = "") +
  theme_light(base_size = 25) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.075,0.85),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=10, angle=45, hjust = 1),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size = 25, colour = 'black')) +
  scale_fill_distiller(palette="Spectral",
                       limits=c(-1,1),
                       na.value = 'white')+
  #scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE)
# dev.off()

# relative error
load(file="./main/metrics.rda")

metrics = get_metrics %>% 
  filter(#sep_num %in% c(0, 10) & 
    noise_level %in% c(0.2, 0.5)) %>% 
  mutate(#sep = ifelse(sep_num == 10, "Distinct Patterns", "Overlapping Patterns"),
    noise = case_when(noise_level == 0.2 ~ "Noise + 20%",
                      noise_level == 0.5 ~ "Noise + 50%", 
                      TRUE ~ "Noise +100%"))

#pdf("./figures/prez_error.pdf", width=9)
metrics %>% 
  filter(matrix == "Pred") %>% 
  mutate(noise = fct_inorder(noise)) %>% 
  ggplot(aes(x = model, y = relerr, fill = model, col = model)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.5, notch=TRUE) +
  #geom_violin(alpha=0.5, scale="count", draw_quantiles=T) +
  facet_grid(.~noise) +
  scale_y_log10() + 
  theme_bw(base_size = 25) + 
  labs(x = "", y = "Relative Prediction Error", fill = "", col = "") + 
  scale_color_discrete(labels = c(expression(paste("B", N^2, "MF")),'FA','NMF-L2', 'NMF-P', 'PCA')) +
  scale_fill_discrete(labels = c(expression(paste("B", N^2, "MF")),'FA','NMF-L2', 'NMF-P', 'PCA')) +
  theme(axis.text.x = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.08), # c(0,0) bottom left, c(1,1) top-right.
        legend.text = element_text(size = 15))
# dev.off()
