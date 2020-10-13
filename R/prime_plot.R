#####
library(ggsci)
library(tidyverse)
library(janitor)

#####
##### Data
#####

load("./HPC/Rout/all_unstand_error.RDA")
load("./HPC/Rout/all_unstand_ssdist.RDA")
load("./HPC/Rout/all_unstand.RDA")
load("./HPC/Rout/over_output_all_un.RDA")

#####
##### Settings
#####

theme_set(theme_bw(base_size = 26) + 
            theme(strip.background = element_rect(fill="white"),
                  axis.text.x = element_blank(), # element_text(angle = 45, hjust = 1, size = 30),
                  #legend.position = "none",
                  legend.direction = "horizontal",
                  legend.position = c(0.5, -0.05), # c(0,0) bottom left, c(1,1) top-right.
                  axis.title.y = element_text(size = 40),
                  plot.margin = unit(c(5.5, 30, 40, 5.5), "points")) # top, then right, bottom and left. 
          )

ggsci <- pal_nejm()(8)

options(
  ggplot2.discrete.colour = ggsci,
  ggplot2.discrete.fill = ggsci)

#####
##### Clean
#####

error <- all_unstand_error %>% 
  filter(name != "nmf_p_pred") %>% 
  mutate(sim_type = paste0(str_to_title(sim_type), " Sims")) %>% 
  mutate(name = case_when(name == "pca_pred" ~ "PCA",
                          name == "fa_pred" ~ "Factor\nAnalysis",
                          name == "nmf_l2_pred" ~ "NMF",
                          name == "bnmf_pred" ~ "BN2MF")) %>% 
  mutate(sim_type = fct_relevel(sim_type, "Distinct Sims",
                                "Overlapping Sims",
                                "Correlated Sims")) %>% 
  mutate(name = fct_relevel(name, "BN2MF",
                            "PCA",
                            "Factor\nAnalysis",
                            "NMF"))

ssdist <- 
  all_unstand_ssdist %>% 
  filter(model != "nmf_p" & model != "NMF_p") %>% 
  mutate(sim_type = paste0(str_to_title(sim_type), " Sims")) %>% 
  mutate(model = case_when(model %in% c("pca","PCA") ~ "PCA",
                           model %in% c("fa","FA") ~ "Factor\nAnalysis",
                           model %in% c("nmf_l2","NMF_l2") ~ "NMF",
                           model %in% c("bnmf","npBNMF") ~ "BN2MF")) %>% 
  mutate(sim_type = fct_relevel(sim_type, "Distinct Sims",
                                "Overlapping Sims",
                                "Correlated Sims")) %>% 
  mutate(model = fct_relevel(model, "BN2MF",
                             "PCA",
                             "Factor\nAnalysis",
                             "NMF"))

#####
##### Symmetric Subspace Distance
#####

pdf("Figures/prime_load.pdf", height = 12)
ssdist %>% 
  filter(type == "Loadings" & Normalized == "No" & 
           sim_type == "Overlapping Sims" & Standardized == "No") %>% 
  ggplot(aes(x = model, y = value, group = model, 
             color = model, fill = model)) +
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", alpha = 0.5) +
  geom_boxplot(alpha = 0.5, width = 1) +
  lims(y = c(0,1)) +
  labs(title = "Distance from True Loadings",
       y = "Symmetric Subspace Distance", x = "",
       color = "", fill = "")
dev.off()

pdf("Figures/prime_score.pdf", height = 12)
ssdist %>% 
  filter(type == "Scores" & Normalized == "No" & 
           sim_type == "Overlapping Sims" & Standardized == "No") %>% 
  ggplot(aes(x = model, y = value, group = model, 
             color = model, fill = model)) +
  geom_hline(yintercept = 0.5, color = "red", alpha = 0.5, linetype = "dashed") +
  geom_boxplot(alpha = 0.5, width = 1, notch = TRUE) +
  lims(y = c(0,1)) +
  labs(title = "Distance from True Scores",
       y = "Symmetric Subspace Distance", x = "",
       color = "", fill = "")
dev.off()

#####
##### Predictive Error
#####

error %>% 
  filter(sim_type == "Overlapping Sims") %>% 
  ggplot(aes(x = name, y = l2_true)) +
  #geom_boxplot(alpha = 0.5, width = 1) +
  geom_jitter(height = 0, width = 0.25, size = 0.5, alpha = 0.5) +
  geom_violin(aes(fill = name), scale = "width", alpha = 0.5, # Scale maximum width to 1 for all violins
              draw_quantiles = c(0.5)) + 
  labs(title = "Distance from Truth",
       y = "Relative Prediction Error", x = "",
       color = "", fill = "") +
  scale_y_log10()

pdf("Figures/prime_error.pdf", height = 12)
error %>% 
  filter(sim_type == "Overlapping Sims") %>% 
  ggplot(aes(x = name, y = l2_true, group = name, 
             color = name, fill = name)) +
  geom_boxplot(alpha = 0.5, width = 1) +
  labs(title = "Distance from Truth",
       y = "Relative Prediction Error", x = "",
       color = "", fill = "") +
  scale_y_log10()
dev.off()

#####
##### Rank
#####

all_unstand %>% 
  filter(sim_type == "overlapping") %>% 
  dplyr::select(seed, grep("rank", colnames(.)), eh) %>%
  mutate(bnmf_rank = map(eh, nrow)) %>% 
  pivot_longer(pca_rank:bnmf_rank) %>% 
  filter(name != "eh") %>% 
  unnest(value) %>% 
  group_by(name, value) %>% 
  summarise(n = n()) %>%
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_all(replace_na, 0) %>% 
  dplyr::select(name, 4, 5, 2, 3) %>% 
  stargazer(summary=FALSE)
