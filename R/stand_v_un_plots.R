#####
library(ggsci)
library(tidyverse)
library(janitor)

load("./HPC/Rout/all_stand_error.RDA")
load("./HPC/Rout/all_unstand_error.RDA")
load("./HPC/Rout/all_stand_ssdist.RDA")
load("./HPC/Rout/all_unstand_ssdist.RDA")

#####
##### Data
#####

all_error <- full_join(all_stand_error, all_unstand_error) %>% 
  mutate(sim_type = paste0(str_to_title(sim_type), " Sims")) %>% 
  mutate(name = case_when(name == "pca_pred" ~ "PCA",
                          name == "fa_pred" ~ "Factor Analysis",
                          name == "nmf_l2_pred" ~ "NMF L2",
                          name == "bnmf_pred" ~ "npBNMF",
                          name == "nmf_p_pred" ~ "NMF Poisson")) %>% 
  mutate(sim_type = fct_relevel(sim_type, "Distinct Sims",
                                "Overlapping Sims",
                                "Correlated Sims")) %>% 
  mutate(name = fct_relevel(name, "npBNMF",
                            "PCA",
                            "Factor Analysis",
                            "NMF L2",
                            "NMF Poisson"))

all_ssdist <- 
  full_join(all_stand_ssdist, all_unstand_ssdist) %>% 
  mutate(sim_type = paste0(str_to_title(sim_type), " Sims")) %>% 
  mutate(model = case_when(model %in% c("pca","PCA") ~ "PCA",
                           model %in% c("fa","FA") ~ "Factor Analysis",
                           model %in% c("nmf_l2","NMF_l2") ~ "NMF L2",
                           model %in% c("bnmf","npBNMF") ~ "npBNMF",
                          model %in% c("nmf_p","NMF_p") ~ "NMF Poisson")) %>% 
  mutate(sim_type = fct_relevel(sim_type, "Distinct Sims",
                                "Overlapping Sims",
                                "Correlated Sims")) %>% 
  mutate(model = fct_relevel(model, "npBNMF",
                            "PCA",
                            "Factor Analysis",
                            "NMF L2",
                            "NMF Poisson"))

#####
##### Symmetric Subspace Distance
#####

#pdf("Figures/sim_ssdist_load_stand_v_un.pdf", width = 10)
all_ssdist %>% 
  filter(type == "Loadings" & Normalized == "No") %>% 
  ggplot(aes(x = model, y = value, group = interaction(model,Standardized), 
             color = model, fill = model, linetype = Standardized
             )) +
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed") +
  geom_boxplot(outlier.size = 0.25, width = 0.5, alpha = 0.5) +
  facet_wrap(~sim_type) + 
  theme_bw(base_size = 20) + 
  scale_y_log10() +
  theme(strip.text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill="white"),
        #legend.position = "bottom",
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.4), # c(0,0) bottom left, c(1,1) top-right.
        axis.title.y = element_text(size = 30)) + 
  labs(title = "Distance from True Pattern Loadings",
       y = "Symmetric Subspace Distance", x = "",
       color = "", fill = "") + scale_color_nejm() + scale_fill_nejm() +
  guides(color = FALSE, fill = FALSE)
#dev.off()

#pdf("Figures/sim_ssdist_score_stand_v_un.pdf", width = 10)
all_ssdist %>% 
  filter(type == "Scores" & Normalized == "No") %>% 
  ggplot(aes(x = model, y = value, group = interaction(model,Standardized), 
             color = model, fill = model, linetype = Standardized
  )) +
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed") +
  geom_boxplot(outlier.size = 0.25, width = 0.5, alpha = 0.5) +
  facet_wrap(~sim_type) + 
  theme_bw(base_size = 20) + 
  scale_y_log10() +
  theme(strip.text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill="white"),
        #legend.position = "bottom",
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.4), # c(0,0) bottom left, c(1,1) top-right.
        axis.title.y = element_text(size = 30)) + 
  labs(title = "Distance from True Individual Scores",
       y = "Symmetric Subspace Distance", x = "",
       color = "", fill = "") + scale_color_nejm() + scale_fill_nejm() +
  guides(color = FALSE, fill = FALSE)
#dev.off()

#####
##### Normalized Symmetric Subspace Distance
#####

#pdf("Figures/sim_ssdist_load_stand_v_un_norm.pdf", width = 10)
all_ssdist %>% 
  filter(type == "Loadings" & Normalized == "Yes") %>% 
  ggplot(aes(x = model, y = value, group = interaction(model,Standardized), 
             color = model, fill = model, linetype = Standardized
  )) +
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed") +
  geom_boxplot(outlier.size = 0.25, width = 0.5, alpha = 0.5) +
  facet_wrap(~sim_type) + 
  theme_bw(base_size = 20) + 
  scale_y_log10() +
  theme(strip.text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill="white"),
        #legend.position = "bottom",
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.4), # c(0,0) bottom left, c(1,1) top-right.
        axis.title.y = element_text(size = 30)) + 
  labs(title = "Distance from True Pattern Loadings",
       y = "Symmetric Subspace Distance", x = "",
       color = "", fill = "") + scale_color_nejm() + scale_fill_nejm() +
  guides(color = FALSE, fill = FALSE)
#dev.off()

#pdf("Figures/sim_ssdist_score_stand_v_un_norm.pdf", width = 10)
all_ssdist %>% 
  filter(type == "Scores" & Normalized == "Yes") %>% 
  ggplot(aes(x = model, y = value, group = interaction(model,Standardized), 
             color = model, fill = model, linetype = Standardized
  )) +
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed") +
  geom_boxplot(outlier.size = 0.25, width = 0.5, alpha = 0.5) +
  facet_wrap(~sim_type) + 
  theme_bw(base_size = 20) + 
  scale_y_log10() +
  theme(strip.text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill="white"),
        #legend.position = "bottom",
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.4), # c(0,0) bottom left, c(1,1) top-right.
        axis.title.y = element_text(size = 30)) + 
  labs(title = "Distance from True Individual Scores",
       y = "Symmetric Subspace Distance", x = "",
       color = "", fill = "") + scale_color_nejm() + scale_fill_nejm() +
  guides(color = FALSE, fill = FALSE)
#dev.off()

#####
##### Predictive Error
#####

#pdf("Figures/sim_error_wnoise_score_stand_v_un.pdf", width = 10)
all_error %>% 
  ggplot(aes(x = name, y = l2_sim, group = interaction(name,Standardized), 
             color = name, fill = name, linetype = Standardized
  )) +
  geom_boxplot(outlier.size = 0.25, width = 0.5, alpha = 0.5) +
  facet_wrap(~sim_type) + 
  theme_bw(base_size = 20) + 
  scale_y_log10(lim = c(0.019, 1), n.breaks = 8) +
  theme(strip.text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.4), # c(0,0) bottom left, c(1,1) top-right.
        axis.title.y = element_text(size = 30)
        ) + 
  labs(title = "Distance from Sims (True Obs + Noise)",
       y = "Relative Prediction Error", x = "",
       color = "", fill = "") + scale_color_nejm() + scale_fill_nejm() +
  guides(color = FALSE, fill = FALSE)
#dev.off()

#pdf("Figures/sim_error_truth_score_stand_v_un.pdf", width = 10)
all_error %>% 
  ggplot(aes(x = name, y = l2_true, group = interaction(name,Standardized), 
             color = name, fill = name, linetype = Standardized)) +
  geom_boxplot(outlier.size = 0.25, width = 0.5, alpha = 0.5) +
  facet_wrap(~sim_type) + 
  theme_bw(base_size = 20) + 
  scale_y_log10(lim = c(0.019,1), n.breaks = 8) +
  theme(strip.text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.4), # c(0,0) bottom left, c(1,1) top-right.
        axis.title.y = element_text(size = 30)
        ) + 
  labs(title = "Distance from Truth (No Noise)",
       y = "Relative Prediction Error", x = "",
       color = "", fill = "") + scale_color_nejm() + scale_fill_nejm() +
  guides(color = FALSE, fill = FALSE)
#dev.off()
