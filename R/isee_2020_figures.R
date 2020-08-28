library(ggsci)
library(MNdata)
library(tidyverse)
library(R.matlab)
library(janitor)

################################
################################

labels <-  c("pcb167",
             "pcb170", "pcb178", "pcb183", "pcb187", "pcb180", "pcb189", "pcb194", "pcb196_203", "pcb199", "pcb206", "pcb209",
             "BDE17", "BDE28", "BDE47", "BDE66", "BDE85", "BDE99", "BDE100", "BDE153", "BDE154", "BDE183", "BDE209", "MECPP",
             "MEHHP", "MEOHP", "MCPP", "MIBP", "MBP", "MBZP", "MEP", "MEHP", "dcp_24", "dcp_25", "b_pb", "bp_3", "m_pb",
             "p_pb", "tcs", "bpa")

create_patterns_cor <- function (seed) {
  set.seed(seed)
  patterns_cor <- rbind(c(runif(12, 1, 2),   runif(11, 0.5, 1), runif(17, 0, .2)),
                        c(runif(12, 0, .2),  runif(11, 1, 2),   runif(9, 0.5, 1),  runif(8, 0, .2)),
                        c(runif(23, 0, .2),                     runif(9, 1, 2),    runif(8, 0.5, 1)),
                        c(runif(12, 0.5, 1), runif(11, 0, .2),  runif(9, 0, .2),   runif(8, 1, 2)))
  
  colnames(patterns_cor) <- labels
  patterns_cor
}

create_patterns_dist <- function (seed) {
  set.seed(seed)
  patterns_dist <- rbind(c(runif(12, 1, 2),  runif(11, 0, .2), runif(17, 0, .2)),
                         c(runif(12, 0, .2), runif(11, 1, 2),  runif(17, 0, .2)),
                         c(runif(23, 0, .2), runif(9, 1, 2),   runif(8, 0, .2)),
                         c(runif(23, 0, .2), runif(9, 0, .2),  runif(8, 1, 2)))
  colnames(patterns_dist) <- labels
  patterns_dist
}

distinct <- create_patterns_dist(1)
overlap <- create_patterns_cor(1)

plot_dist <- distinct %>% as_tibble() %>% 
  mutate(Pattern = 1:4) %>% dplyr::select(Pattern, everything()) %>% 
  gather(key = Chemicals, value = Loadings, -Pattern) %>% 
  mutate(Chemicals = fct_inorder(Chemicals)) %>% 
  mutate(Group = ifelse(str_detect(Chemicals, "BDE"), "PBDEs",
                        ifelse(str_detect(Chemicals, "^M"), "Phthalates",
                               ifelse(str_detect(Chemicals, "pcb"), "PCBs",
                                      "Phenols")))) %>% 
  mutate(Pattern = paste0("Pattern ", Pattern))

plot_over <- overlap %>% as_tibble() %>% 
  mutate(Pattern = 1:4) %>% dplyr::select(Pattern, everything()) %>% 
  gather(key = Chemicals, value = Loadings, -Pattern) %>% 
  mutate(Chemicals = fct_inorder(Chemicals)) %>% 
  mutate(Group = ifelse(str_detect(Chemicals, "BDE"), "PBDEs",
                        ifelse(str_detect(Chemicals, "^M"), "Phthalates",
                               ifelse(str_detect(Chemicals, "pcb"), "PCBs",
                                      "Phenols")))) %>% 
  mutate(Pattern = paste0("Pattern ", Pattern))

# NEJM "DeepCerulean" = "#0072B5"

#pdf("Figures/isee_2020_patterns.pdf", width = 12)
plot_over %>% 
  ggplot(aes(x = Chemicals, y = Loadings)) + 
  geom_col(color = "#0072B5", width = 0.75) +
  facet_grid(.~Pattern) + theme_bw(base_size = 40) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") +
  labs(title = "Overlapping Simulated Patterns") +
  scale_fill_nejm()
#dev.off()

#pdf("Figures/isee_2020_patterns_dist.pdf", width = 12)
plot_dist %>% 
  ggplot(aes(x = Chemicals, y = Loadings)) + 
  geom_col(color = "#0072B5", width = 0.75) +
  facet_grid(.~Pattern) + theme_bw(base_size = 40) +
  theme(axis.text.x = element_blank(),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  labs(title = "Distinct Simulated Patterns") +
  scale_fill_nejm()
#dev.off()

################################
################################

### Patterns
edc_cc <- 
  mn_edc %>% 
  dplyr::select(grep("sg_adj", colnames(.))) %>% 
  drop_na() %>% 
  clean_names()
names(edc_cc) <- str_replace(names(edc_cc), "_sg_adj", "")
eh_cc <- readMat(here::here("./MATLAB/Output/isee_2020_eh3.mat"))[[1]]
labels <-  toupper(names(edc_cc))
colnames(eh_cc) <- labels
eh_sd <- apply(eh_cc, 2, function(x) x/sum(x))

applied <- eh_sd %>% 
  as_tibble() %>% 
  mutate(Pattern = 1:3) %>% 
  gather(key = Chemicals, value = Loadings, -Pattern) %>%
  mutate(Chemicals = fct_inorder(Chemicals)) %>% 
  mutate(Class = case_when(Chemicals == "TCS" | Chemicals == "BPA" ~ "Phenols", 
                           grepl("PB", Chemicals) ~ "Parabens",
                           grepl("_", Chemicals) ~ "Phenols",
                           grepl("^M", Chemicals) == TRUE ~ "Phthalates")) %>% 
  mutate(Pattern = paste0("Pattern ", Pattern))

#pdf("Figures/isee_2020_applied.pdf", width = 18)
applied %>% 
  ggplot(aes(x = Chemicals, y = Loadings)) + geom_col(aes(fill = Class)) +
  facet_wrap(.~Pattern) + theme_bw(base_size = 40) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        legend.position = c(0.08, 0.65), # c(1,0) bottom left, c(1,1) top-right.
        legend.background = element_rect(fill = "#ffffffaa", colour = NA),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20)) +
  scale_fill_nejm() +
  labs(title = "Identified Patterns of EDC Exposure during Pregnancy")
#dev.off()
