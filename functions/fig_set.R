## Figure Settings
library(ggsci)
library(tidyverse)

# ggsci_db$"nejm"$"default" <- c(
#   "TallPoppy" = "#BC3C29", "DeepCerulean" = "#0072B5",
#   "Zest" = "#E18727", "Eucalyptus" = "#20854E",
#   "WildBlueYonder" = "#7876B1", "Gothic" = "#6F99AD",
#   "Salomie" = "#FFDC91", "FrenchRose" = "#EE4C97"
# )

theme_set(theme_bw(base_size = 26) + 
            theme(strip.background = element_rect(fill="white"),
                  axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
                  axis.title.x = element_blank(),
                  legend.direction = "horizontal",
                  legend.title = element_blank(),
                  legend.position = "none",
                  plot.margin = unit(c(5.5, 30, 50, 5.5), "points")) # top, then right, bottom and left. 
)

ggsci <- pal_nejm()(8)

options(
  ggplot2.discrete.color = ggsci,
  ggplot2.discrete.fill = ggsci)