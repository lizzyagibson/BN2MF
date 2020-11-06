## Figure Settings
library(ggsci)
library(tidyverse)

theme_set(theme_bw(base_size = 26) + 
            theme(strip.background = element_rect(fill="white"),
                  # axis.text.x = element_blank(), 
                  axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
                  axis.title.x = element_blank(),
                  legend.direction = "horizontal",
                  legend.title = element_blank(),
                  legend.position = "none",
                  # legend.position = c(0.5, -0.05), # c(0,0) bottom left, c(1,1) top-right.
                  axis.title.y = element_text(size = 40),
                  plot.margin = unit(c(5.5, 30, 40, 5.5), "points")) # top, then right, bottom and left. 
)

ggsci <- pal_nejm()(8)

options(
  ggplot2.discrete.colour = ggsci,
  ggplot2.discrete.fill = ggsci)