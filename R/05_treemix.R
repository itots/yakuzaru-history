library(tidyverse)
library(patchwork)
library(ggplotify)
source("plotting_funcs.R")

# Input paths
input_path <- list(
  m0 = "../stats/treemix/treemix.0",
  m1 = "../stats/treemix/treemix.1",
  m2 = "../stats/treemix/treemix.2",
  m3 = "../stats/treemix/treemix.3",
  m4 = "../stats/treemix/treemix.4",
  m5 = "../stats/treemix/treemix.5"
)

plots$tree <- list(
  m0 = as.ggplot(~plot_tree(input_path$m0,
                            cex = 0.7, scale = FALSE)),
  m1 = as.ggplot(~plot_tree(input_path$m1,
                            cex = 0.7, scale = FALSE)),
  m2 = as.ggplot(~plot_tree(input_path$m2,
                            cex = 0.7, scale = FALSE)),
  m3 = as.ggplot(~plot_tree(input_path$m3,
                            cex = 0.7, scale = FALSE)),
  m4 = as.ggplot(~plot_tree(input_path$m4,
                            cex = 0.7, scale = FALSE)),
  m5 = as.ggplot(~plot_tree(input_path$m5,
                            cex = 0.7, scale = FALSE))
)

plots$wrap <- wrap_plots(
  plots$tree,
  ncol = 2
) & 
  coord_cartesian(clip = "off")

ggsave(
  filename = "output/treemix.pdf",
  plot = plots$wrap,
  width = 300,
  height = 320,
  units = "mm"
  )