library(tidyverse)
library(admixtools)
library(furrr)
library(purrr)
library(patchwork)

# Input paths
input_path <- list(
  f2 = list(
    A = "../stats/admixtools/f2/A",
    X = "../stats/admixtools/f2/X"
  ),
  best_edges = list(
    k0 = "../stats/admixtools/runs/best_edges_k0.rds",
    k1 = "../stats/admixtools/runs/best_edges_k1.rds",
    k2 = "../stats/admixtools/runs/best_edges_k2.rds",
    k3 = "../stats/admixtools/runs/best_edges_k3.rds",
    k4 = "../stats/admixtools/runs/best_edges_k4.rds",
    k5 = "../stats/admixtools/runs/best_edges_k5.rds"
  )
)

# Define values
config <- list(
  pops_fus = c(
    "fus_Yakushima",
    "fus_Kojima", "fus_Takasakiyama", "fus_Nakatosa", "fus_Okazaki", "fus_Minoh", "fus_Tsubaki",
    "fus_Hagachizaki", "fus_Jigokudani", "fus_Yamagata", "fus_Shimokita"
  ),
  pops_fus_sub = c(
    "fus_Yakushima",
    "fus_Kojima", "fus_Nakatosa", "fus_Minoh", "fus_Tsubaki",
    "fus_Hagachizaki", "fus_Jigokudani", "fus_Yamagata", "fus_Shimokita"
  ), 
  pops_out = c(
    "fas_Indonesia", "cyc_Taiwan", "mul_India", 
    "mul_China", "mul_Sichuan", "mul_Hainan", "mul_Hubei", "mul_Yunnan" 
  ),
  pops_chm = c("cyc_Taiwan", "mul_China", "mul_Sichuan", "mul_Hainan", "mul_Hubei", "mul_Yunnan", "mul_India")
)

# f2 blocks
f2 <- list(
  A = read_f2(input_path$f2$A),
  X = read_f2(input_path$f2$X)
)

# Import data
results <- list()
results$best_edges <- list(
  k0 = readRDS(input_path$best_edges$k0),
  k1 = readRDS(input_path$best_edges$k1),
  k2 = readRDS(input_path$best_edges$k2),
  k3 = readRDS(input_path$best_edges$k3),
  k4 = readRDS(input_path$best_edges$k4),
  k5 = readRDS(input_path$best_edges$k5)
)

# Calculate fit statistics for each graph 
results$fit$A <- list(
  k0 = qpgraph(f2$A, results$best_edges$k0, return_fstats = TRUE),
  k1 = qpgraph(f2$A, results$best_edges$k1, return_fstats = TRUE),
  k2 = qpgraph(f2$A, results$best_edges$k2, return_fstats = TRUE),
  k3 = qpgraph(f2$A, results$best_edges$k3, return_fstats = TRUE),
  k4 = qpgraph(f2$A, results$best_edges$k4, return_fstats = TRUE),
  k5 = qpgraph(f2$A, results$best_edges$k5, return_fstats = TRUE)
)

results$fit$X <- list(
  k0 = qpgraph(f2$X, results$best_edges$k0, return_fstats = TRUE),
  k1 = qpgraph(f2$X, results$best_edges$k1, return_fstats = TRUE),
  k2 = qpgraph(f2$X, results$best_edges$k2, return_fstats = TRUE),
  k3 = qpgraph(f2$X, results$best_edges$k3, return_fstats = TRUE),
  k4 = qpgraph(f2$X, results$best_edges$k4, return_fstats = TRUE),
  k5 = qpgraph(f2$X, results$best_edges$k5, return_fstats = TRUE)
)


# qpWave analysis
results$qpwave <- qpwave(
  data = f2$A, 
  left = config$pops_fus_sub, 
  right = config$pops_out
  )

write_csv(
  results$qpwave$rankdrop,
  "output/qpwave_rankdrop.csv"
)

# f3 statistics
results$f3_IntSp$A <- f3(
  data = f2$A, 
  pop1 = "fas_Indonesia", 
  pop2 = config$pops_chm, 
  pop3 = config$pops_fus
  )

results$f3_Yak$A <- f3(
  data = f2$A, 
  pop1 = "fas_Indonesia", 
  pop2 = config$pops_fus[config$pops_fus != "fus_Yakushima"], 
  pop3 = "fus_Yakushima"
)


# f4 statistics
results$f4_Yak$A <- f4(
  data = f2$A,
  pop1 = "fas_Indonesia",
  pop2 = "fus_Yakushima",
  pop3 = config$pops_fus[config$pops_fus != "fus_Yakushima" & config$pops_fus != "fus_Shimokita"],
  pop4 = "fus_Shimokita"
)

results$f4_Yak$X <- f4(
  data = f2$X,
  pop1 = "fas_Indonesia",
  pop2 = "fus_Yakushima",
  pop3 = config$pops_fus[config$pops_fus != "fus_Yakushima" & config$pops_fus != "fus_Shimokita"],
  pop4 = "fus_Shimokita"
)

# Worst residual
results$z_score_summary <- imap_dfr(results$fit, function(chr_list, chr_name) {
  tibble(
    chr            = chr_name,                             
    k              = names(chr_list),                      
    worst_residual = map_dbl(chr_list, ~ .x$worst_residual),
    score          = map_dbl(chr_list, ~ .x$score)
  )
})

write_csv(
  results$z_score_summary,
  "output/z_score_summary.csv"
)

# Plot f3 statistics
plots <- list(
  f3_IntSp = results$f3_IntSp$A %>%
  mutate(
    pop2 = factor(pop2, levels = config$pops_chm),
    pop3 = factor(pop3, levels = config$pops_fus)
  ) %>%
  ggplot(aes(x = pop2, y = pop3, fill = est)) +
  geom_tile(color = "white") +          
  geom_text(
    aes(
      label = sprintf("%.4f", est), 
      color = est < 0.053 
    ), 
    size = 3, 
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "white")) +
  scale_fill_viridis_c(option = "plasma") + 
  theme_minimal() +                      
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) +
  labs(
    x = "Population 2",
    y = "Population 3",
    fill = "f3"
  )
)

plots$f3_Yak$A <- 
results$f3_Yak$A %>%
  mutate(
    pop2 = factor(pop2, levels = config$pops_fus[config$pops_fus != "fus_Yakushima"])
  ) %>%
  ggplot(aes(x = pop2, y = est)) +
  geom_point(stat = "identity") + 
  geom_errorbar(aes(ymin = est - se, ymax = est + se), width = 0.1) +
  theme_minimal() +                      
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) +
  labs(
    x = "Population 2",
    y = "f3",
    fill = "f3"
  )

plots$f4_Yak <- 
bind_rows(
  results$f4_Yak$A %>% mutate(chr = "A"),
  results$f4_Yak$X %>% mutate(chr = "X")
) %>%
  mutate(
    pop3 = factor(pop3, levels = config$pops_fus[config$pops_fus != "fus_Yakushima" & config$pops_fus != "fus_Shimokita"])
  ) %>%
  ggplot(aes(x = pop3, y = est, color = chr)) +
  geom_point(stat = "identity", position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = est - se, ymax = est + se), width = 0.1, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("A" = "black", "X" = "gray60"), labels = c("Autosomes", "X chromosome")) +
  theme_minimal() +                      
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) +
  labs(
    x = "Population 3",
    y = "f4",
    color = "Chromosome"
  )


plots$graph <- 
results$fit %>%
  modify_depth(2, function(x) {
    plot_graph(x$edges, textsize = 2) 
  })

plots$graph_wrap <-  wrap_plots(plots$graph$A, ncol = 2)


ggsave(
  filename = "output/wrap_graph.pdf",
  plot = plots$graph_wrap,
  width = 260,
  height = 300,
  units = "mm"
)

ggsave(
  filename = "output/f3_IntSp.pdf",
  plot = plots$f3_IntSp,
  width = 180,
  height = 140,
  units = "mm"
)

ggsave(
  filename = "output/f3_Yak.pdf",
  plot = plots$f3_Yak$A,
  width = 180,
  height = 140,
  units = "mm"
)

ggsave(
  filename = "output/f4_Yak.pdf",
  plot = plots$f4_Yak,
  width = 180,
  height = 140,
  units = "mm"
)