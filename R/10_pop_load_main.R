library(tidyverse)
library(purrr)
library(cowplot)

source("10_pop_load_function.R")

# Input paths
input_paths <- list(
  chromosome_sizes = "../list/chromosome.txt",
  samples_genome = "../list/samples_genome.csv",
  population = "../list/population.csv",
  load_dir = "../stats/genetic_load/population"
)

# Settings
config <- list(
  effects = c("lof","missense","synonymous","intergenic")
)

# Load data
tbs <- list(
  info = read_csv(input_paths$samples_genome),
  pop = read_csv(input_paths$population) %>%
    mutate(population = fct_inorder(population)),
  chr_lens = fnc$read_genome_sizes(input_paths$chromosome_sizes)
)

tbs$load <- bind_rows(
  map(
    list.files(
      input_paths$load_dir,
      pattern = "\\.txt",
      recursive = TRUE,
      full.names = TRUE
    ),
    ~ {
      path_parts <- strsplit(.x, .Platform$file.sep)[[1]]
      n <- length(path_parts)
      type <- ifelse(n >= 2, path_parts[n-1], NA)
      if (!identical(type, "SNP")) return(NULL)

      effect <- ifelse(n >= 3, path_parts[n-2], NA)
      population <- ifelse(n >= 4, path_parts[n-3], NA)
      read_table(
        .x,
        col_names = c("chrom_no", "position", "ac", "an"),
        col_types = cols(
          chrom_no = col_character(),
          position = col_integer(),
          ac = col_integer(),
          an = col_integer()
        )
      ) %>%
        mutate(
          population = population,
          effect = effect,
          type = type,
          chrom = basename(.x),
          af = ac / an
        )
    }
  )
)


# Yakushima
stats <- list()
stats$f <- setNames(
  map(config$effects, ~ fnc$make_f(tbs$load, .x)),
  config$effects
)
stats$res_lof <- fnc$lobo_rxy(
  stats$f$lof, stats$f$intergenic, B = 100, chr_lens = tbs$chr_lens
)
stats$res_missense <- fnc$lobo_rxy(
  stats$f$missense, stats$f$intergenic, B = 100, chr_lens = tbs$chr_lens
)
stats$res_synonymous <- fnc$lobo_rxy(
  stats$f$synonymous, stats$f$intergenic, B = 100, chr_lens = tbs$chr_lens
)

plots <- list()
plots$rxy <- 
bind_rows(
  stats$res_lof$summary        %>% mutate(effect = "lof"),
  stats$res_missense$summary   %>% mutate(effect = "missense"),
  stats$res_synonymous$summary %>% mutate(effect = "synonymous")
) %>%
  mutate(
    pair   = factor(pair, levels = c("ew","ey","wy")),
    effect = factor(effect, levels = config$effects),
    delta  = Rxy_full - 1,
    ci_lo  = ci_lower - 1,
    ci_hi  = ci_upper - 1
  ) %>% ggplot(aes(x = pair, y = delta, fill = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_col(position = position_dodge(width = 0.7), width = 0.7, linewidth = 0.1, color = "black",key_glyph = "polygon") +
  geom_errorbar(
    aes(ymin = ci_lo, ymax = ci_hi, group = effect),
    width = 0.12,
    position = position_dodge(width = 0.7),
    linewidth = 0.3,
    color = "grey50",
    alpha = 1,
    show.legend = FALSE
  ) +
  scale_x_discrete(labels = c(
    ew = "East/\nWest",
    ey = "East/\nYakushima",
    wy = "West/\nYakushima"
  )) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = function(b) scales::number(b + 1, accuracy = 0.01)
  ) +
  scale_fill_manual(
    values = c("grey0", "grey60", "grey90", "grey100"),
    limits = c("lof","missense","synonymous", "intergenic")
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        colour = "black", 
        shape = 22, 
        size = 5, 
        linetype = 1 
      ) 
    )
  ) +
  labs(x = "Cluster pair", y = "Rxy", fill = "Effect") +
  theme_minimal()

# SFS
tbs$sfs_summary <- fnc$project_and_summarise_sfs(
  tbs$load %>% filter(population != "Yakushima_wild", type == "SNP"), 
  k = 14, 
  reps = 100, 
  seed = 1234
)

plots$sfs <- 
tbs$sfs_summary %>%
  mutate(effect = factor(effect, levels = config$effects)) %>%
  ggplot(aes(x = dac, y = mean, fill = effect)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8, linewidth = 0.1, color = "black") +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                position = position_dodge(width = 0.8), width = 0.3, color = "grey50") +
  facet_grid(rows = vars(population), scales = "fixed") +
  scale_x_continuous(breaks = 1:13) +
  scale_fill_manual(
    values = c("grey0", "grey60", "grey90", "grey100")
  ) +
  labs(x = "Derived allele count",
       y = "Proportion of sites",
       fill = "Effect") +
  theme_minimal()

plots$rxy_sfs <- 
plot_grid(
  plot_grid(
    plots$rxy + theme(legend.position = "none"),
    get_legend(plots$sfs + theme(legend.title = element_blank())),
    ncol = 1,
    rel_heights = c(1,0.5)
  ),
  plots$sfs + theme(legend.position = "none"),
  ncol = 2,
  rel_widths = c(1,1.7)
)


ggsave(
  "output/rxy_sfs.pdf",
  plots$rxy_sfs,
  width = 260,
  height = 160,
  units = "mm"
)

save.image()
