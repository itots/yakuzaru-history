library(tidyverse)
library(cowplot)

# Input paths
input_path <- list(
  pop_dir = "/Volumes/research/yakuzaru_history/stats//msmc/pop",
  ccc_dir = "/Volumes/research/yakuzaru_history/stats/msmc/combined_cc",
  pop = "../../list/population.csv",
  info = "../../list/samples_genome.csv",
  smcpp = "/Volumes/research/yakuzaru_history/stats/smc++/estimate/plot.csv"
)

# Define variables
config <- list(
  mu = 7.7e-9,
  gen = 10.4
)

# Color setting
color <- list(
  three = c("#1784bf", "#d1990d", "#07360c")
)

# Load data
tbs <- list(
  popinfo = read_csv(input_path$pop) %>%
    mutate(population = fct_inorder(population))
)

tbs$info <- read_csv(input_path$info) %>%
  mutate(
    population = factor(population, levels = tbs$popinfo$population)
  )

tbs$pops <- bind_rows(
  map(
    list.files(
      path = input_path$pop_dir,
      pattern = "final.txt",
      full.names = TRUE
    ),
    function(file){
      fname = basename(file)
      population = str_extract(fname, "^[^_]+")
      sample_size = str_extract(fname, "(?<=_)[0-9]+(?=samples)")
      
      read_tsv(file, comment = "#") %>%
        mutate(
          years = left_time_boundary / config$mu * config$gen,
          ne = 1 / (2 * lambda * config$mu),
          population = population,
          sample_size = sample_size
        )
    }
  )
) %>%
  left_join(
    tibble(
      population = c(
        "Yakushima", 
        "Hagachizaki", "Jigokudani", "Shimokita", "Yamagata", 
        "Kojima", "Minoh", "Nakatosa", "Tsubaki"
        ),
      cluster = c("Yakushima", rep("East", 4), rep("West", 4))
    )
  ) %>%
  mutate(
    population = factor(population, levels = tbs$popinfo$population)
  )

tbs$divs <- bind_rows(
  map(
    list.files(
      path = input_path$ccc_dir,
      pattern = ".txt",
      full.names = TRUE
    ),
    function(file){
      fname = basename(file)
      pop1 = str_split(fname, "_")[[1]][1]
      pop2 = str_split(fname, "_")[[1]][2] %>% str_remove("\\.txt")
      
      read_tsv(
        file, comment = "#"
      ) %>%
        mutate(
          pop1 = pop1,
          pop2 = pop2,
          years = left_time_boundary / config$mu * config$gen,
          rccr = (2 * lambda_01) / (lambda_00 + lambda_11)
        )
    }
  )
) %>%
  mutate(
    pair = case_when(
      pop1 == "Yakushima" & pop2 %in% c("Hagachizaki", "Jigokudani", "Shimokita", "Yamagata") ~ "East-Yakushima",
      pop1 == "Yakushima" & pop2 %in% c("Kojima", "Minoh", "Nakatosa", "Tsubaki") ~ "West-Yakushima",
      TRUE ~ "East-West"
    )
  ) %>%
  mutate(
    pair = factor(pair, levels = c("East-Yakushima", "West-Yakushima", "East-West"))
  )

tbs$smcpp <- read_csv(
  input_path$smcpp
)


# plotting
plots <- list(
  pops = tbs$pops %>% 
    filter(
      time_index >= 2,
      time_index <= 31,
      sample_size >= 2
    ) %>%
    ggplot(
      aes(
        x = years,
        y = ne,
        group = interaction(cluster, sample_size, population),
        color = cluster,
        linetype = sample_size,
        shape =  population
      )
    ) +
    geom_step(linewidth = 0.2) +
    geom_point(stroke = 0.2) +
    scale_x_log10(breaks = c(1e3, 1e4, 1e5, 1e6),
                  labels = c(1e3, 1e4, 1e5, 1e6),
                  minor_breaks = unlist(lapply(1:6, function(i) seq(10^i, 10^(i+1)-10^i, by=10^i))),
                  limits = c(2e2, 2e6)) +
    scale_y_log10(minor_breaks = unlist(lapply(1:6, function(i) seq(10^i, 10^(i+1)-10^i, by=10^i))),
                  limits = c(1.5e3, 2.5e5)) +
    scale_linetype_manual(
      values = c("2" = "solid", "5" = "dashed"),
      labels = c("2 samples", "5 samples")
    ) +
    scale_shape_manual(values = c(0:4,6:11))  +
    scale_color_manual(
      values = color$three,
      labels = c("East", "West", "Yakushima")
    ) +
    labs(
      y = "Ne"
    ) +
    theme_minimal() +
    theme(
      legend.key.size = unit(0.4, "cm"),   
      legend.text = element_text(size = 7),
      legend.key.spacing.y = unit(0, "pt"),
      legend.box = "horizontal",
      axis.title.x = element_blank()
    )
)


plots$divs <- tbs$divs %>%
  filter(
    time_index >= 2,
    time_index <= 31,
  ) %>%
  mutate(
    population_pair = paste(pop1, pop2, sep = "-"),
    cluster_pair = pair
  ) %>%
  filter(years > 0) %>%
  filter(time_index <= 25) %>%
  ggplot(
    aes(
      x = years,
      y = rccr,
      group = interaction(pop1, pop2, cluster_pair),
      color = cluster_pair,
      shape = population_pair
    )
  ) +
  geom_step(linewidth = 0.2) +
  geom_point(stroke = 0.2) +
  scale_x_log10(
                minor_breaks = unlist(lapply(1:6, function(i) seq(10^i, 10^(i+1)-10^i, by=10^i))),
                breaks = c(1e3, 1e4, 1e5, 1e6),
                labels = c(1e3, 1e4, 1e5, 1e6),
                limits = c(2e2, 2e6)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_shape_manual(values = 0:12) +
  scale_color_manual(values = c(
    "East-Yakushima" = color$three[1],
    "West-Yakushima" = color$three[2],
    "East-West" = "grey40")) +
  labs(
    x = "Years before present (g = 10.4, μ = 7.7e-9)",
    y = "rCCR"
  )  +
  theme_minimal() +
  theme(
    legend.key.size = unit(0.4, "cm"),   
    legend.text = element_text(size = 7),
    legend.key.spacing.y = unit(0, "pt"),
    legend.box = "horizontal"
  )

plots$pops_divs <- plot_grid(
  plots$pops, plots$divs, ncol = 1, align = "hv"
)

plots$smcpp <- tbs$smcpp %>%
  slice(1:(n() - 1)) %>%
  filter(x > 10) %>%
  ggplot(
    aes(
      x = x,
      y = y
    )
  ) +
  geom_step(linewidth = 0.5, color = color$three[3]) +
  scale_x_log10(breaks = c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5),
                labels = c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5),
                minor_breaks = unlist(lapply(1:6, function(i) seq(10^i, 10^(i+1)-10^i, by=10^i))),
                limits = c(100, 2e6)) +
  scale_y_log10(minor_breaks = unlist(lapply(1:6, function(i) seq(10^i, 10^(i+1)-10^i, by=10^i))),
                limits = c(1e2, 2e4)) +
  labs(
    x = "Years before present (g = 10.4, μ = 7.7e-9)",
    y = "Ne"
  )  +
  theme_minimal()

ggsave(
  "output/pops_divs.pdf", 
  plots$pops_divs,
  width = 260,
  height = 130,
  units = "mm",
  device = cairo_pdf)

ggsave(
  "output/smcpp.pdf", 
  plots$smcpp,
  width = 160,
  height = 80,
  units = "mm",
  device = cairo_pdf
)


# Divergence time and Harmonic mean Ne
tbs$divs %>% 
  mutate(
    generation = years / config$gen 
  ) %>%
  filter(
    time_index >= 2,
    time_index <= 31,
    pop1 == "Yakushima",
    rccr <= 0.5,
  ) %>% 
  arrange(-generation) %>%
  slice(1)

tbs$pops %>% 
  mutate(
    generation = years / config$gen 
  ) %>%
  filter(
    time_index >= 2,
    time_index <= 31,
    sample_size == 2,
    population == "Yakushima",
    generation <= 22437
  ) %>% 
  mutate(
    time_interval = right_time_boundary - left_time_boundary
  ) %>%
  dplyr::summarise(
    harmonic_ne = sum(time_interval) / sum(time_interval / ne),
    arithmetic_ne = sum(time_interval * ne) / sum(time_interval), 
    ne = mean(ne)
  )