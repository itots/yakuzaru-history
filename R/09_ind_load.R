library(tidyverse)
library(cowplot)

# Input paths
input_path <- list(
  samples_genome = "../list/samples_genome.csv",
  population = "../list/population.csv",
  fai = "../refseq/GCF_003339765.1_Mmul_10_genomic.fna.fai",
  roh = "../stats/plink_roh/roh.hom",
  load_dir = "../stats/genetic_load/individual"
)

# Color palettes
color <- list(
  three = c("#1784bf", "#d1990d", "#07360c")
)

# Variables
config <- list(
  ci_level = 0.95,
  snp_cols = c(
    "r_mis_sites", "r_mis_homo", "r_mis_alleles",
    "r_lof_sites", "r_lof_homo", "r_lof_alleles",
    "r_syn_sites", "r_syn_homo", "r_syn_alleles"
  ),
  indel_cols = c("r_lof_sites", "r_lof_homo", "r_lof_alleles")
)

# Functions
fnc <- list(
  
  predict_ci = function(mod, grid, level) {
    as_tibble(predict(mod, newdata = grid, interval = "confidence", level = level))
  },
  
  make_load_data = function(summary, vtype, cols) {
    summary %>%
      filter(variant_type == vtype) %>%
      pivot_longer(cols = all_of(cols)) %>%
      separate(col = name, into = c("r", "effect", "count_type")) %>%
      mutate(count_type = factor(count_type, levels = c("sites", "homo", "alleles")))
  },
  
  make_load_pred = function(load_data, ci_level) {
    load_data %>%
      group_by(effect, count_type) %>%
      nest() %>%
      mutate(
        lin = map(data, ~ lm(value ~ froh_total, data = .x)),
        quad = map(data, ~ lm(value ~ froh_total + I(froh_total^2), data = .x)),
        lin_bic = map_dbl(lin, BIC),
        quad_bic = map_dbl(quad, BIC),
        best = if_else(quad_bic < lin_bic, quad, lin),
        summary_best = map(best, summary),
        grid = map(data, ~ tibble(
          froh_total = seq(min(.x$froh_total, na.rm = TRUE),
            max(.x$froh_total, na.rm = TRUE),
            length.out = 200
          )
        )),
        pred_df = map2(best, grid, ~ bind_cols(.y, fnc$predict_ci(.x, .y, ci_level)))
      ) %>%
      select(effect, count_type, pred_df, best, summary_best) %>%
      unnest(pred_df) %>%
      ungroup()
  },
  
  plot_roh_load = function(load, pred_unnest, color, ncol = 3, legend = TRUE) {
    p <- load %>%
      ggplot(
        aes(
          x = froh_total,
          y = value,
          shape = population,
          color = cluster
        )
      ) +
      geom_ribbon(
        data = pred_unnest,
        mapping = aes(x = froh_total, ymin = lwr, ymax = upr),
        inherit.aes = FALSE,
        fill = "grey40", alpha = 0.2
      ) +
      geom_line(
        data = pred_unnest,
        mapping = aes(x = froh_total, y = fit),
        inherit.aes = FALSE,
        color = "black", linewidth = 0.3
      ) +
      geom_point() +
      scale_shape_manual(values = 0:10) +
      scale_color_manual(values = color$three) +
      #facet_wrap(
      #  count_type ~ effect,
      #  scale = "free",
      #  ncol = ncol
      #) +
      facet_wrap(
        ~ interaction(effect, count_type, sep = " - "),
       scale = "free",
       ncol = ncol
       ) +
      theme_minimal() +
      labs(
        x = "FROH",
        y = "Count"
      )
    if (!legend) {
      p <- p + theme(legend.position = "none")
    } else {
      p <- p + theme(
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 7),
        legend.key.spacing.y = unit(0, "pt")
      )
    }
    p
  }
)

# Load data
tbs <- list(
  info = read_csv(input_path$samples_genome),
  pop = read_csv(input_path$population) %>%
    mutate(population = fct_inorder(population)),
  fai = read_tsv(input_path$fai,
                 col_names = c("chr", "length", "offset", "linebases", "linewidth"), 
                 col_types = cols(
                   length = col_double()
                 )),
  roh = read_table(input_path$roh,
                   col_types = cols(
                     IID = col_character()
                   ))
)

# Calculate total autosome length
config$autosome_length <- tbs$fai %>%
  filter(chr %in% sprintf("NC_0417%02d.1", 54:73)) %>%
  summarize(total_length = sum(as.numeric(length))) %>%
  pull(total_length) / 1000

# Load genetic load data
tbs$load <- bind_rows(
  map(
    list.files(
      input_path$load_dir,
      pattern = "\\.txt",
      recursive = TRUE,
      full.names = TRUE
    ),
    ~ {
      path_parts <- strsplit(.x, .Platform$file.sep)[[1]]
      n <- length(path_parts)
      variant_type <- ifelse(n >= 2, path_parts[n-1], NA)
      effect <- ifelse(n >= 3, path_parts[n-2], NA)
      read_table(
        .x,
        col_names = c("sample_id", "stats", "value")
      ) %>%
        mutate(
          effect = effect,
          variant_type = variant_type,
          chrom = basename(.x)
        )
    }
  )
)

# Calculate mean load per type
tbs$load_mean <- tbs$load %>%
  pivot_wider(
    names_from = c("effect", "stats"),
    values_from = value,
    values_fill = NA
  ) %>%
  group_by(
    variant_type, sample_id
  ) %>%
  dplyr::summarise(
    dplyr::across(where(is.numeric) & !dplyr::any_of("sample_id"), 
                  \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  group_by(variant_type) %>%
  dplyr::summarise(
    dplyr::across(where(is.numeric), 
                  \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::rename_with(~ paste0(.x, "_mean"), where(is.numeric))

# Summarize load data and merge with roh data
tbs$summary <- tbs$info %>% 
  filter(
    !is.na(cluster)
  ) %>%
  left_join(
    tbs$load %>%
      pivot_wider(
        names_from = c("effect", "stats"),
        values_from = value,
        values_fill = NA
      ) %>%
      group_by(
        variant_type, sample_id
      ) %>%
      summarize(
        across(4:ncol(.)-2, \(x) sum(x, na.rm = TRUE)),
        .groups = "drop"
        ) %>%
      left_join(
        tbs$load_mean,
        by = "variant_type"
      ) %>%
      mutate(
        r_lof_sites = lof_sites / lof_called_sites * lof_called_sites_mean,
        r_lof_homo = lof_homozygotes / lof_called_sites * lof_called_sites_mean,
        r_lof_alleles = lof_alleles / lof_called_sites * lof_called_sites_mean,
        r_mis_sites = missense_sites / missense_called_sites * missense_called_sites_mean,
        r_mis_homo = missense_homozygotes / missense_called_sites * missense_called_sites_mean,
        r_mis_alleles = missense_alleles / missense_called_sites * missense_called_sites_mean,
        r_syn_sites = synonymous_sites / synonymous_called_sites * synonymous_called_sites_mean,
        r_syn_homo = synonymous_homozygotes / synonymous_called_sites * synonymous_called_sites_mean,
        r_syn_alleles = synonymous_alleles/ synonymous_called_sites * synonymous_called_sites_mean
      ),
    by = c("seq_id" = "sample_id")
  ) %>%
  left_join(
    tbs$roh %>%
      filter(
        CHR <= 20
      ) %>%
      group_by(
        IID
      ) %>%
      summarise(
        total_length = sum(KB),
        long_length = sum(KB[KB >= 2000]),
        froh_total = total_length / config$autosome_length,
        froh_long = long_length / config$autosome_length
      ),
    by = c("seq_id" = "IID")
  ) %>%
  filter(
    !is.na(total_length)
  ) %>%
  mutate(
    population = factor(population, levels = tbs$pop$population)
  )

# Calculate missense / synonymous ratio
plots$NS <- tbs$summary %>%
  filter(variant_type == "SNP") %>%
  mutate(mis_syn = missense_heterozygotes / synonymous_heterozygotes) %>%
  ggplot(aes(x = population, y = mis_syn, color = cluster)) +
  geom_boxplot(linewidth = 0.3) +
  scale_color_manual(values = color$three) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  ) +
  labs(
    y = "Missense / Synonymous"
  )

# Summarize and plot the relationship between Load and ROH
tbs$snp_load <- fnc$make_load_data(tbs$summary, "SNP", config$snp_cols)
tbs$snp_load_pred <- fnc$make_load_pred(tbs$snp_load, config$ci_level)

tbs$indel_load <- fnc$make_load_data(tbs$summary, "INDEL", config$indel_cols)
tbs$indel_load_pred <- fnc$make_load_pred(tbs$indel_load, config$ci_level)

plots <- list(
  snp_roh_load = fnc$plot_roh_load(tbs$snp_load, tbs$snp_load_pred, color, ncol = 3, legend = TRUE),
  indel_roh_load = fnc$plot_roh_load(tbs$indel_load, tbs$indel_load_pred, color, ncol = 1, legend = FALSE)
)

plots$snv_indel_ns <- 
  plot_grid(
    plot_grid(
      plots$snp_roh_load + theme(legend.position = "none"),
      NULL,
      plots$indel_roh_load,
      ncol = 3,
      rel_widths = c(1, 0.05,0.34)
  ),
  NULL,
  plot_grid(
    plots$NS + theme(legend.position = "none"),
    get_legend(plots$snp_roh_load + theme(legend.box = "horizontal")),
    ncol = 2,
    rel_widths = c(1,0.7)
    ),
  ncol = 1,
  rel_heights = c(1,0.05,0.5)
)

ggsave(
  "output/load.pdf", 
  plots$snv_indel_ns,
  width = 200,
  height = 200,
  units = "mm"
  )

