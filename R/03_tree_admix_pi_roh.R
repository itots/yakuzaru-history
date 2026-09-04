library(tidyverse)
library(ape)
library(tidygenomes)
library(ggtree)
library(bnpsd)
library(paletteer)

# Input paths
input_path <- list(
  pop = "../list/population.csv",
  info = "../list/samples_genome.csv",
  pixy = "../stats/pixy/pi_ind",
  roh = "../stats/plink_roh/roh.hom",
  seqid = "../stats/data/plink/out/jm_ld_pruned.fam",
  iqtree = "../stats/iqtree/main/thined.min4.variants.treefile",
  admixture_cv = "../stats/admixture/cross_validation.txt",
  admixture_q = "../stats/admixture/jm_ld_pruned_autosome.3.Q",
  fai = "../stats/refseq/GCF_003339765.1_Mmul_10_genomic.fna.fai"
)

# Color palettes
color <- list(
  three = c("#1784bf", "#d1990d", "#07360c")
)

# Load data
tbs <- list(
  pop = read_csv(input_path$pop) %>%
    mutate(population = fct_inorder(population))
)

tbs$info <- read_csv(input_path$info) %>%
  mutate(
    population = factor(population, levels = tbs$pop$population)
  ) %>%
  arrange(population)

tbs$fai <- read_tsv(input_path$fai,
                    col_names = c("chr", "length", "offset", "linebases", "linewidth"), 
                    col_types = cols(
                      length = col_double()
                    ))

config <- list(
  autosome_length = tbs$fai %>%
    filter(chr %in% sprintf("NC_0417%02d.1", 54:73)) %>%
    summarize(total_length = sum(as.numeric(length))) %>%
    pull(total_length) /1000
)

# Pixy (heterozygosity)
tbs$pixy <- list.files(input_path$pixy, pattern = "*\\.txt", full.names = TRUE) %>%
  map_dfr(~ read_tsv(.x, col_types = cols()), .id = "file") %>%
  select(-file) %>%
  rename(
    seq_id = pop
  ) %>%
  left_join(
    tbs$info,
    by = "seq_id"
  ) %>%
  group_by(
    seq_id,
    population,
    species
  ) %>%
  summarise(
    genome_avg_pi = sum(count_diffs) / sum(count_comparisons),
    .groups = "drop"
  ) %>%
  select(
    seq_id, genome_avg_pi
  ) %>%
  left_join(
    tbs$info,
    by = "seq_id"
  )

# ROH 
tbs$roh <- read_table(input_path$roh,
                      col_types = cols(
                        IID = col_character(),
                        SNP1 = col_character(),
                        SNP2 = col_character()
                        )
                      ) %>%
  left_join(
    tbs$info,
    by = c("IID" = "seq_id")
  )

tbs$roh_summary <- tbs$roh %>%
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
  ) %>%
  left_join(
    tbs$info,
    by = c("IID" = "seq_id")
  )


# Load tree
tbs$seqid <- read.table(input_path$seqid) %>%
  as_tibble() %>%
  select(seq_id = V2)

tree <- list(
  iqtree = read.tree(input_path$iqtree)
)

tree$iqtree_rooted <- root(tree$iqtree, outgroup = c("SRR24738440", "SRR24738437"), resolve.root = TRUE)
tree$iqtree_rooted_jm <- drop.tip(tree$iqtree_rooted, tbs$info$seq_id[is.na(tbs$info$population)])
tree$iqtree_rooted_jm_rotated <- rotateConstr(tree$iqtree_rooted_jm, rev(tbs$info %>% right_join(tbs$seqid) %>% .$seq_id))

# Load admixture data
tbs$cv <- read.table(input_path$admixture_cv, header = FALSE, col.names = c("V1", "V2", "V3", "V4")) %>%
  filter(V1 == "CV")
tbs$cv$K <-gsub("[\\(\\)]", "", regmatches(tbs$cv$V3, gregexpr("\\(.*?\\)", tbs$cv$V3))) 
tbs$cv <- tbs$cv %>% 
  select(
    CV = V4,
    K
  )

tbs$admix <- 
bind_cols(
  tbs$seqid %>%
    left_join(
      tbs$info,
      by = "seq_id"
    ),
  read.table(input_path$admixture_q) %>%
    as_tibble() %>%
    select(K1 = V1, K2 = V2, K3 = V3)
) %>%
  pivot_longer(
    cols = starts_with("K"),
    names_to = "K",
    values_to = "ancestry"
  ) 

tbs$admix_info <- tbs$info %>% 
right_join(tbs$seqid) %>%
mutate(
  population = fct_inorder(population),
  cluster = fct_inorder(cluster)
) %>%
 droplevels()
 
# Plotting
plots <- list(
  tree = ggtree(tree$iqtree_rooted_jm_rotated, linewidth = 0.2) %>%
    rotate(43) %<+% 
    tbs$admix_info +
    scale_y_reverse() +
    geom_tiplab(aes(label=""), align = TRUE, linetype = "dashed", color = "grey50") +
    geom_tippoint(aes(x = max(x, na.rm = TRUE), shape = population, color = cluster), size = 2, stroke = 0.4) +
    scale_shape_manual(values = 0:11) +
    scale_colour_manual(values = color$three) +
    theme_tree2() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
)

plots$tree_admix <- facet_plot(
  plots$tree,
    panel = 'Admixture',
    data = tbs$admix,
    mapping = aes(x = ancestry, fill = factor(K, levels = c("K2", "K3", "K1"))),
    geom = geom_bar,
    orientation = 'y', width = 0.8, stat='identity'
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  scale_fill_manual(values = color$three)

plots$tree_admix_pi <- facet_plot(
  plots$tree_admix,
  panel = "Heterozygosity",
  data = tbs$pixy,
  mapping = aes(x = genome_avg_pi),
  geom = geom_bar,
  orientation = 'y', width = 0.8, stat='identity'
) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

plots$tree_admix_pi_roh <- facet_plot(
  plots$tree_admix_pi,
  panel = "FROH",
  data = tbs$roh_summary %>%
    mutate(
      froh_short = froh_total - froh_long
    ) %>%
    pivot_longer(
      cols = c("froh_short", "froh_long"),
      names_to = "type",
      values_to = "froh"
    ) %>%
    mutate(
      type = factor(type, levels = c("froh_short", "froh_long"))
    ),
  mapping = aes(x = froh, fill = type),
  geom = geom_bar,
  orientation = 'y', width = 0.8, stat='identity'
) +
  scale_fill_manual(values = c(color$three, "grey75", "grey30")) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_line(color = "grey70", linewidth = 0.4), 
    panel.grid.minor.x = element_line(color = "grey85", linewidth = 0.2),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_blank(),        
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) 

# Admixture CV
plots$admixture_cv <- tbs$cv %>%
  mutate(
    K = parse_number(K)
  ) %>%
  ggplot(aes(x = K, y = CV)) +
  geom_line(group = 1) +
  geom_point() +
  theme_minimal()


# Cumulative ROH
plots$cum_roh <- tbs$roh %>%
  filter(
    species == "fuscata"
  ) %>%
  mutate(
    time_coal = 100 / (2 * 0.433 * KB/1000)
  ) %>%
  mutate(
    population = factor(population, levels = tbs$pop$population)
  ) %>%
  arrange(cluster, IID, KB) %>%
  group_by(cluster, IID) %>%
  mutate(Cumulative_Proportion = cumsum(KB)/config$autosome_length) %>%
  ungroup() %>%
  ggplot(
    aes(
      x = KB,
      y = Cumulative_Proportion,
      color = cluster,
      group = IID
    )
  ) +
  geom_step(
    linewidth = 0.2
  ) +
  geom_point(
    data = . %>% 
      group_by(cluster, IID) %>% 
      filter(KB == max(KB)) %>% 
      ungroup(),
    aes(
      x = KB,
      y = Cumulative_Proportion,
      color = cluster,
      shape = population
    ),
  ) +
  scale_x_log10(
    name = "Length of ROH (kb)",
    minor_breaks = unlist(lapply(1:6, function(i) seq(10^i, 10^(i+1)-10^i, by=10^i))),
    
    sec.axis = sec_axis(
      transform = ~ (10.4 * 100000 / (2 * 0.433)) / ., 
      name = "Time of coalescence (years before present)",
      breaks = c(100, 200, 500, 200, 1000, 2000, 5000, 10000)
    )
  ) +
  scale_color_manual(values = color$three) +
  scale_shape_manual(values = 0:10) +
  theme_minimal() +
  labs(
    y = "Cumulative proportion\n of genome in ROH"
  ) +
  theme(
    legend.position = "none"
  ) 


ggsave("output/genpop.pdf", plots$tree_admix_pi_roh,
       width = 200, height = 120, units = "mm", device = cairo_pdf)

ggsave("output/admixture_cv.pdf", plots$admixture_cv,
       width = 90, height = 60, units = "mm", device = cairo_pdf)

ggsave(
  "output/p_cum_roh.pdf", 
  plots$cum_roh,
  width = 90,
  height = 80,
  units = "mm"
)

# Statistics
tbs$roh_summary %>%
  filter(
    species == "fuscata"
  ) %>%
  t.test(
    froh_total ~ subspecies,
    data = .
  )

tbs$roh_summary %>%
  filter(
    species == "fuscata"
  ) %>%
  t.test(
    froh_long ~ subspecies,
    data = .
  )

tbs$pixy_king_pruned <- tbs$roh_summary %>%
  select(IID) %>%
  left_join(
    tbs$pixy,
    by = c("IID" = "seq_id")
  )

tbs$pixy_king_pruned %>%
  filter(
    subspecies == "yakui"
  ) %>%
  mutate(
    status = ifelse(str_starts(IID, "Mfy"), "captive", "wild")
  ) %>%
  t.test(
    genome_avg_pi ~ status,
    data = .
  )

tbs$pixy_king_pruned %>%
  filter(
    species == "fuscata"
  ) %>%
  t.test(
    genome_avg_pi ~ subspecies,
    data = .
  )

tbs$roh_summary %>%
  filter(
    subspecies == "yakui"
  ) %>%
  mutate(
    status = ifelse(str_starts(IID, "Mfy"), "captive", "wild")
  ) %>%
  t.test(
    froh_total ~ status,
    data = .
  )