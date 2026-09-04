library(tidyverse)
library(cowplot)

# Input path
input_path <- list(
  tree_a = "../stats/dsuite/a/DTparallel_sample_dsuite_a_combined_tree.txt",
  tree_x = "../stats/dsuite/x/sample_dsuite_x_tree.txt",
  fbranch_a = "../stats/dsuite/a/Fbranch.txt",
  fbranch_x = "../stats/dsuite/x/Fbranch.txt"
)

tbs <- list(
  tree_a = read_tsv(input_path$tree_a),
  tree_x = read_tsv(input_path$tree_x),
  fbranch_a = read_tsv(input_path$fbranch_a),
  fbranch_x = read_tsv(input_path$fbranch_x)
)

params <- list(
  fus_pops = c("fus_Yakushima", "fus_Kojima", "fus_Takasakiyama", "fus_Nakatosa", "fus_Okazaki", "fus_Minoh",
              "fus_Tsubaki", "fus_Hagachizaki", "fus_Jigokudani", "fus_Yamagata", "fus_Shimokita")
)

tbs$tree <- 
bind_rows(
  tbs$tree_a %>%
    mutate(
      chromosome = "A"
    ) %>%
    relocate(
      chromosome
    ),
  tbs$tree_x %>%
    mutate(
      chromosome = "X"
    ) %>%
    relocate(
      chromosome
    )
) 


tbs$fbranch <- bind_rows(
  tbs$fbranch_a %>% 
    mutate(
      chromosome = "A"
    ) %>%
    relocate(
      chromosome
    ),
  tbs$fbranch_x %>%
    mutate(
      chromosome = "X"
    ) %>%
    relocate(
      chromosome
    )
) %>%
  pivot_longer(
    cols = -c(1:3),
    names_to = "P3",
    values_to = "fbranch"
  ) %>%
  mutate(
    chromosome = factor(chromosome, levels = c("X", "A"))
  )

plots <- list()
plots$RMC_TM_JM <-
tbs$fbranch %>%
  filter(
    (branch == "b8" | branch == "b5"), # the ancestral lineage of Chinese rhesus
    str_starts(P3, "fus_")
  ) %>%
  mutate(
    branch = if_else(
      branch == "b8",
      branch,
      branch_descendants
    )
  ) %>%
  mutate(
    P3 = factor(
      P3,
      levels = params$fus_pops
    )
  ) %>%
  ggplot(
    aes(
      x = P3,
      y = fbranch,
      fill = chromosome
    )
  ) + geom_bar(
    stat = "identity",
    position = "dodge",
    width = 0.6
  ) +
  scale_fill_manual(
    values = c(
      "A" = "grey70",
      "X" = "grey10"
    ),
    breaks = c("A", "X")
  ) +
  theme_minimal() + 
  coord_flip() +
  facet_wrap(
    ~ branch,
    ncol = 2,
    scales = "fixed"
  )

plots$Yaku <-
  tbs$fbranch %>%
  filter(
    P3 == "fus_Yakushima",
    str_starts(branch_descendants, "fus")
  ) %>%
  group_by(branch, branch_descendants, P3) %>%
  filter(any(fbranch > 0)) %>%
  ungroup() %>%
  ggplot(
    aes(
      x = branch,
      y = fbranch,
      fill = chromosome
    )
  ) + geom_bar(
    stat = "identity",
    position = "dodge",
    width = 0.6
  ) +
  scale_fill_manual(
    values = c(
      "A" = "grey70",
      "X" = "grey10"
    ),
    breaks = c("A", "X")
  ) +
  theme_minimal() + 
  coord_flip()




tbs$fbranch %>%
  filter(
    branch == "b8", # the ancestral lineage of Chinese rhesus
    str_starts(P3, "fus_")
  ) %>%
  group_by(
    chromosome
  ) %>%
  summarize(
    min = min(fbranch),
    max = max(fbranch)
  )

tbs$fbranch %>%
  filter(
    P3 == "fus_Yakushima"
  ) %>%
  arrange(
    -fbranch
  )

ggsave(
  "output/yaku.pdf",
  plots$Yaku + theme(
    legend.position = "none"
  ),
  width = 70,
  height = 45,
  units = "mm"
)
ggsave(
  "output/RMC_TM.pdf",
  plots$RMC_TM_JM + theme(
    legend.position = "none"
  ),
  width = 140,
  height = 100,
  units = "mm"
)
ggsave(
  "output/legend.pdf",
  get_legend(
    plots$RMC_TM_JM + theme(
      legend.position = "right",
      legend.box.margin = margin(0,0,0,0)
    )
  ),
  width = 20,
  height = 45,
  units = "mm"
)

