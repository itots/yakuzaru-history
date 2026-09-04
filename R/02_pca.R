library(tidyverse)
library(forcats)
library(paletteer)
library(AssocTests)

# Input paths
input_path <- list(
  info = "../list/samples_genome.csv",
  pop = "../list/population.csv",
  eigenval = "../stats/data/plink/out/pca.eigenval",
  eigenvec = "../stats/data/plink/out/pca.eigenvec"
)

# Color palettes
color <- list(
  three = c("#1784bf", "#d1990d", "#07360c")
)

# Load data
tbs <- list(
  info = read_csv(input_path$info),
  pop = read_csv(input_path$pop),
  eigenval = read_tsv(input_path$eigenval, col_names = "eigenvalue")
)
tbs$eigenvec <- read_delim(input_path$eigenvec, 
                          delim = " ",
                          col_names = c("no", "seq_id", paste0("PC", 1:10))
                          ) %>%
  left_join(
    tbs$info %>%
      mutate(
        population = factor(population, 
                            levels = tbs$pop$population),
        cluster = factor(cluster, 
                         levels = c("East", "West", "Yakushima"))
      ),
    by = "seq_id"
  )

config <- list(
  pc_contribution = tbs$eigenval/sum(tbs$eigenval)
)

# Statistics
results <- 
  list(
    tw_005 = tw(
      eigenvalues = tbs$eigenval$eigenvalue,
      eigenL = length(tbs$eigenval$eigenvalue),
      criticalpoint = 0.9793
    ),
    tw_001 = tw(
      eigenvalues = tbs$eigenval$eigenvalue,
      eigenL = length(tbs$eigenval$eigenvalue),
      criticalpoint = 2.0234
    )
  )

# a numeric value corresponding to the significance level. 
# If the significance level is 0.05, 0.01, 0.005, or 0.001, 
# the criticalpoint should be set to be 0.9793, 2.0234, 2.4224, or 3.2724, accordingly. The default is 2.0234.


# PCA plot
plots<- list(
  pca = tbs$eigenvec %>%
    ggplot(
      aes(
        x = PC1, 
        y = PC2,
        colour = cluster,
        shape = population
      )
    ) + 
    geom_point(size = 3, stroke = 0.5) +
    scale_shape_manual(
      values = 0:11
    ) +
    scale_color_manual(
      values = color$three
    ) +
    labs(
      x = paste0("PC1 (", round(config$pc_contribution[1,1]*100, 1), "%)"),
      y = paste0("PC2 (", round(config$pc_contribution[2,1]*100, 1), "%)")
    ) +
    theme_minimal() +
    theme(
      legend.position = "none"
    )
)

ggsave("output/pca.pdf",
       plot = plots$pca,
       width = 110, 
       height = 80, 
       units = "mm", device = cairo_pdf)

writeLines(capture.output(results$tw_001), "output/tw_001.txt")
writeLines(capture.output(results$tw_005), "output/tw_005.txt")

