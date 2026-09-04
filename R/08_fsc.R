library(tidyverse)

# Input paths
input_path <- list(
  model_comp = "../stats/fsc/model_comp",
  best_var = "../stats/fsc/best",
  bs = "../stats/fsc/bs"
)

# Define variables
config <- list(
  generation = 10.4
)

# Load data
# Model comparisons
tbs <- list(
  model_comp = list.files(
    input_path$model_comp, 
    pattern = "\\.bestlhoods$", 
    full.names = TRUE,
    recursive = TRUE
  ) %>%
    map_df(~{
      data <- read_table(
        .x, 
        show_col_types = FALSE
      )
      
      metadata <- str_match(.x, "model(\\d+)/run(\\d+)/")
      
      data %>%
        mutate(
          model = as.integer(metadata[, 2]),
          run = as.integer(metadata[, 3])
        )
    }) 
)

tbs$model_comp_summary <- tbs$model_comp %>%
  mutate(
    k = case_when(
      model == 0 ~ 3,
      model == 1 ~ 5,
      model == 2 ~ 5,
      model == 3 ~ 6,
      model == 4 ~ 7,
      TRUE ~ NA_integer_
    ),
    aic = 2 * k - 2 * MaxEstLhood/log10(exp(1)), # convert log10 likelihood to natural log likelihood for AIC calculation
    delta_lhood = MaxObsLhood - MaxEstLhood
  ) %>%
  group_by(model) %>%
  slice_max(order_by = MaxEstLhood, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    delta_aic = aic - min(aic)
  ) %>%
  mutate(
    across(starts_with("N"), ~ . / 2), # convert to diploid individuals
    across(starts_with("T"), ~ . * config$generation) # convert to years
  )

# Likelihood distributions
tbs$best_var <- list.files(
  input_path$best_var,
  pattern = "\\_maxL.lhoods$",
  full.names = TRUE,
  recursive = TRUE
) %>%
  map_df(~{
    data <- read_table(
      .x,
      show_col_types = FALSE
    )
    
    metadata <- str_match(.x, "model(\\d+)/run(\\d+)/")
    
    data %>%
      mutate(
        model = as.integer(metadata[, 2]),
        run = as.integer(metadata[, 3])
      )
  }) 

plots <- list(
  lhood =  tbs$best_var %>%
    ggplot(
      aes(
        x = factor(model), y = DAFLHood_1
      )
    ) +
    geom_boxplot() +
    labs(
      x = "Model",
      y = "log Likelihood"
    ) +
    theme_minimal()
)

# Bootstrap
tbs$bs <- list.files(
  input_path$bs, 
  pattern = "\\.bestlhoods$", 
  full.names = TRUE,
  recursive = TRUE
) %>%
  map_df(~{
    data <- read_table(
      .x, 
      show_col_types = FALSE
    )
    
    metadata <- str_match(.x, "bs_(\\d+)/run(\\d+)/")
    
    data %>%
      mutate(
        bs = as.integer(metadata[, 2]),
        run = as.integer(metadata[, 3])
      )
  }) 

# Bootstrap summary for the best model (model 2)
tbs$bs_summary <- tbs$bs %>%
  group_by(bs) %>%
  slice_max(order_by = MaxEstLhood, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    across(starts_with("N"), ~ . / 2), # convert to diploid individuals
    across(starts_with("T"), ~ . * config$generation), # convert to years
    delta_lhood = MaxObsLhood - MaxEstLhood
  ) %>%
  pivot_longer(
    cols = c(starts_with("N"), starts_with("T"), starts_with("M"), starts_with("d")),
    names_to = "parameter",
    values_to = "estimate"
  ) %>%
  group_by(
    parameter
  ) %>%
  summarise(
    mean = mean(estimate),
    sd = sd(estimate),
    ci_lower = quantile(estimate, 0.025),
    ci_upper = quantile(estimate, 0.975),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  left_join(
    tbs$model_comp_summary %>%
      filter(
        model == 2
      ) %>%
      select(
        NCUR, TANC, NANC, MaxEstLhood, TENL, NENL, delta_lhood
      ) %>%
      pivot_longer(
        cols = c(NCUR, TANC, NANC, MaxEstLhood, MaxEstLhood, TENL, NENL, delta_lhood),
        names_to = "parameter",
        values_to = "estimate"
      ),
    by = "parameter"
  )

write_csv(
  tbs$bs_summary,
  "output/bs_summary.csv"
)

write_csv(
  tbs$model_comp_summary,
  "output/model_comp_summary.csv"
)

ggsave(
  "output/lhood.pdf", 
  plots$lhood,
  width = 160,
  height = 80,
  units = "mm",
  device = cairo_pdf
)