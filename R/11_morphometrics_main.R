library(geomorph)
library(Morpho)
library(tidyverse)
library(LOST)
library(purrr)
library(ggbeeswarm)
library(brms)
library(vcvComp)
library(car)
library(cowplot)
library(evolqg)
library(VCA)
library(furrr)
library(ggrepel)

source("11_morphometrics_functions.R")

# Input path
input_path <- list(
  coords = "../stats/data/coords",
  original = "../list/samples_morph.csv",
  pop = "../list/population_morph.csv",
  fst = "../stats/fst/out/fst_wc.fst.summary"
)

# Color palettes
color <- list(three = c("#1784bf", "#d1990d", "#07360c"))

# Analysis parameters (Magic numbers)
params <- list(
  beta_pc_n = list(female = 22, male = 23),
  pst_c_h2_list = seq(0, 2, by = 0.1),
  pst_n_boot = 1000,
  t = 22437,
  ne = 7893
)

# Prepare cache directory
if (!dir.exists("output/cache")) dir.create("output/cache", recursive = TRUE)

# Import associated information
tbs <- list(
  original = read_csv(input_path$original),
  pop = read_csv(input_path$pop)
)

tbs$original_pop <- tbs$original %>% 
  left_join(
    tbs$pop %>%
      select(-region),
    by = "population"
  )

# Read landmarks of 4 replicates
coords <- list(raw = imap(set_names(1:4, paste0("rep_", 1:4)), ~ {
  functions$read_morphoj(paste0(input_path$coords, "/land_", .x, ".txt"))
}))

# Estimate missing landmarks for each replicate
coords$fix <- map(coords$raw,
                  ~ MissingGeoMorph(.x, method = "BPCA", original.scale = TRUE))

# Split coordinates by sex
coords$sex <- map(set_names(c("female", "male")), ~ {
  morph_ids <- tbs$original %>% filter(sex == .x) %>% .$morph_id
  map(set_names(1:4, paste0("rep_", 1:4)), ~ coords$fix[[paste0("rep_", .x)]][, , morph_ids])
})

# Evaluate measurement error and bilateral symmetry
coords$merged <- map(set_names(c("female", "male")), ~ {
  bindArr(
    coords$sex[[.x]]$rep_1,
    coords$sex[[.x]]$rep_2,
    coords$sex[[.x]]$rep_3,
    coords$sex[[.x]]$rep_4,
    along = 3
  )
})

gdfs <- list(merged = map(set_names(c("female", "male")), ~ {
  coords_sub <- coords$merged[[.x]]
  morph_ids <- dimnames(coords_sub)[[3]]
  n_inds <- length(morph_ids)
  n_replicates <- dim(coords_sub)[3] / n_inds
  
  geomorph.data.frame(
    shape = coords_sub,
    ind = morph_ids,
    replicate = rep(1:n_replicates, each = n_inds)
  )
}))

results <- list(bilat_sym = map(set_names(c("female", "male")), ~ {
  bilat.symmetry(
    A = shape,
    ind = ind,
    replicate = replicate,
    data = gdfs$merged[[.x]],
    land.pairs = matrix(12:45, ncol = 2),
    object.sym = TRUE,
    Parallel = TRUE
  )
}))

# Centroid size
results$csize <- map(set_names(c("female", "male")), ~ {
  functions$calculate_mean_cs(coords$merged[[.x]])
})

# Principal component analysis
results$pca$sym <- map(set_names(c("female", "male")), ~ {
  gm.prcomp(A = results$bilat_sym[[.x]]$symm.shape)
})
results$pca$sym$female$sdev ^2 / sum(results$pca$sym$female$sdev ^2)
results$pca$sym$female$sdev ^2 %>% cumsum / sum(results$pca$sym$female$sdev ^2)

results$pca$asym <- map(set_names(c("female", "male")), ~ {
  gm.prcomp(A = results$bilat_sym[[.x]]$asymm.shape)
})
results$pca$sym$male$sdev ^2 / sum(results$pca$sym$male$sdev ^2)
results$pca$sym$male$sdev ^2 %>% cumsum / sum(results$pca$sym$male$sdev ^2)

# Allometry
gdfs$sym <- map(set_names(c("female", "male")), ~ {
  shape <- results$bilat_sym[[.x]]$symm.shape
  morph_ids <- dimnames(shape)[[3]]
  csize <- results$csize[[.x]][match(morph_ids, names(results$csize[[.x]]))]
  subspecies <- tbs$original_pop$subspecies[match(morph_ids, tbs$original$morph_id)]
  cluster <- tbs$original_pop$cluster[match(morph_ids, tbs$original$morph_id)]
  
  geomorph.data.frame(
    shape = shape,
    centroid_size = as.vector(csize),
    subspecies = as.factor(subspecies),
    cluster = as.factor(cluster)
  )
})

results$allometry$stats <- map(set_names(c("female", "male")), ~ {
  procD.lm(
    f1 = shape ~ centroid_size + subspecies,
    Parallel = TRUE,
    data = gdfs$sym[[.x]]
  )
})

results$allometry$plots <- map(set_names(c("female", "male")), ~ {
  plotAllometry(
    results$allometry$stats[[.x]],
    size = gdfs$sym[[.x]]$centroid_size,
    logsz = TRUE,
    method = "RegScore"
  )
})

# Count NA and anomaly
tbs$na_anomaly <- tbs$original_pop %>%
  # count NA
  left_join(
    coords$raw$rep_1[, 1, ] %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var = "morph_id") %>%
      as_tibble() %>%
      mutate(
        na_count = rowSums(
          is.na(select(., -morph_id)) 
        )
      ) %>%
      select(morph_id, na_count),
    by = "morph_id"
  ) %>%
  # count anomaly
  mutate(anomaly_count = rowSums(
    is.na(select(
      ., 
      all_of(c(
        "supernumerary_teeth",
        "missing_teeth",
        "rotated_or_abnormaly_oriented_teeth",
        "bone_swelling",
        "bone_depression",
        "bone_deformity_likely_due_to_trauma",
        "broken"
      ))
    )) == FALSE 
  )) %>%
  ungroup()


# Fluctuating asymmetry
results$unsigned_ai <- map(set_names(c("all", "sub")), function(condition) {
  
  # Each sex-specific model
  res_list <- map(set_names(c("female", "male")), function(y) {
    data <- results$bilat_sym[[y]]$unsigned.AI %>%
      enframe(name = "morph_id", value = "unsigned_ai") %>%
      left_join(tbs$na_anomaly, by = "morph_id") %>%
      mutate(unsigned_ai = (unsigned_ai - 1) * 1000) %>%
      { if (condition == "sub") filter(., na_count == 0, anomaly_count == 0) else . } %>%
      group_by(population) %>%
      filter(n() > 3) %>%
      ungroup()
    
    model <- brm(
      formula = unsigned_ai ~ subspecies + (1 | population),
      data = data,
      family = lognormal(),
      prior = c(
        set_prior("student_t(3, 0, 5)", class = "b"),
        set_prior("student_t(3, 0, 5)", class = "sd", lb = 0),
        set_prior("student_t(3, 0, 5)", class = "sigma", lb = 0)
      ),
      cores = 10,
      iter = 4000,
      warmup = 2000,
      control = list(adapt_delta = 0.9999, max_treedepth = 30),
      file = paste0("output/cache/brm_ai_", condition, "_", y) # brms built-in cache
    )
    
    list(model = model, summary = summary(model))
  })
  
  # Sex-combined model
  res_list$combined <- map_df(set_names(c("female", "male")), function(y) {
    results$bilat_sym[[y]]$unsigned.AI %>%
      enframe(name = "morph_id", value = "unsigned_ai") %>%
      left_join(tbs$na_anomaly, by = "morph_id") %>%
      mutate(unsigned_ai = (unsigned_ai - 1) * 1000, sex = y)
  }) %>%
    { if (condition == "sub") filter(., na_count == 0, anomaly_count == 0) else . } %>%
    group_by(sex, population) %>% 
    filter(n() > 3) %>% 
    ungroup() %>%
    {
      model <- brm(
        formula = unsigned_ai ~ sex * subspecies + (1 | population),
        data = .,
        family = lognormal(),
        prior = c(
          set_prior("student_t(3, 0, 5)", class = "b"),
          set_prior("student_t(3, 0, 5)", class = "sd", lb = 0),
          set_prior("student_t(3, 0, 5)", class = "sigma", lb = 0)
        ),
        cores = 10,
        iter = 4000,
        warmup = 2000,
        control = list(adapt_delta = 0.99995, max_treedepth = 30),
        file = paste0("output/cache/brm_ai_", condition, "_combined") # brms built-in cache
      )

      list(model = model, summary = summary(model))
    }
  
  return(res_list)
})


# Within-group PCA
results$wpca <- map(set_names(c("female", "male")), ~ {
  tb <- gdfs$sym[[.x]]$shape %>%
    two.d.array() %>%
    as.data.frame() %>%
    rownames_to_column(var = "morph_id") %>%
    left_join(tbs$original, by = "morph_id")
  
  tb_w <- tb %>%
    group_by(population) %>%
    filter(n() > 20) %>% # Shimokita, Tone, and Yakushima pooled
    ungroup()
  
  tb_b <- tb %>%
    group_by(population) %>%
    filter(n() > 3) %>% # only the population with N > 3 used
    summarize_at(paste0("V", 1:135), mean) %>%
    ungroup()
  
  data_b <- tb_b %>%
    select(paste0("V", 1:135)) %>%
    as.matrix()
  
  data_w <- tb_w %>%
    select(paste0("V", 1:135)) %>%
    as.matrix()
  
  cov_w <- cov.W(X = data_w, groups = tb_w$population)
  
  eigen_w <- eigen(cov_w)
  
  pc_w <- data_w %*% eigen_w$vectors
  pc_b <- data_b %*% eigen_w$vectors
  
  list(
    values_w = eigen_w$values,
    vectors_w = eigen_w$vectors,
    cov_w = cov_w,
    tb_w = tb_w,
    tb_b = tb_b,
    pc_w = pc_w,
    pc_b = pc_b,
    data_w = data_w,
    data_b = data_b
  )
})

# Beta test
round(results$wpca$female$values_w %>% cumsum() * 100 / sum(results$wpca$female$values_w), 1)
round(results$wpca$female$values_w * 100 / sum(results$wpca$female$values_w), 1)
round(results$wpca$male$values_w %>% cumsum() * 100 / sum(results$wpca$male$values_w), 1)
round(results$wpca$male$values_w * 100 / sum(results$wpca$male$values_w), 1)

# The first 22 (female) or 23 (male) PCs, which explain more than 90% of total variance, used
results$beta_test <- map(set_names(c("female", "male")), ~ {
  pc_n <- params$beta_pc_n[[.x]]
  
  var_b <- results$wpca[[.x]]$pc_b %>% apply(2, var) %>% .[1:pc_n] %>% log()
  var_w <- results$wpca[[.x]]$values_w %>% .[1:pc_n] %>% log()
  
  model_fit <- lm(var_b ~ var_w)
  summary <- summary(model_fit)
  slope_test <- linearHypothesis(model_fit, "var_w = 1")
  
  list(
    var_w = var_w,
    var_b = var_b,
    model_fit = model_fit,
    summary = summary,
    slope_test = slope_test
  )
})


# MultivDriftTest
round(results$pca$sym$female$sdev^2 * 100 / sum(results$pca$sym$female$sdev^2), 1)
round(results$pca$sym$female$sdev^2 %>% cumsum() * 100  / sum(results$pca$sym$female$sdev^2), 1)
round(results$pca$sym$male$sdev^2 * 100 / sum(results$pca$sym$male$sdev^2), 1)
round(results$pca$sym$male$sdev^2 %>% cumsum() * 100  / sum(results$pca$sym$male$sdev^2), 1)

results$multivDriftTest <- map(set_names(c("female", "male")), ~ {
  
  cache_file <- paste0("output/cache/multivdrift_", .x, ".rds")
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  
  pc_score <- results$pca$sym[[.x]]$x %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "morph_id") %>% 
    left_join(tbs$original_pop, by = "morph_id")
  
  fuscata <- pc_score %>% 
    filter(subspecies == "fuscata") %>%
    select(paste0("Comp",1:25))
  
  yakui <- pc_score %>% 
    filter(subspecies == "yakui") %>%
    select(paste0("Comp",1:25))
  
  W_data <- pc_score %>%
    group_by(population) %>%
    filter(n() > 20) %>%
    ungroup()
  
  W <- cov.W(X = W_data %>% select(paste0("Comp", 1:25)), 
             groups = W_data$population) 

  G <- W
  
  res <- MultivDriftTest(
    population1 = fuscata,
    population2 = yakui,
    G = G,
    Ne = params$ne,
    generations = params$t,
    iterations = 10000
  )
  
  saveRDS(res, cache_file)
  return(res)
})

results$multivDriftTest$female
results$multivDriftTest$male

# PST-FST comparison (with relative path & explicit cache)
tbs$fst <- read_tsv(input_path$fst,
         skip = 1,
         col_names = c("pop1", "pop2", "fst"))

plan(multisession, workers = parallel::detectCores() - 1)

results$pst <- map(set_names(c("female", "male")), ~ {
  sex_id <- .x
  c_h2_list <- params$pst_c_h2_list
  n_boot <- params$pst_n_boot
  
  data <- results$pca$sym[[sex_id]]$x %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "morph_id") %>% 
    left_join(
      results$csize[[sex_id]] %>%
        enframe(name = "morph_id", value = "centroid_size"),
      by = "morph_id"
    ) %>%
    left_join(tbs$original_pop, by = "morph_id") %>%
    filter(cluster == "West" | cluster == "Yakushima") %>% rename_with(~ str_replace(., "Comp", "PC"))

  original_vc <- functions$get_vc_core(
    data,
    target_vars = c("centroid_size", paste0("PC", 1:4))
    )
  
  # Bootstrap caching
  cache_boot_file <- paste0("output/cache/boot_vc_list_", sex_id, ".rds")
  if (file.exists(cache_boot_file)) {
    boot_vc_list <- readRDS(cache_boot_file)
  } else {
    boot_vc_list <- future_map(1:n_boot, function(i) {
      boot_data <- data %>% 
        dplyr::group_by(cluster) %>% 
        dplyr::sample_frac(replace = TRUE) %>% 
        dplyr::ungroup()
      functions$get_vc_core(
        boot_data,
        target_vars = c("centroid_size", paste0("PC", 1:4))
      )
    }, .options = furrr_options(
      seed = TRUE, 
      packages = c("tidyverse"), 
      globals = TRUE)) %>% 
      dplyr::bind_rows()
      
    saveRDS(boot_vc_list, cache_boot_file)
  }
  
  results_pst <- map_dfr(c_h2_list, function(c_h2) {
    
    boot_pst_summary <- boot_vc_list %>%
      dplyr::mutate(pst = vc_b * c_h2 / (vc_b * c_h2 + 2 * vc_w)) %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(
        pst_lower = quantile(pst, 0.025, na.rm = TRUE),
        pst_upper = quantile(pst, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
    
    original_vc %>%
      dplyr::mutate(
        c_h2 = c_h2,
        pst_orig = vc_b * c_h2 / (vc_b * c_h2 + 2 * vc_w)
      ) %>%
      dplyr::left_join(boot_pst_summary, by = "variable")
  })
  
  return(
    results_pst
    )
})


# bgPCA — between subspecies
results$bg_pca <- map(set_names(c("female", "male")), ~ {
  groupPCA(
    dataarray = gdfs$sym[[.x]]$shape,
    groups = gdfs$sym[[.x]]$subspecies,
    cv = FALSE
  )
})

# bgPCA — among populations
results$bg_pca_pop <- map(set_names(c("female", "male")), ~ {
  tb <- gdfs$sym[[.x]]$shape %>%
    two.d.array() %>%
    as.data.frame() %>%
    rownames_to_column(var = "morph_id") %>%
    left_join(tbs$original, by = "morph_id") %>%
    group_by(population) %>%
    filter(n() > 2) %>%
    ungroup()
  
  groupPCA(
    dataarray = gdfs$sym[[.x]]$shape[,,tb$morph_id],
    groups = as.factor(tb$population),
    cv = FALSE
  )
})


# Save results
results$bilat_sym$female$shape.anova %>% 
  write_csv("output/female_bilat_sym_shape_anova.csv")

results$bilat_sym$male$shape.anova %>% 
  write_csv("output/male_bilat_sym_shape_anova.csv")

results$allometry$stats$female$aov.table %>%
  write_csv("output/female_allometry_aov_table.csv")

results$allometry$stats$male$aov.table %>%
  write_csv("output/male_allometry_aov_table.csv")

bind_rows(
  results$unsigned_ai$all$female$summary$random$population,
  results$unsigned_ai$all$female$summary$fixed,
  results$unsigned_ai$all$female$summary$spec_pars
) %>% 
  write_csv("output/female_unsigned_ai_all_summary.csv")

bind_rows(
  results$unsigned_ai$all$combined$summary$random$population,
  results$unsigned_ai$all$combined$summary$fixed,
  results$unsigned_ai$all$combined$summary$spec_pars
) %>% 
  write_csv("output/combined_unsigned_ai_all_summary.csv")

bind_rows(
  results$unsigned_ai$sub$female$summary$random$population,
  results$unsigned_ai$sub$female$summary$fixed,
  results$unsigned_ai$sub$female$summary$spec_pars
) %>% 
  write_csv("output/female_unsigned_ai_sub_summary.csv")

bind_rows(
  results$unsigned_ai$all$male$summary$random$population,
  results$unsigned_ai$all$male$summary$fixed,
  results$unsigned_ai$all$male$summary$spec_pars
) %>% 
  write_csv("output/male_unsigned_ai_all_summary.csv")

bind_rows(
  results$unsigned_ai$sub$male$summary$random$population,
  results$unsigned_ai$sub$male$summary$fixed,
  results$unsigned_ai$sub$male$summary$spec_pars
) %>% 
  write_csv("output/male_unsigned_ai_sub_summary.csv")

bind_rows(
  results$unsigned_ai$sub$combined$summary$random$population,
  results$unsigned_ai$sub$combined$summary$fixed,
  results$unsigned_ai$sub$combined$summary$spec_pars
) %>% 
  write_csv("output/combined_unsigned_ai_sub_summary.csv")


