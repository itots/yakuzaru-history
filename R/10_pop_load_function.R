library(tidyverse)
library(purrr)
library(readr)

fnc <- list()

# Helper function to read genome sizes and calculate cumulative offsets
fnc$read_genome_sizes <- function(file_path) {
  read_tsv(file_path, show_col_types = FALSE) %>%
    # Autosomes only
    filter(`Sequence-Name` %in% paste0("chr", 1:20)) %>%
    select(
      chrom_no = `RefSeq-Accn`,
      length = `Sequence-Length`,
      seq_name = `Sequence-Name`
    ) %>%
    mutate(seq_name = fct_relevel(seq_name, paste0("chr", 1:20))) %>%
    arrange(seq_name) %>%
    mutate(offset = lag(cumsum(length), default = 0)) %>%
    mutate(chrom_no = as.character(chrom_no)) %>%
    select(-seq_name)
}


# 1. Calculate f values
fnc$make_f <- function(tb, eff) {
  wide <- tb %>%
    filter(type == "SNP", effect == eff) %>%
    select(chrom_no, position, af, population) %>%
    pivot_wider(values_from = af, names_from = "population", values_fill = list(af = 0))

  wide %>%
    mutate(
      f_ew = East * (1 - West),
      f_ey = East * (1 - Yakushima),
      f_we = West * (1 - East),
      f_wy = West * (1 - Yakushima),
      f_ye = Yakushima  * (1 - East),
      f_yw = Yakushima  * (1 - West)
    )
}


# 2. Leave-One-Block-Out Rxy estimation
fnc$lobo_rxy <- function(eff, int, B = 100, chr_lens = NULL) {
  sumcols <- c("f_ew","f_we","f_ey","f_ye","f_wy","f_yw")

  eff <- eff %>% mutate(chrom_no = as.character(chrom_no))
  int <- int %>% mutate(chrom_no = as.character(chrom_no))

  len_df <- chr_lens %>% mutate(chrom_no = as.character(chrom_no))

  eff <- eff %>% 
    left_join(len_df %>% select(chrom_no, offset), by = "chrom_no") %>%
    mutate(cum_pos = position + offset)
    
  int <- int %>% 
    left_join(len_df %>% select(chrom_no, offset), by = "chrom_no") %>%
    mutate(cum_pos = position + offset)

  total_genome_size <- sum(len_df$length)
  block_size <- total_genome_size / B
  
  eff$block_raw <- as.integer(ceiling(eff$cum_pos / block_size))
  int$block_raw <- as.integer(ceiling(int$cum_pos / block_size))
  
  valid_blocks <- sort(unique(c(eff$block_raw, int$block_raw)))
  B0 <- length(valid_blocks)
  
  eff$block <- match(eff$block_raw, valid_blocks)
  int$block <- match(int$block_raw, valid_blocks)
  
  as_numvec <- function(one_row_tbl) {
    v <- as.numeric(one_row_tbl[1, ])
    names(v) <- names(one_row_tbl)
    v
  }
  
  seff_full <- eff %>% summarise(across(all_of(sumcols), ~ sum(.x, na.rm = TRUE))) %>% select(all_of(sumcols)) %>% as_numvec()
  sint_full <- int %>% summarise(across(all_of(sumcols), ~ sum(.x, na.rm = TRUE))) %>% select(all_of(sumcols)) %>% as_numvec()

  sint_full[sint_full <= 0 | is.na(sint_full)] <- NA_real_
  L_full <- seff_full / sint_full
  rxy_full <- tibble(
    pair = c("ew","ey","wy"),
    Rxy  = c(L_full["f_ew"]/L_full["f_we"],
             L_full["f_ey"]/L_full["f_ye"],
             L_full["f_wy"]/L_full["f_yw"])
  )
  
  rxy_lobo <- map_dfr(seq_len(B0), function(b) {
    seff <- eff %>% filter(block != b) %>% summarise(across(all_of(sumcols), ~ sum(.x, na.rm = TRUE))) %>% select(all_of(sumcols)) %>% as_numvec()
    sint <- int %>% filter(block != b) %>% summarise(across(all_of(sumcols), ~ sum(.x, na.rm = TRUE))) %>% select(all_of(sumcols)) %>% as_numvec()
    sint[sint <= 0 | is.na(sint)] <- NA_real_
    Lb <- seff / sint
    tibble(
      rep = b, pair = c("ew","ey","wy"),
      Rxy = c(Lb["f_ew"]/Lb["f_we"], Lb["f_ey"]/Lb["f_ye"], Lb["f_wy"]/Lb["f_yw"]) 
    )
  })
  
  summary <- rxy_lobo %>%
    group_by(pair) %>%
    summarise(
      jk_mean = mean(Rxy, na.rm = TRUE),
      jk_var  = ((B0 - 1)/B0) * sum((Rxy - jk_mean)^2, na.rm = TRUE),
      se      = sqrt(jk_var),
      .groups = "drop"
    ) %>%
    left_join(rxy_full, by = "pair") %>%
    transmute(pair, Rxy_full = Rxy, se,
              ci_lower = Rxy_full - 1.96*se,
              ci_upper = Rxy_full + 1.96*se,
              B = B0)
              
  list(summary = summary, replicates = rxy_lobo, full = rxy_full)
}


# 3. Summarise SFS
fnc$project_and_summarise_sfs <- function(tb, k = 14, reps = 100, seed = NULL) {
  stopifnot(k > 0, reps >= 1)

  df <- tb %>%
    filter(!is.na(ac), !is.na(an), an >= k) %>%
    mutate(
      m_white = pmax(as.integer(ac), 0),
      n_black = pmax(as.integer(an - ac), 0)
    )

  if (!is.null(seed)) set.seed(seed)

  out_counts <- map_dfr(seq_len(reps), function(r) {
    sampled_ac <- rhyper(nn = nrow(df), m = df$m_white, n = df$n_black, k = k)
    
    tibble(
      population = df$population,
      effect     = df$effect,
      dac        = as.integer(sampled_ac) # downsampled ac
    ) %>%
    count(population, effect, dac, name = "count") %>%
    mutate(rep = as.integer(r))
  })

  summary_stats <- out_counts %>%
    filter(dac > 0, dac < k) %>%
    group_by(rep, population, effect) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup() %>%
    group_by(population, effect, dac) %>%
    summarise(
      mean = mean(prop, na.rm = TRUE),
      lo   = quantile(prop, 0.025, na.rm = TRUE),
      hi   = quantile(prop, 0.975, na.rm = TRUE),
      nrep = dplyr::n_distinct(rep),
      .groups = "drop"
    )
    
  return(summary_stats)
}
