plots <- list()

# PCA
plots$pca <- map(set_names(c("sym", "asym")), ~ {  
  pca_type <- .x  
  map(set_names(c("female", "male")), ~ {  
    sex <- .x  
    map(list(Comp1_2 = c("Comp1", "Comp2"), Comp2_3 = c("Comp2", "Comp3")), ~ {  
      comp <- .x  
      ggplot(
        as.data.frame(results$pca[[pca_type]][[sex]]$x) %>%
          rownames_to_column(var = "morph_id") %>%
          left_join(tbs$original_pop, by = "morph_id") %>%
          filter(sex == !!sex) %>%
          mutate(
            subspecies = factor(subspecies, levels = c("yakui", "fuscata")),
          ),
        aes(
          x = .data[[comp[1]]],
          y = .data[[comp[2]]]
        )
      ) +
        geom_point(aes(
          colour = .data[["cluster"]],
          shape = .data[["subspecies"]]
        )) +
        stat_ellipse(
          aes(
          group = .data[["subspecies"]], 
          linetype = .data[["subspecies"]]
        ),
        linewidth = 0.3) +
        scale_color_manual(
          values = color$three
        ) +
        theme_minimal()
    })
  })
})


# Allometry
plots$allometry <- map(
  set_names(c("female", "male")), ~ {
    tibble(
      morph_id = rownames(results$allometry$plots[[.x]]$RegScore),
      lnCS = results$allometry$plots[[.x]]$size.var,
      reg_score = results$allometry$plots[[.x]]$RegScore
    ) %>%
      left_join(
        tbs$original_pop,
        by = "morph_id"
      ) %>%
      mutate(
        subspecies = factor(subspecies, levels = c("yakui", "fuscata"))
      ) %>%
      ggplot(
        aes(x = lnCS, y = reg_score, shape = subspecies)
      ) +
      geom_point(aes(colour = cluster)) +
      stat_ellipse(
        aes(
          group = .data[["subspecies"]], 
          linetype = .data[["subspecies"]]
        ),
        linewidth = 0.3) +
      scale_color_manual(
        values = color$three
      ) +
      labs(
        x = "ln Centroid Size",
        y = "Regression score"
      ) +
      theme_minimal()
  }
)

# Unsigned asymmetry index
plots$unsigned_ai <- map(
  set_names(c("all", "sub")), ~ {
    data <- bind_rows(
      results$bilat_sym$female$unsigned.AI %>%
        enframe(name = "morph_id", value = "unsigned_ai"),
      results$bilat_sym$male$unsigned.AI %>%
        enframe(name = "morph_id", value = "unsigned_ai")
    ) %>%
      left_join(tbs$na_anomaly %>%
                  mutate(
                    population = factor(population, levels = na.omit(tbs$pop$population))
                  ), by = "morph_id") %>%
      mutate(
        unsigned_ai = (unsigned_ai - 1) * 1000
      ) 
    
    if (.x == "sub") {
      data <- data %>% filter(na_count == 0, anomaly_count == 0)
    }
    
    fuscata_mean <- data %>%   
      filter(subspecies == "fuscata") %>% 
      group_by(sex) %>%
      dplyr::summarize(mean_value = mean(unsigned_ai)) 
    
    fuscata_groupmean <- data %>%   
      filter(subspecies == "fuscata") %>% 
      group_by(sex, population) %>%
      filter(n() > 3) %>%
      dplyr::summarize(pop_mean = mean(unsigned_ai, na.rm = TRUE), .groups = "drop") %>%
      group_by(sex) %>%
      dplyr::summarize(mean_value = mean(pop_mean, na.rm = TRUE))
    
    yakui_median <- data %>%   
      filter(subspecies == "yakui") %>% 
      group_by(sex) %>%
      dplyr::summarize(median_value = median(unsigned_ai, na.rm = TRUE))
 
    boxplot <- ggplot(data, aes(x = population, y = unsigned_ai, colour = cluster, shape = subspecies)) +
      geom_boxplot(outlier.shape = NA, fill = "white", linewidth = 0.3, fatten = 1) +  
      geom_beeswarm(aes(group = cluster), dodge.width = 0.75, alpha = 0.75) +
      #geom_hline(data = fuscata_mean, aes(yintercept = mean_value), colour = "#F8766D") +
      #geom_hline(data = fuscata_groupmean, aes(yintercept = mean_value), colour = "#F8766D", linetype = "dashed") +
      geom_hline(data = yakui_median, aes(yintercept = median_value), colour = color$three[3], linetype = "dashed", linewidth = 0.4) +
      scale_color_manual(
        values = color$three
      ) +
      facet_wrap(. ~ sex, ncol = 1, scale = "free_y") +
      labs(
        x = "population",
        y = "Scaled unsigned asymmetry index"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right")
    
    histogram <- ggplot(data, aes(x = unsigned_ai, fill = cluster)) +
      geom_histogram() +
      scale_fill_manual(
        values = color$three
      ) +
      facet_wrap(. ~ sex, ncol = 1) 
    
    loghistogram <- ggplot(data, aes(x = log(unsigned_ai), fill = cluster)) +
      geom_histogram() +
      scale_fill_manual(
        values = color$three
      ) +
      facet_wrap(. ~ sex, ncol = 1) 
    
    list(
      boxplot = boxplot,
      histogram = histogram
    )
  }
)

plots$pca_allometry <- plot_grid(
  plot_grid(
    plots$pca$sym$female$Comp1_2 + theme(legend.position = "none") + 
      labs(
        x = "PC1 (20.4%)",
        y = "PC2 (12.9%)"
      ),
    plots$pca$sym$male$Comp1_2 + theme(legend.position = "none") + 
      labs(
        x = "PC1 (17.2%)",
        y = "PC2 (15.9%)"
      ),
    plots$allometry$female + theme(legend.position = "none"),
    plots$allometry$male + theme(legend.position = "none"),
    ncol = 4
  ),
  get_legend(plots$pca$sym$female$Comp1_2 + 
                 theme(legend.position = "bottom",
                       legend.box = "horizontal",
                       legend.key.size = unit(0.6, "cm"),   
                       legend.text = element_text(size = 8),
                       legend.key.spacing.y = unit(0, "pt")
                 )
  ),
  ncol = 1,
  rel_heights = c(6,1)
)

plots$pca_asym <- plot_grid(

    plots$pca$asym$female$Comp1_2 + theme(legend.position = "none") + 
      labs(
        x = "PC1 (14.0%)",
        y = "PC2 (10.8%)"
      ),
    plots$pca$asym$male$Comp1_2 + theme(legend.position = "none") + 
      labs(
        x = "PC1 (15.8%)",
        y = "PC2 (10.3%)"
      ) ,
  get_legend(plots$pca$asym$female$Comp1_2 + 
               theme(legend.position = "right",
                     legend.box = "vertical",
                     legend.key.size = unit(0.6, "cm"),   
                     legend.text = element_text(size = 8),
                     legend.key.spacing.y = unit(0, "pt")
               )
  ),
  ncol = 3,
  rel_widths = c(2,2,1)
)

results$pca$sym$female$sdev^2/sum(results$pca$sym$female$sdev^2)
results$pca$sym$male$sdev^2/sum(results$pca$sym$male$sdev^2)

results$pca$asym$female$sdev^2/sum(results$pca$asym$female$sdev^2)
results$pca$asym$male$sdev^2/sum(results$pca$asym$male$sdev^2)


# bgPCA — subspecies
plots$bgpca_ssp$female <- results$bg_pca$female$Scores %>%
  as.data.frame() %>%
  rownames_to_column(var = "morph_id") %>%
  left_join(tbs$original_pop, by = "morph_id") %>%
  ggplot() +
  geom_histogram(aes(x = V1, fill = subspecies), position = "identity", alpha = 0.7, bins = 30) +
  scale_fill_manual(values = c("lightgray", "black")) +
  theme_minimal() +
  labs(
    x = "bgPC1 (ssp.)",
    y = "Count"
  )

plots$bgpca_ssp$male <- results$bg_pca$male$Scores %>%
  as.data.frame() %>%
  rownames_to_column(var = "morph_id") %>%
  left_join(tbs$original_pop, by = "morph_id") %>%
  ggplot() +
  geom_histogram(aes(x = V1, fill = subspecies), position = "identity", alpha = 0.7, bins = 30) +
  scale_fill_manual(values = c("darkgray", "black")) +
  theme_minimal() +
  labs(
    x = "bgPC1 (ssp.)",
    y = "Count"
  )

# bgPCA — populations
plots$bgpca_pop$female <- results$bg_pca_pop$female$Scores %>%
  as.data.frame() %>%
  rownames_to_column(var = "morph_id") %>%
  left_join(tbs$original_pop, by = "morph_id") %>%
  ggplot(aes(x = V1, y = V2, color = cluster, shape = population)) +
  geom_point() +
  scale_color_manual(values = color$three) +
  scale_shape_manual(values = 1:15) +
  labs(
    x = "bgPC1",
    y = "bgPC2"
  ) +
  theme_minimal()

plots$bgpca_pop$male <- results$bg_pca_pop$male$Scores %>%
  as.data.frame() %>%
  rownames_to_column(var = "morph_id") %>%
  left_join(tbs$original_pop, by = "morph_id") %>%
  ggplot(aes(x = V1, y = V2, color = cluster, shape = population)) +
  geom_point() +
  scale_color_manual(values = color$three) +
  scale_shape_manual(values = 1:15) +
  labs(
    x = "bgPC1",
    y = "bgPC2"
  ) +
  theme_minimal()


# Beta test
plots$beta_test <- 
  bind_rows(
    tibble(
      sex = "female",
      component = 1:22,
      log_w = results$beta_test$female$var_w,
      log_b = results$beta_test$female$var_b
    ),
    tibble(
      sex = "male",
      component = 1:23,
      log_w = results$beta_test$male$var_w,
      log_b = results$beta_test$male$var_b
    )
  ) %>%
  {
    data <- .
    
    abline_data <- data %>%
      group_by(sex) %>%
      dplyr::summarise(
        mean_log_b = mean(log_b),
        mean_log_w = mean(log_w),
        intercept = mean_log_b - mean_log_w
      )
    
    ggplot(
      data = data,
      aes(x = log_w, y = log_b)
    ) +
      geom_smooth(formula = "y ~ x", method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
      geom_point() +
      geom_text_repel(
        aes(label = component),
        size = 3,
        box.padding = 0.1, 
        point.padding = 0.1,
        max.overlaps = Inf 
      ) +
      geom_abline(
        data = abline_data,
        aes(
          slope = 1,
          intercept = intercept
        ),
        linetype = "dashed"
      ) +
      theme_minimal() +
      labs(
        x = "Log within-population variance",
        y = "Log between-population variance"
      ) +
      facet_wrap(~sex, ncol = 1)
  }
  


# PST-FST
plots$pst <- bind_rows(
  results$pst, .id = "sex"
  ) %>%
  ggplot(aes(x = c_h2, y = pst_orig)) +
  geom_hline(yintercept = as.numeric(tbs$fst[3,3]), color = "orange") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = pst_lower, ymax = pst_upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(sex ~ variable) +
  labs(
    x = expression(italic(c)/italic(h)^2),  
    y = "PST or FST"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


  
# Save
ggsave(
  "output/pca_allometry.pdf",
  plots$pca_allometry,
  width = 260,
  height = 70,
  units = "mm"
)

ggsave(
  "output/pca_asym.pdf",
  plots$pca_asym,
  width = 180,
  height = 70,
  units = "mm"
)

ggsave(
  "output/unsigned_ai_all.pdf",
  plots$unsigned_ai$all$boxplot,
  width = 240,
  height = 120,
  units = "mm"
)

ggsave(
  "output/unsigned_ai_sub.pdf",
  plots$unsigned_ai$sub$boxplot,
  width = 240,
  height = 120,
  units = "mm"
)

ggsave(
  "output/beta_test.pdf",
  plots$beta_test,
  width = 80,
  height = 140,
  units = "mm"
)

ggsave(
  "output/pst.pdf",
  plots$pst,
  width = 140,
  height = 140,
  units = "mm"
)