library(tidyverse)

# Input paths
input_path <- list(
  xchrom_female = "../mask/detect_nonpar/sexchrom/xchrom_female.bed.gz",
  xchrom_male = "../mask/detect_nonpar/sexchrom/xchrom_male.bed.gz",
  ychrom_female = "../mask/detect_nonpar/sexchrom/ychrom_female.bed.gz",
  ychrom_male = "../mask/detect_nonpar/sexchrom/ychrom_male.bed.gz"
)

# Load data
depth_data <- list(
  xchrom_female = read_delim(gzfile(input_path$xchrom_female), 
                              col_names=c("chromosome", "start", "end", "depth"), delim=" "),
  xchrom_male = read_delim(gzfile(input_path$xchrom_male), 
                            col_names=c("chromosome", "start", "end", "depth"), delim=" "),
  
  ychrom_female = read_delim(gzfile(input_path$ychrom_female), 
                              col_names=c("chromosome", "start", "end", "depth"), delim=" "),
  ychrom_male = read_delim(gzfile(input_path$ychrom_male), 
                            col_names=c("chromosome", "start", "end", "depth"), delim=" ")
)

depth_data$xchrom_stats <- depth_data$xchrom_female %>%
  left_join(
    depth_data$xchrom_male,
    by=c("chromosome", "start", "end"),
    suffix=c("_female", "_male")
  ) %>%
  mutate(
    depth_sub = depth_female - depth_male,
    depth_ratio = depth_female / depth_male
  ) 

plots <- list(
  depth_sub = depth_data$xchrom_stats %>% 
    ggplot(
      aes(x=start, y = depth_sub)
    ) +
    geom_point() +
    theme_minimal()
)


ggsave("output/depth_sub.pdf",
       plot = plots$depth_sub,
       width = 200, 
       height = 200, 
       units = "mm", device = cairo_pdf)

depth_data$xchrom_stats %>% 
  write_csv("output/xchrom_stats.csv")