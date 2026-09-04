functions <- list()

# Read MorphoJ format landmark data
functions$read_morphoj <- function(file_path){
  tb_data <- read_tsv(file_path)
  id_list <- tb_data$id.string
  array_landmarks <- arrayspecs(
    as.matrix(
      tb_data[,-1]
    ),
    p = 45,
    k = 3
  )
  dimnames(array_landmarks)[[3]] <- id_list
  return(array_landmarks)
}

# Calculate specimen mean centroid size
functions$calculate_mean_cs <- function(landmark_array){
  gpa <- gpagen(landmark_array)
  centroid_size <- gpa$Csize
  mean_centroid_size <- tapply(centroid_size, names(centroid_size), mean)
  return(mean_centroid_size)
}

# Function for calculating VC 
functions$get_vc_core <- function(df, target_vars) {
  map_dfr(set_names(target_vars), function(var_name) {
    form <- as.formula(paste(var_name, "~ cluster"))
    vca_out <- suppressMessages(VCA::anovaVCA(form, Data = as.data.frame(df)))
    tibble(vc_b = vca_out$VCoriginal[1], vc_w = vca_out$VCoriginal[2])
  }, .id = "variable")
}