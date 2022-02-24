#' @export
#'
#' @name transfer
#' @title Create, apply and graph the results of a transfer function
#'
#' @param x Name of your own data to apply the transfer function to
#' @param tf The continent for which you want the transfer function data for (na or eu)
#' @param country The country within the continental data you wish to choose. From the eu dataset these are:
#' "denmark"     "england"     "estonia"     "faroeis."
#' "finland"     "germany"     "ireland"     "n.ireland"
#' "netherlands" "poland"      "russia"      "scotland"
#' "sweden"      "switzerland" "turkey"
#' From the na dataset these are :
#' "canada" "usa"
#' @param save Choose to save the output to disk
#' @param age_file The name of the optional age file to graph the results with
#' @param depth_start The first depth you sampled
#' @param depth_int The interval in cm between depths of each sample
#' @param boot_size The numbe of bootstrap cycles to run
#' @return Performance of the model, reconstruction, zscores and error, graphs of WTD and Zscores
#' @examples
#' transfer(testate_Data = "lh1_tests", tf = "eu", country = "england", save = T, age_file = "lh1_ages", depth_start = 1, depth_int = 1, boot_size = 1000)

transfer <- function(testate_data = "lh1_tests",
                     tf = "eu",
                     country = "england",
                     save = T,
                     age_file = "lh1_ages",
                     depth_start = 1,
                     depth_int = 1,
                     boot_size = 1000) {
  require(rioja)
  require(dplyr)
  require(ggpubr)
  require(vegan)
  require(effectsize)
  name <- testate_data
  site_nme <-  gsub("\\_.*", "", name)
  mydir <- paste0("tf_runs/")
  if (!dir.exists(mydir)){
    dir.create(mydir)}
  site_dir <- paste0(mydir, "/", gsub("\\_.*", "", name))
  if (!dir.exists(site_dir)){
    dir.create(site_dir)}
  country <- country
  countries <- paste0(country, collapse = "_")
  tf <- tf
  age_file <- age_file
  #read age file loop
  if (!(age_file == F)) {
    age_file <- list.files(pattern = age_file)
  }
  if (age_file == F) {
    age_file
  } else if (grepl(".csv", age_file)) {
    age_file <- read.csv(age_file)
  } else{
    age_file <- read.delim(age_file)
  }
  # read data
  if(grepl(".csv", testate_data)){testate_data <- read.csv(testate_data)}
  else{testate_data <- read.csv(paste0(testate_data, ".csv"))}
  #add depth
  if(!(age_file == F)) {depth <- seq(
    from = depth_start,
    by = depth_int,
    length.out = nrow(age_file))
  age_file$depth <- depth}
  
  depth <- seq(
    from = depth_start,
    by = depth_int,
    length.out = nrow(testate_data))
  
  csv <- paste0(site_dir, "/", name, "_", tf, "_", countries)
  csv_f <- paste0(site_dir, "/", name, "_", tf)
  data(eu, envir = environment())
  data(na, envir = environment())
  boot_size <- boot_size
  
  if (tf == "eu") {
    EuroTF <- eu
  } else{
    EuroTF <-  na
  }
  
  
  
  #To filter by country, change country from F and write countries wanted
  EuroTF <- if (country == F) {
    EuroTF
  }
  else{
    filter(EuroTF, COUNTRY %in% c(country))
  }
  
  #Extract species from tf file
  Spec <- EuroTF[, 9:55]
  
  #Remove sum 0 columns
  Spec <- Spec[, which(colSums(Spec) != 0)]
  
  #Extract WTD
  WT <- EuroTF$WTD
  
  #Run Model
  EuroTF_model <- WA(Spec, WT, tolDW = T)
  
  EuroTF_model.cv <- crossval(EuroTF_model, cv.method = "loo")
  Sys.sleep(20)
  #Check performance
  print(performance(EuroTF_model.cv))
  
  model_stats <-
    as.data.frame(do.call(rbind, performance(EuroTF_model.cv)))
  model_stats$site <- site_nme 
  
  #Run model on data
  EuroTF_recon <-
    predict(EuroTF_model.cv,
            testate_data,
            sse = TRUE,
            nboot = boot_size)
  
  #Create data
  recon <<- data.frame(EuroTF_recon$fit.boot) %>%
    mutate(depth = depth,
           type = "reconstuction")
  
  error <<- data.frame(EuroTF_recon$SEP.boot) %>%
    mutate(depth = depth,
           type = "error")
  
  zscores <<- standardize(data.frame(EuroTF_recon$fit.boot)) %>%
    mutate(depth = depth,
           type = "zscores")
  
  # if loop for age wtd plots
  if (!(age_file == F)) {
    recon <- left_join(recon, age_file)
    
    error <- left_join(error, age_file)
    
    zscores <- left_join(zscores, age_file)
  }
  
  
  # create age or depth plots
  scat_plot <- function(data, var_x, var_y, x, y) {
    lab_y <- max(data$WA.inv.tol)
    data <- data
    var_x <- var_x
    var_y <- var_y
    x <- x
    y <- y
    ggscatter(
      data = data,
      x = var_x,
      y = var_y,
      add = "reg.line",
      conf.int = T,
      add.params = list(color = "blue",
                        fill = "lightgray")
    ) +
      xlab(x) +
      ylab(y) +
      stat_cor(aes(
        label = paste(..r.label.., ..rr.label.., ..p.label.., sep = "~`,`~")
      ),
      label.y = lab_y)
  }
  
  all_data <- rbind(recon, error, zscores)
  all_data$site <- site_nme
  all_data <<- all_data
  
  if (age_file == F) {
    recon_plot <<- scat_plot(
      recon,
      var_x = "depth",
      var_y = "WA.inv.tol",
      x = "depth(cm)",
      y = "WTD (cm)"
    ) + ggtitle(paste0(name, " WTD Reconstruction"))
    
    z_plot <<- scat_plot(
      zscores,
      var_x = "depth",
      var_y = "WA.inv.tol",
      x = "depth(cm)",
      y = "z_score"
    ) + ggtitle(paste0(name, " WTD Deviation"))
  }
  
  else {
    recon_plot <<- scat_plot(
      recon,
      var_x = "mean",
      var_y = "WA.inv.tol",
      x = "Age(y/bp)",
      y = "WTD (cm)"
    ) +
      ggtitle(paste0(name, " WTD Reconstruction"))
    
    z_plot <<- scat_plot(
      zscores,
      var_x = "mean",
      var_y = "WA.inv.tol",
      x = "Age(y/bp)",
      y = "z_score"
    ) + ggtitle(paste0(name, " WTD zscores"))
    
  }
  print(recon_plot)
  print(z_plot)
  
  # save
  if (save == T & country == F) {
    
    write_csv(all_data, file = paste0(csv_f, "_reconstruction.csv"))
    
    write_csv(model_stats, file = paste0(csv_f, "_model_performance.csv"))
    
    ggsave(
      filename = paste0(csv_f, "_z_plot.png"),
      plot = z_plot,
      units = "cm",
      height = 20,
      width = 20
    )
    
    ggsave(
      filename = paste0(csv_f, "_recon_plot.png"),
      plot = recon_plot,
      units = "cm",
      height = 20,
      width = 20
    )
  } else if (save == T){
    
    write.csv(all_data, file = paste0(csv, "_reconstruction.csv"))
    
    write.csv(model_stats, file = paste0(csv, "_model_performance.csv"))
    
    ggsave(
      filename = paste0(csv, "_z_plot.png"),
      plot = z_plot,
      units = "cm",
      height = 20,
      width = 20
    )
    
    ggsave(
      filename = paste0(csv, "_recon_plot.png"),
      plot = recon_plot,
      units = "cm",
      height = 20,
      width = 20)
  }
}





