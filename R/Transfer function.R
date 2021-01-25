trans_fun <- function(testate_data = "sq_tests",
                      tf = "eu",
                      country = "scotland",
                      save = T,
                      age_file = F,
                      depth_start = 1,
                      depth_int = 1) {
  require(rioja)
  require(tidyverse)
  require(ggpubr)
  require(vegan)
  require(effectsize)
  mydir <- paste0("tf_runs/")
  if (!dir.exists(mydir))
    dir.create(mydir)
  name <- testate_data
  country <- country
  tf <- tf
  age_file <- age_file
  #read age file loop
  if (age_file == F) {
    age_file
  }
  else{
    age_file <- list.files(pattern = age_file)
  }
  if (age_file == F) {
    age_file
  }
  else if (".csv" %in% age_file) {
    age_file <- read.csv(age_file)
  }
  else {
    age_file <- read.delim(age_file)
  }

  # add depth to testates if it doesn't exists
  testate_data <- read.csv(paste0(testate_data, ".csv"))
  if ("depth" %in% colnames(testate_data)) {
    depth <- testate_data$depth
  }
  else{
    depth <- seq(
      from = depth_start,
      by = depth_int,
      length.out = nrow(testate_data)
    )
  }

  csv <- paste0(mydir, "/", name, "_", tf, "_", country)
  data(eu, envir = environment())
  data(na, envir = environment())

  if (tf == "eu") {
    EuroTF <<- eu
  } else{
    EuroTF <<-  na
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
  Spec <-  Spec[, colSums(Spec != 0) > 0]

  #Extract WTD
  WT <- EuroTF$WTD

  #Run Model
  EuroTF_model <- WA(Spec, WT, tolDW = T)

  EuroTF_model.cv <- crossval(EuroTF_model, cv.method = "loo")
  #Check performance
  print(performance(EuroTF_model.cv))

  model_stats <-
    as.data.frame(do.call(rbind, performance(EuroTF_model.cv)))


  if (save == T) {
    write.csv(model_stats, file = paste0(csv, "_performance.csv"))
  }

  #Run model on data
  EuroTF_recon <-
    predict(EuroTF_model.cv,
            testate_data,
            sse = TRUE,
            nboot = 1000)

  #Create data
  recon <- data.frame(EuroTF_recon$fit.boot) %>%
    mutate(depth = depth,
           type = "reconstuction")

  error <- data.frame(EuroTF_recon$SEP.boot) %>%
    mutate(depth = depth,
           type = "error")

  zscores <- standardize(data.frame(EuroTF_recon$fit.boot)) %>%
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

  all_data <<- rbind(recon, error, zscores)

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
      x = "age(y/bp)",
      y = "WTD (cm)"
    ) +
      ggtitle(paste0(name, " WTD Reconstruction"))

    z_plot <<- scat_plot(
      zscores,
      var_x = "mean",
      var_y = "WA.inv.tol",
      x = "depth(cm)",
      y = "z_score"
    ) + ggtitle(paste0(name, " WTD zscores"))

  }
  print(recon_plot)
  print(z_plot)

  # save
  if (save == T) {
    write.csv(all_data, file = paste0(csv, "_reconstruction.csv"))


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
      width = 20
    )
  }
}

