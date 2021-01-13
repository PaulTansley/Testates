
trans_fun <- function(test_data,
                       region_name = F){
require(rioja)
require(tidyverse)
require(ggpubr)
require(vegan)
require(effectsize)
mydir <- paste0(test_data, "/")
if(!dir.exists(mydir)) dir.create(mydir)
name <- test_data
test_data <- read.csv(paste0(test_data, ".csv"))
EuroTF <- data(eurowood)
region_name <- region_name


#To extract by region change reigion to T and region name to region wanted
EuroTF <-if(region_name == F){EuroTF}
  else{filter(EuroTF, COUNTRY %in% c(region_name))}

#Extract species from euro tf
Spec <- EuroTF[,9:55]

#Extract coulmns that have 0 in every row
Spec <-  Spec[, colSums(Spec != 0) > 0]

#Extract WTD
WT <- EuroTF$WTD

#Run Model
EuroTF_model <- WA(Spec, WT, tolDW = T)

EuroTF_model.cv <- crossval(EuroTF_model, cv.method="loo")
#Check performance
print(performance(EuroTF_model.cv))

perf <- as.data.frame(do.call(rbind,performance(EuroTF_model.cv)))
perf <- perf[-1,]

write.csv(perf, file =paste0(name, "/", name,"_performance.csv"))

#Run model on data
EuroTF_recon <- predict(EuroTF_model.cv, test_data, sse=TRUE, nboot=1000)

#Create data
recon <- data.frame(EuroTF_recon$fit.boot)
error <- EuroTF_recon$SEP.boot
z <- standardize(recon)

write.csv(recon, file = paste0(name,"/",name, "_reconstruction.csv"))

write.csv(error, file = paste0(name,"/",name, "_error.csv"))

write.csv(z, file = paste0(name,"/",name, "_z_scores.csv"))
}

