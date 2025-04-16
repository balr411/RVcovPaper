#This script takes as input a list of samples (n = 139,974) used for mean sphered cell volume
#during the SKAT and weighted burden analysis and outputs 20 meta analysis study
#sample subsets of size 6998

library("stringr")

samples <- read.table("../SKAT_weightedBurden/testSize/meanSpheredCellVolume/139974.meanSpheredCellVolume.sample.test")[[1]]

n_test <- 6998

for(i in 1:20){
  set.seed(i)
  samples_curr <- sample(samples, n_test)
  samples_curr <- as.character(samples_curr)
  samples <- samples[!(samples %in% samples_curr)]
  
  study_curr <- paste0("S", i)
  file_curr <- str_glue("testSize/meanSpheredCellVolume/{study_curr}.6998.meanSpheredCellVolume.sample.test")
  
  fileConn <- file(file_curr)
  writeLines(samples_curr, fileConn)
  close(fileConn)
}

