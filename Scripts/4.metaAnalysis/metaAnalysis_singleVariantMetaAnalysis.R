#Script that calls raremetal on each gene to get then single-variant summary statistics
#from the single-variant meta-analysis

library("stringr")

#Read in all gene files 
dat <- read.table("../Genes/c1.genes", header = F)
dat[[4]] <- "c1"

for(i in 2:22){
  chr_current <- paste0("c", i)
  dat_current <- read.table(paste0("../Genes/", chr_current, ".genes"), header = F)
  dat_current[[4]] <- chr_current
  dat <- rbind(dat, dat_current)
}

#First create the directories 
chrs <- paste0("c", 1:22)

for(chr in chrs){
  dir_to_make <- str_glue("singleVariant/{chr}")
  if(!dir.exists(dir_to_make)){
    dir.create(dir_to_make)
  }
}

for(chr in chrs){
  dir_to_make <- str_glue("singleVariant/{chr}/6998")
  if(!dir.exists(dir_to_make)){
    dir.create(dir_to_make)
  }
}

for(chr in chrs){
  dir_to_make <- str_glue("singleVariant/{chr}/6998/meanSpheredCellVolume")
  if(!dir.exists(dir_to_make)){
    dir.create(dir_to_make)
  }
}



test <- 6998
trait <- "meanSpheredCellVolume"
for(i in 1:dim(dat)[1]){
  if(i %% 100 == 0){
    print(i)
  }
  
  chr <- dat[i,4]
  gene <- dat[i,1]
  
  summary_files <- str_glue("raremetalFilesNew/summary/{chr}/{test}/{trait}/{chr}.{test}.{trait}.{gene}.summaryFiles")
  full_prefix <- str_glue("singleVariant/{chr}/{test}/{trait}/{chr}.{test}.{trait}.{gene}")
  
  cm <- str_glue("/net/snowwhite/home/welchr/projects/raremetal/build/release/bin/raremetal --summaryFiles {summary_files} --prefix {full_prefix}")
  system(cm)
  
  
}


