#This script creates a summaryFile for each gene containing the file names of the 
#20 studies to use in raremetal

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

#Make sure directories are created
chrs <- paste0("c", 1:22)

for(chr in chrs){
  dir_to_make <- str_glue("raremetalFilesNew/summary/{chr}/")
  if(!dir.exists(dir_to_make)){
    dir.create(dir_to_make)
  }
}

for(chr in chrs){
  dir_to_make <- str_glue("raremetalFilesNew/summary/{chr}/6998")
  if(!dir.exists(dir_to_make)){
    dir.create(dir_to_make)
  }
}

for(chr in chrs){
  dir_to_make <- str_glue("raremetalFilesNew/summary/{chr}/6998/meanSpheredCellVolume")
  if(!dir.exists(dir_to_make)){
    dir.create(dir_to_make)
  }
}

#Now all directories should be made so create the files in the form "raremetalFilesNew/summary/{chr}/{test}/{trait}/{chr}.{test}.{trait}.{gene}.summaryFiles"
study_list <- function(s, chr, gene){
  str_glue("rmw/{chr}/{s}/6998/meanSpheredCellVolume/{chr}.{s}.6998.{gene}.meanSpheredCellVolume.singlevar.score.txt.gz")  
}


test <- 6998
trait <- "meanSpheredCellVolume"
for(i in 1:dim(dat)[1]){
  if(i %% 100 == 0){
    print(i)
  }
  chr <- dat[i,4]
  gene <- dat[i,1]
  
  file_curr <- str_glue("raremetalFilesNew/summary/{chr}/{test}/{trait}/{chr}.{test}.{trait}.{gene}.summaryFiles")
  
  studies <- paste0("S", 1:20)
  to_write <- unname(sapply(studies, function(x) study_list(x, chr, gene)))
  to_write <- to_write[file.exists(to_write)] 
  
  fileConn <- file(file_curr)
  writeLines(to_write, fileConn)
  close(fileConn)
}



