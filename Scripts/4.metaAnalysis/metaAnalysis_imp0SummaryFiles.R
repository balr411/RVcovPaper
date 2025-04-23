#Script that calls raremetal on each gene to get then single-variant summary statistics
#from the single-variant meta-analysis

library("stringr")

#Read in all gene files 
dat <- read.table("../Genes/c1.genes.metaAnalysis", header = F)
dat[[4]] <- "c1"

for(i in 2:22){
  chr_current <- paste0("c", i)
  dat_current <- read.table(paste0("../Genes/", chr_current, ".genes.metaAnalysis"), header = F)
  dat_current[[4]] <- chr_current
  dat <- rbind(dat, dat_current)
}

#Make sure directories are created
chrs <- paste0("c", 1:22)

for(chr in chrs){
  dir_to_make <- str_glue("raremetalFilesNew/impCov/{chr}/")
  if(!dir.exists(dir_to_make)){
    dir.create(dir_to_make)
  }
}

for(chr in chrs){
  dir_to_make <- str_glue("raremetalFilesNew/impCov/{chr}/6998")
  if(!dir.exists(dir_to_make)){
    dir.create(dir_to_make)
  }
}

refs <- c(100, 1000, 9999)
for(chr in chrs){
  for(ref in refs){
    dir_to_make <- str_glue("raremetalFilesNew/impCov/{chr}/6998/{ref}")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    } 
  }
}

for(chr in chrs){
  for(ref in refs){
    dir_to_make <- str_glue("raremetalFilesNew/impCov/{chr}/6998/{ref}/meanSpheredCellVolume")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    } 
  }
}

#Now all directories should be made so create the files in the form "raremetalFilesNew/impCov/{chr}/6998/{ref}/meanSpheredCellVolume/{chr}.6998.{ref}.meanSpheredCellVolume.{gene}.imputed0.est_covFiles"
study_list <- function(s, ref, chr, gene){
  str_glue("estCovImp0New/{chr}/{s}/{test}/{ref}/{trait}/{chr}.{s}.{test}.{ref}.{trait}.{gene}.imp0.estimated.cov.gz")  
}


test <- 6998
trait <- "meanSpheredCellVolume"
studies <- paste0("S", 1:20)

for(i in 1:dim(dat)[1]){
  chr <- dat[i,4]
  gene <- dat[i,1]
  
  if(i %% 100 == 0){
    print(i)
  }
  
  for(ref in refs){
    file_curr <- str_glue("raremetalFilesNew/impCov/{chr}/6998/{ref}/meanSpheredCellVolume/{chr}.6998.{ref}.meanSpheredCellVolume.{gene}.imputed0.est_covFiles")
    
    to_write <- unname(sapply(studies, function(x) study_list(x, ref, chr, gene)))
    
    fileConn <- file(file_curr)
    writeLines(to_write, fileConn)
    close(fileConn) 
  }
}











