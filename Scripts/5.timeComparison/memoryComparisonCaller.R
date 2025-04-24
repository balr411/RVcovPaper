#Script that does the memory comparison between two-stage and full approaches
#for chromosome 2, mean sphered cell volume, and pLOF + missense(broad) mask
#using GNU time

library("stringr")

#First read in all of the genes 
genes <- read.table("../Genes/c2.genes", header = F)
names(genes) <- c("Gene", "Start", "End")
chr <- "c2"

test <- 139974
ref <- 1000
traitname <- "meanSpheredCellVolume"
trait <- traitname

group_file_prefix <- c("BROADMISSENSE")

#For each desired group file, read in first line to see if it actually has any variants
empty_gene <- c()
non_empty_gene <- c()
prefix <- group_file_prefix[1]
for(gene in genes[[1]]){
  group_file_path <- paste("../SKAT_weightedBurden/groupFilesNew", prefix, chr, test, traitname, "", sep = "/")
  group_file <- paste(c(group_file_path, paste(chr, test, traitname, gene, prefix, "group.file", sep = ".")), collapse = "")
  file_current <- readLines(group_file, n = 1)
  l <- unlist(strsplit(file_current, split = "\t"))
  
  if(length(l) == 1){
    empty_gene <- c(empty_gene, gene)
  }else{
    non_empty_gene <- c(non_empty_gene, gene)
  }
}

length(non_empty_gene) #1,157
genes <- genes[genes$Gene %in% non_empty_gene,]

#First create the directories to write the covariances and raremetal outputs to
iters <- c(1:5)
for(iter in iters){
  print(iter)
  for(gene_curr in genes[[1]]){
    dir_to_make <- str_glue("memoryComparison/twoStageCovariance/{iter}/c2/{gene_curr}/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
    dir_to_make <- str_glue("memoryComparison/twoStageRaremetal/{iter}/c2/{gene_curr}/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
    dir_to_make <- str_glue("memoryComparison/fullCovariance/{iter}/c2/{gene_curr}/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
    dir_to_make <- str_glue("memoryComparison/fullRaremetal/{iter}/c2/{gene_curr}/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
    dir_to_make <- str_glue("memoryComparison/timeOutputs/twoStage/{iter}/c2/{gene_curr}/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
    dir_to_make <- str_glue("memoryComparison/timeOutputs/trad/{iter}/c2/{gene_curr}/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
  }
}

#Now iterate through and call 
iters <- 1:5
for(iter in iters){
  for(i in 1:length(non_empty_gene)){
    if(i %% 100 == 0){
      print(iter)
      print(i)
    }
    
    idx <- which(genes[[1]] == non_empty_gene[i])
    
    gene <- genes[idx,1]
    gene_start <- as.numeric(genes[idx,2])
    gene_end <- as.numeric(genes[idx,3])
    
    
    score_stat_file <- str_glue("../SKAT_weightedBurden/rmw/c2/139974/meanSpheredCellVolume/c2.139974.{gene}.meanSpheredCellVolume.singlevar.score.txt.gz")
    vcf_file <- str_glue("../SKAT_weightedBurden/refVCFs/{chr}/{ref}/{trait}/{chr}.{ref}.{trait}.{gene}.ref.vcf.gz")
    group_file <- str_glue("../SKAT_weightedBurden/groupFilesNew/BROADMISSENSE/c2/139974/meanSpheredCellVolume/c2.139974.meanSpheredCellVolume.{gene}.BROADMISSENSE.group.file")
    altCovariancePath <- str_glue("memoryComparison/twoStageCovariance/{iter}/c2/{gene}/")
    altRaremetalPath <- str_glue("memoryComparison/twoStageRaremetal/{iter}/c2/{gene}/")
    outputFile <- str_glue("memoryComparison/timeOutputs/twoStage/{iter}/c2/{gene}/time.txt")
    
    #Now compute and time the covariances 
    cm <- str_glue("{{ /usr/bin/time -v Rscript Scripts/5.timeComparison/memoryComparison.R --twoStage TRUE --scoreStat {score_stat_file} --vcf {vcf_file} --groupFile {group_file} --altCovariancePath {altCovariancePath} --altRaremetalPath {altRaremetalPath} --gene {gene} ; }} 2> {outputFile}")
    system(cm)
    
    altCovariancePath <- str_glue("memoryComparison/fullCovariance/{iter}/c2/{gene}/")
    altRaremetalPath <- str_glue("memoryComparison/fullRaremetal/{iter}/c2/{gene}/")
    outputFile <- str_glue("memoryComparison/timeOutputs/trad/{iter}/c2/{gene}/time.txt")
    
    cm <- str_glue("{{ /usr/bin/time -v Rscript Scripts/5.timeComparison/memoryComparison.R --twoStage FALSE --scoreStat {score_stat_file} --vcf {vcf_file} --groupFile {group_file} --altCovariancePath {altCovariancePath} --altRaremetalPath {altRaremetalPath} --gene {gene} ; }} 2> {outputFile}")
    system(cm)
  }
}




