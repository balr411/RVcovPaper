#Script that reads in the results of the memory comparison and makes some plots

library("stringr")
library("data.table")
library("ggplot2")

#First read in all of the genes 
genes <- read.table("../Genes/c2.genes", header = F)
names(genes) <- c("Gene", "Start", "End")
chr <- "c2"

test <- 139974
ref <- 9999
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

#Now read in all of the results
two_stage_mem_df <- data.frame(Gene = genes$Gene, 
                               iter1 = rep(NA, nrow(genes)), 
                               iter2 = rep(NA, nrow(genes)),
                               iter3 = rep(NA, nrow(genes)), 
                               iter4 = rep(NA, nrow(genes)),
                               iter5 = rep(NA, nrow(genes)))

two_stage_user_df <- data.frame(Gene = genes$Gene, 
                                iter1 = rep(NA, nrow(genes)), 
                                iter2 = rep(NA, nrow(genes)),
                                iter3 = rep(NA, nrow(genes)), 
                                iter4 = rep(NA, nrow(genes)),
                                iter5 = rep(NA, nrow(genes)))

two_stage_system_df <- data.frame(Gene = genes$Gene, 
                                  iter1 = rep(NA, nrow(genes)), 
                                  iter2 = rep(NA, nrow(genes)),
                                  iter3 = rep(NA, nrow(genes)), 
                                  iter4 = rep(NA, nrow(genes)),
                                  iter5 = rep(NA, nrow(genes)))

trad_mem_df <- data.frame(Gene = genes$Gene, 
                          iter1 = rep(NA, nrow(genes)), 
                          iter2 = rep(NA, nrow(genes)),
                          iter3 = rep(NA, nrow(genes)), 
                          iter4 = rep(NA, nrow(genes)),
                          iter5 = rep(NA, nrow(genes)))

trad_user_df <- data.frame(Gene = genes$Gene, 
                           iter1 = rep(NA, nrow(genes)), 
                           iter2 = rep(NA, nrow(genes)),
                           iter3 = rep(NA, nrow(genes)), 
                           iter4 = rep(NA, nrow(genes)),
                           iter5 = rep(NA, nrow(genes)))

trad_system_df <- data.frame(Gene = genes$Gene, 
                             iter1 = rep(NA, nrow(genes)), 
                             iter2 = rep(NA, nrow(genes)),
                             iter3 = rep(NA, nrow(genes)), 
                             iter4 = rep(NA, nrow(genes)),
                             iter5 = rep(NA, nrow(genes)))

for(iter in 1:5){
  for(i in 1:nrow(genes)){
    #Keep track of how far along
    if(i %% 100 == 0){
      print(paste0("Iter = ", iter, " i = ", i))
    }
    
    gene <- genes[[1]][i] 
    time_file_twoStage <- str_glue("memoryComparison/timeOutputs/twoStage/{iter}/9999/c2/{gene}/time.txt")
    df_curr <- fread(time_file_twoStage)
    
    #Add to data frame
    two_stage_user_df[i, (iter + 1)] <- as.numeric(unlist(strsplit(unlist(df_curr[2,2]), split = ":", fixed = TRUE))[2])
    two_stage_system_df[i, (iter + 1)] <- as.numeric(unlist(strsplit(unlist(df_curr[3,2]), split = ":", fixed = TRUE))[2])
    two_stage_mem_df[i, (iter + 1)] <- as.numeric(unlist(strsplit(unlist(df_curr[10,2]), split = ":", fixed = TRUE))[2])
    
    time_file_trad <- str_glue("memoryComparison/timeOutputs/trad/{iter}/9999/c2/{gene}/time.txt")
    df_curr <- fread(time_file_trad)
    
    #Add to data frame
    trad_user_df[i, (iter + 1)] <- as.numeric(unlist(strsplit(unlist(df_curr[2,2]), split = ":", fixed = TRUE))[2])
    trad_system_df[i, (iter + 1)] <- as.numeric(unlist(strsplit(unlist(df_curr[3,2]), split = ":", fixed = TRUE))[2])
    trad_mem_df[i, (iter + 1)] <- as.numeric(unlist(strsplit(unlist(df_curr[10,2]), split = ":", fixed = TRUE))[2])
    
  }
}


#Look at summaries
df_comb <- data.frame(two_stage_user_median = apply(two_stage_user_df[,2:6], 1, median),
                      two_stage_user_mean = rowMeans(two_stage_user_df[,2:6]),
                      two_stage_system_median = apply(two_stage_system_df[,2:6], 1, median),
                      two_stage_system_mean = rowMeans(two_stage_system_df[,2:6]),
                      two_stage_mem_median =  apply(two_stage_mem_df[,2:6], 1, median),
                      two_stage_mem_mean = rowMeans(two_stage_mem_df[,2:6]),
                      two_stage_elapsed_median = apply(two_stage_system_df[,2:6] + two_stage_user_df[,2:6], 1, median),
                      two_stage_elapsed_mean = rowMeans(two_stage_system_df[,2:6] + two_stage_user_df[,2:6]),
                      trad_user_median = apply(trad_user_df[,2:6], 1, median),
                      trad_user_mean = rowMeans(trad_user_df[,2:6]),
                      trad_system_median = apply(trad_system_df[,2:6], 1, median),
                      trad_system_mean = rowMeans(trad_system_df[,2:6]),
                      trad_mem_median =  apply(trad_mem_df[,2:6], 1, median),
                      trad_mem_mean = rowMeans(trad_mem_df[,2:6]),
                      trad_elapsed_median = apply(trad_system_df[,2:6] + trad_user_df[,2:6], 1, median),
                      trad_elapsed_mean = rowMeans(trad_system_df[,2:6] + trad_user_df[,2:6]))

#Memory summaries in GB:
summary(df_comb$two_stage_mem_median)/1000000
summary(df_comb$trad_mem_median)/1000000
