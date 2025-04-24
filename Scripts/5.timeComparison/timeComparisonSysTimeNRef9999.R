#Script that does the time comparison between two-stage and full approaches
#for chromosome 2, mean sphered cell volume, and pLOF + missense(broad) mask

library("stringr")
library("RVcov")

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

#First create the directories to write the covariances and raremetal outputs to
iters <- c(1:5)

for(iter in iters){
  print(iter)
  for(gene_curr in genes[[1]]){
    dir_to_make <- str_glue("twoStageCovariance/9999/{iter}/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
    dir_to_make <- str_glue("twoStageRaremetal/9999/{iter}/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
    dir_to_make <- str_glue("fullCovariance/9999/{iter}/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
    dir_to_make <- str_glue("fullRaremetal/9999/{iter}/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
  }
}

for(iter in iters){
  print(iter)
  for(gene_curr in genes[[1]]){
    dir_to_make <- str_glue("twoStageCovariance/9999/{iter}/c2/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
    dir_to_make <- str_glue("twoStageRaremetal/9999/{iter}/c2/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
    dir_to_make <- str_glue("fullCovariance/9999/{iter}/c2/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
    dir_to_make <- str_glue("fullRaremetal/9999/{iter}/c2/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
  }
}

for(iter in iters){
  print(iter)
  for(gene_curr in genes[[1]]){
    dir_to_make <- str_glue("twoStageCovariance/9999/{iter}/c2/{gene_curr}/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
    dir_to_make <- str_glue("twoStageRaremetal/9999/{iter}/c2/{gene_curr}/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
    dir_to_make <- str_glue("fullCovariance/9999/{iter}/c2/{gene_curr}/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
    dir_to_make <- str_glue("fullRaremetal/9999/{iter}/c2/{gene_curr}/")
    if(!dir.exists(dir_to_make)){
      dir.create(dir_to_make)
    }
    
  }
}

two_stage_user_list <- list()
two_stage_elapsed_list <- list()
trad_user_list <- list()
trad_elapsed_list <- list()

for(iter in iters){
  two_stage_user <- vector(length = length(non_empty_gene))
  two_stage_elapsed <- vector(length = length(non_empty_gene))
  trad_user <- vector(length = length(non_empty_gene))
  trad_elapsed <- vector(length = length(non_empty_gene))
  
  for(i in 1:length(non_empty_gene)){
    if(i %% 10 == 0){
      print(iter)
      print(i)
      print(mean(two_stage_elapsed[1:i]))
      print(mean(trad_elapsed[1:i]))
    }
    
    idx <- which(genes[[1]] == non_empty_gene[i])
    
    gene <- genes[idx,1]
    gene_start <- NULL
    gene_end <- NULL
    
    
    score_stat_file <- str_glue("../SKAT_weightedBurden/rmw/c2/139974/meanSpheredCellVolume/c2.139974.{gene}.meanSpheredCellVolume.singlevar.score.txt.gz")
    vcf_file <- str_glue("../SKAT_weightedBurden/refVCFs/{chr}/{ref}/{trait}/{chr}.{ref}.{trait}.{gene}.ref.vcf.gz")
    chr_num <- 2
    burden <- wburden <- SKAT <- TRUE
    anno_file <- NULL
    anno <- NULL
    two_stage_threshold <- 3
    group_file <- c(broadMissense = str_glue("../SKAT_weightedBurden/groupFilesNew/BROADMISSENSE/c2/139974/meanSpheredCellVolume/c2.139974.meanSpheredCellVolume.{gene}.BROADMISSENSE.group.file"))
    pLOF <- pLOF_narrowMissense <- pLOF_broadMissense <- FALSE
    altGroupFilePath <- NULL 
    altCovariancePath <- str_glue("twoStageCovariance/9999/{iter}/c2/{gene}/")
    altRaremetalPath <- str_glue("twoStageRaremetal/9999/{iter}/c2/{gene}/")
    altRaremetalName <- NULL 
    mafThreshold <- 0.01
    hwe <- 0.000001
    
    
    #Now compute and time the covariances 
    time_two_stage <- system.time(agg_test(score_stat_file, vcf_file, chr = chr_num, burden, wburden,
                                           SKAT, anno_file,
                                           anno, two_stage = TRUE, two_stage_threshold,
                                           group_file, pLOF,
                                           pLOF_narrowMissense, pLOF_broadMissense,
                                           altGroupFilePath, altCovariancePath,
                                           altRaremetalPath, altRaremetalName,
                                           mafThreshold, gene,
                                           gene_start, gene_end, hwe))
    
    
    two_stage_user[i] <- time_two_stage[1] + time_two_stage[4]
    two_stage_elapsed[i] <- time_two_stage[3]
    
    altCovariancePath <- str_glue("fullCovariance/9999/{iter}/c2/{gene}/")
    altRaremetalPath <- str_glue("fullRaremetal/9999/{iter}/c2/{gene}/")
    
    time_trad <- system.time(agg_test(score_stat_file, vcf_file, chr = chr_num, burden, wburden,
                                      SKAT, anno_file,
                                      anno, two_stage = FALSE, two_stage_threshold,
                                      group_file, pLOF,
                                      pLOF_narrowMissense, pLOF_broadMissense,
                                      altGroupFilePath, altCovariancePath,
                                      altRaremetalPath, altRaremetalName,
                                      mafThreshold, gene,
                                      gene_start, gene_end, hwe))
    
    trad_user[i] <- time_trad[1] + time_trad[4]
    trad_elapsed[i] <- time_trad[3]
    
  }
  
  two_stage_user_list[[iter]] <- two_stage_user
  two_stage_elapsed_list[[iter]] <- two_stage_elapsed
  trad_user_list[[iter]] <- trad_user
  trad_elapsed_list[[iter]] <- trad_elapsed
  
}

#Write the results
#Write out the other runs as well 
#Make a data frame to write out 
df_run1 <- data.frame(gene = non_empty_gene,
                      trad_user = trad_user_list[[1]], 
                      trad_elapsed = trad_elapsed_list[[1]], 
                      two_stage_user = two_stage_user_list[[1]],
                      two_stage_elapsed = two_stage_elapsed_list[[1]])

write.table(df_run1, file = "run1_9999.res", quote = FALSE, row.names = FALSE)

df_run2 <- data.frame(gene = non_empty_gene,
                      trad_user = trad_user_list[[2]], 
                      trad_elapsed = trad_elapsed_list[[2]], 
                      two_stage_user = two_stage_user_list[[2]],
                      two_stage_elapsed = two_stage_elapsed_list[[2]])

write.table(df_run2, file = "run2_9999.res", quote = FALSE, row.names = FALSE)

df_run3 <- data.frame(gene = non_empty_gene,
                      trad_user = trad_user_list[[3]], 
                      trad_elapsed = trad_elapsed_list[[3]], 
                      two_stage_user = two_stage_user_list[[3]],
                      two_stage_elapsed = two_stage_elapsed_list[[3]])

write.table(df_run3, file = "run3_9999.res", quote = FALSE, row.names = FALSE)

df_run4 <- data.frame(gene = non_empty_gene,
                      trad_user = trad_user_list[[4]], 
                      trad_elapsed = trad_elapsed_list[[4]], 
                      two_stage_user = two_stage_user_list[[4]],
                      two_stage_elapsed = two_stage_elapsed_list[[4]])

write.table(df_run4, file = "run4_9999.res", quote = FALSE, row.names = FALSE)

df_run5 <- data.frame(gene = non_empty_gene,
                      trad_user = trad_user_list[[5]], 
                      trad_elapsed = trad_elapsed_list[[5]], 
                      two_stage_user = two_stage_user_list[[5]],
                      two_stage_elapsed = two_stage_elapsed_list[[5]])

write.table(df_run5, file = "run5_9999.res", quote = FALSE, row.names = FALSE)

#First do overlapping histograms of the median and mean times 
#Make a data frame for each different method
df_two_stage_user <- data.frame(run1 = two_stage_user_list[[1]],
                                run2 = two_stage_user_list[[2]],
                                run3 = two_stage_user_list[[3]],
                                run4 = two_stage_user_list[[4]],
                                run5 = two_stage_user_list[[5]])

df_two_stage_elapsed <- data.frame(run1 = two_stage_elapsed_list[[1]],
                                   run2 = two_stage_elapsed_list[[2]],
                                   run3 = two_stage_elapsed_list[[3]],
                                   run4 = two_stage_elapsed_list[[4]],
                                   run5 = two_stage_elapsed_list[[5]])

df_trad_user <- data.frame(run1 = trad_user_list[[1]],
                           run2 = trad_user_list[[2]],
                           run3 = trad_user_list[[3]],
                           run4 = trad_user_list[[4]],
                           run5 = trad_user_list[[5]])

df_trad_elapsed <- data.frame(run1 = trad_elapsed_list[[1]],
                              run2 = trad_elapsed_list[[2]],
                              run3 = trad_elapsed_list[[3]],
                              run4 = trad_elapsed_list[[4]],
                              run5 = trad_elapsed_list[[5]])

#Now make a data frame of the median and means
df_comb <- data.frame(two_stage_user_median = apply(df_two_stage_user, 1, median),
                      two_stage_user_mean = rowMeans(df_two_stage_user),
                      two_stage_elapsed_median = apply(df_two_stage_elapsed, 1, median),
                      two_stage_elapsed_mean = rowMeans(df_two_stage_elapsed),
                      trad_user_median = apply(df_trad_user, 1, median),
                      trad_user_mean = rowMeans(df_trad_user),
                      trad_elapsed_median = apply(df_trad_elapsed, 1, median),
                      trad_elapsed_mean = rowMeans(df_trad_elapsed))

#Do a table of the summaries for each
summary(df_comb$two_stage_elapsed_median)
summary(df_comb$trad_elapsed_median)
sum(df_comb$two_stage_elapsed_median)
sum(df_comb$trad_elapsed_median)

