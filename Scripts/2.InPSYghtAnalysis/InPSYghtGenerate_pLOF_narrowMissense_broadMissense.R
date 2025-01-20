library("argparser")
library("data.table")
library("dplyr")
library("stringr")

p <- arg_parser("Program that takes set of test variants, reference variants, test allele frequency file, 
                and VEP annotation to produce group file containing pLOF variants, pLOF + missense(narrow) 
                and pLOF + missense(broad)") 
p <- add_argument(p, "--allelefreq", help = "Allele frequency file in form CHROM,POS,REF,ALT,AF,AC,PVAL")
p <- add_argument(p, "--anno", help = "Annotation file to create masks")
p <- add_argument(p, "--multiAllelic", help = "Name of file containing multi-allelic InPSYght variants for that chromosome")
p <- add_argument(p, "--InPSYghtHWE", help = "Name of file containing variants failing HWE for InPSYght")
p <- add_argument(p, "--gene", help = "Name of gene")
p <- add_argument(p, "--geneCoord", help = "Name of file containing gene coordinates")

p <- add_argument(p, "--groupStats", help = "File to output number of variants and allele count in each group")

p <- add_argument(p, "--outputPTV", help = "File name to output group file containing only PTVs")
p <- add_argument(p, "--outputNarrowMissense", help = "File name to output group file containing only pLOF + missense(narrow) variants")
p <- add_argument(p, "--outputBroadMissense", help = "File name to output group file containing only pLOF + missense(broad) variants")

argv <- parse_args(p)

allele_freq <- read.table(argv$allelefreq, header = FALSE)
names(allele_freq) <- c("CHROM", "POS", "REF", "ALT", "AF", "AC", "PVAL")
anno <- read.table(argv$anno, header = TRUE) 

#Delete multi-alleleic variants 
allele_freq <- anti_join(allele_freq, allele_freq[duplicated(allele_freq[,2]),], by = "POS")

#Delete multi-allelic variants from reference set
ref_dup <- fread(argv$multiAllelic, header = F, data.table = F)[[2]]

idx <- which(allele_freq$POS %in% ref_dup)
if(length(idx) > 0){
  allele_freq <- allele_freq[-idx,]
}

#Now remove variants failing HWE in InPSYght
ref_hwe <- fread(argv$InPSYghtHWE, header = F, data.table = F)[[2]]

idx <- which(allele_freq$POS %in% ref_hwe)
if(length(idx) > 0){
  allele_freq <- allele_freq[-idx,]
}

#Delete variants with AF > 0.01 and AF = 0
allele_freq <- allele_freq[allele_freq[,5] > 0 & allele_freq[,5] <= 0.01,]

#Function converting chr, pos, ref, alt to one string
pos_to_loc <- function(chr, pos, ref, alt){
  return(paste(chr, pos, ref, alt, sep = ":"))
}

allele_freq$UploadedVariation <- apply(allele_freq, 1, function(x) pos_to_loc(as.numeric(x[1]), as.numeric(x[2]), x[3], x[4]))


#Create normal group file containing only PTVs with MAF < 0.01 - note for single gene this could give back an empty data frame
allele_freq_ptv <- allele_freq[allele_freq[,5] <= 0.01,]

#Create normal group file containing only missense variants with MAF <0.01
allele_freq_narrowMissense <- allele_freq[allele_freq[,5] <= 0.01,]

#Create normal group file containing only non-synomynous variants with MAF < 0.01 
allele_freq_broadMissense <- allele_freq[allele_freq[,5] <= 0.01,]

#Match variants in annotation file to those in allele frequency file - note ptv allele freq has all variants in test set w/ AF < 0.01 right now
idx_common_anno <- which(anno[[1]] %in% allele_freq_ptv$UploadedVariation)
anno <- anno[idx_common_anno,]

#Find pLOF variants and missense variants - note that every variant here should be a missense variant or pLOF variant
cons <- unique(anno[[2]])
imp_vars <- c()
imp_vars_missense <- c()
bad_cons <- c("stop_gained", "splice_acceptor_variant", "stop_lost", "splice_donor_variant", "frameshift_variant", "start_lost")
if(length(cons) > 0){
  for(i in 1:length(cons)){
    temp <- unlist(strsplit(cons[i], split = ","))
    if(sum(bad_cons %in% temp) > 0){
      imp_vars <- c(imp_vars, cons[i])
    }
    
    if((sum(bad_cons %in% temp) == 0) & (sum(temp == "missense_variant") > 0)){ 
      imp_vars_missense <- c(imp_vars_missense, cons[i])
    }
  }
}

idx_annoPTV <- which(anno[[2]] %in% imp_vars)
idx_anno_narrowMissense <- c(which(anno[[2]] %in% imp_vars), which((anno[[2]] %in% imp_vars_missense) & (anno$Sum_Algs == 1)))
idx_anno_narrowMissense <- idx_anno_narrowMissense[order(idx_anno_narrowMissense)]
idx_anno_broadMissense <- c(which(anno[[2]] %in% imp_vars), which(anno[[2]] %in% imp_vars_missense))
idx_anno_broadMissense <- idx_anno_broadMissense[order(idx_anno_broadMissense)]

anno_ptv <- anno[idx_annoPTV,]
anno_narrowMissense <- anno[idx_anno_narrowMissense,]
anno_broadMissense <- anno[idx_anno_broadMissense,]

#Table that contains number of variants for each group
genes_unique <- argv$gene
mask_list <- as.list(genes_unique)

group_table <- data.frame(Genes = genes_unique, ptv_vars = 0, ptv_ac = 0, 
                          narrowMissense_vars = 0, narrowMissense_ac = 0, 
                          broadMissense_vars = 0, broadMissene_ac = 0)

###################################################################################
#PTVs 
to_add <- unique(anno_ptv[anno_ptv$SYMBOL == genes_unique,][[1]])
if(length(to_add) > 0){
  if(sum(is.na(to_add)) == 0){
    mask_list[[1]] <- c(mask_list[[1]], to_add)
  }
}

group_table[1,2] <- length(mask_list[[1]]) - 1


to_write <- vector(length = 1)
to_write[1] <- paste(mask_list[[1]], collapse = "\t")

#If there are any variants, put their allele counts in the proper column, otherwise they stay as 0
if(length(mask_list[[1]]) > 1){
  idx_temp <- which(allele_freq_ptv$UploadedVariation %in% mask_list[[1]])
  group_table[1,3] <- group_table[1,3] + sum(allele_freq_ptv[idx_temp,6])
}

fileConn <- file(argv$outputPTV)
writeLines(to_write, fileConn)
close(fileConn)

####################################################################################################
#Create pLOF + missense(narrow) variant masks
mask_list <- as.list(genes_unique)

to_add <- unique(anno_narrowMissense[anno_narrowMissense$SYMBOL == genes_unique,][[1]]) 

if(length(to_add) > 0){
  if(sum(is.na(to_add)) == 0){
    mask_list[[1]] <- c(mask_list[[1]], to_add)
  }
}

group_table[1,4] <- length(mask_list[[1]]) - 1

to_write <- vector(length = 1)
to_write[1] <- paste(mask_list[[1]], collapse = "\t")

if(length(mask_list[[1]]) > 1){
  idx_temp <- which(allele_freq_narrowMissense$UploadedVariation %in% mask_list[[1]])
  group_table[1,5] <- group_table[1,5] + sum(allele_freq_narrowMissense[idx_temp,6])
}

fileConn <- file(argv$outputNarrowMissense)
writeLines(to_write, fileConn)
close(fileConn)

############################################################################################
#Create non-synomynous variant masks
mask_list <- as.list(genes_unique)

to_add <- unique(anno_broadMissense[anno_broadMissense$SYMBOL == genes_unique,][[1]]) 

if(length(to_add) > 0){
  if(sum(is.na(to_add)) == 0){
    mask_list[[1]] <- c(mask_list[[1]], to_add)
  }
}

group_table[1,6] <- length(mask_list[[1]]) - 1

to_write <- vector(length = 1)
to_write[1] <- paste(mask_list[[1]], collapse = "\t")

if(length(mask_list[[1]]) > 1){
  idx_temp <- which(allele_freq_broadMissense$UploadedVariation %in% mask_list[[1]])
  group_table[1,7] <- group_table[1,7] + sum(allele_freq_broadMissense[idx_temp,6])
}

fileConn <- file(argv$outputBroadMissense)
writeLines(to_write, fileConn)
close(fileConn)

#Write out group stat files
write.table(group_table, argv$groupStats)








