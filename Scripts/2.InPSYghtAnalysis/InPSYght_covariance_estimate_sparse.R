library("argparser")
library("VariantAnnotation")
library("dplyr")
library("bigmemory")
library("data.table")
library("Matrix")
library("stringr")

p <- arg_parser("Program that computes estimated covariance matrix between variants for a specific gene. Note that this will fail if 
                the VCF of the gene to be read in is too large for R. Also outputs RAREMETAL covFiles file with path to file.")
p <- add_argument(p, "--refallelefreq", help = "Reference panel allele frequency file for gene in form CHROM,POS,REF,ALT,AN,AC")
p <- add_argument(p, "--testallelefreq", help = "Test set allele frequency file for gene in form CHROM,POS,REF,ALT,AF,AC,PVAL")
p <- add_argument(p, "--vcf", help = "Reference set vcf to read in")
p <- add_argument(p, "--testSizeResidualVariance", help = "File containing test sample size in column 1 and residual variance in column 2")
p <- add_argument(p, "--geneCoord", help = "File containing genomic coordinates of genes")
p <- add_argument(p, "--gene", help = "Name of gene")
p <- add_argument(p, "--chr", help = "Chromosome to be worked on")
p <- add_argument(p, "--refsamples", help = "File of reference samples")
p <- add_argument(p, "--refsize", help = "Size of reference panel")
p <- add_argument(p, "--testsize", help = "Size of test set")
p <- add_argument(p, "--multiAllelic", help = "Name of file containing multi-allelic InPSYght variants for that chromosome")
p <- add_argument(p, "--InPSYghtHWE", help = "Name of file containing variants failing HWE for InPSYght")
p <- add_argument(p, "--output", help = "Name of file to output to")


argv <- parse_args(p)
chr <- as.numeric(unlist(strsplit(argv$chr, split = "c"))[2])
n_test <- as.numeric(argv$testsize)
n_ref <- as.numeric(argv$refsize)

allele_freq_reference <- read.table(argv$refallelefreq, header=F)

original_vars <- allele_freq_reference[[2]]

#Ignore the fact that this introduces NAs
allele_freq_reference[[7]] <- allele_freq_reference[[6]]/allele_freq_reference[[5]]
names(allele_freq_reference) <- c("CHROM", "POS", "REF", "ALT", "AN", "AC", "AF")

#read in allele frequency file taken directly from raremetalworker output
allele_freq_test <- read.table(argv$testallelefreq, header=F) 
names(allele_freq_test) <- c("CHROM", "POS", "REF", "ALT", "AF", "AC", "PVAL")

#Remove multiallelic variants in test set 
allele_freq_test <- anti_join(allele_freq_test, allele_freq_test[duplicated(allele_freq_test[,2]),], by = "POS")

#Delete multi-allelic variants from reference set
ref_dup <- fread(argv$multiAllelic, header = F, data.table = F)[[2]]

idx <- which(allele_freq_test$POS %in% ref_dup)
if(length(idx) > 0){
  allele_freq_test <- allele_freq_test[-idx,]
}

#Now remove variants failing HWE in InPSYght
ref_hwe <- fread(argv$InPSYghtHWE, header = F, data.table = F)[[2]]

idx <- which(allele_freq_test$POS %in% ref_hwe)
if(length(idx) > 0){
  allele_freq_test <- allele_freq_test[-idx,]
}

#Now remove all variants from the test set that are multi-allelic in the reference panel 
multi_allelic_ref <- allele_freq_reference[duplicated(allele_freq_reference[,2]),2]
idx <- which(allele_freq_test[[2]] %in% multi_allelic_ref)
if(length(idx) > 0){
  allele_freq_test <- allele_freq_test[-idx,]
}

#remove variants with 0 AF in test set
allele_freq_test <- allele_freq_test[allele_freq_test[[5]] > 0,]

#Now reduce ref allele freq to only include variants found in test set (by SNP not position) 
#And switch alleles that have been switched

#First match on position 
allele_freq_reference <- allele_freq_reference[allele_freq_reference$POS %in% allele_freq_test$POS,]

#Create copy of test allele freq that only contains the variants from the reference panel
allele_freq_test_red_temp <- allele_freq_test[allele_freq_test$POS %in% allele_freq_reference$POS,]

#Find allele switches
idx_switched <- which((allele_freq_reference$REF == allele_freq_test_red_temp$ALT) & (allele_freq_reference$ALT == allele_freq_test_red_temp$REF))

#Track the allele switches (to change them in the genotype matrix later)
allele_freq_reference$ALLELE_SWITCH <- 0

if(length(idx_switched) > 0){
  allele_freq_reference$ALLELE_SWITCH[idx_switched] <- 1
  
  #Fix the allele switches 
  allele_freq_reference$REF[idx_switched] <- allele_freq_test_red_temp$REF[idx_switched]
  allele_freq_reference$ALT[idx_switched] <- allele_freq_test_red_temp$ALT[idx_switched]
  
  allele_freq_reference$AC[idx_switched] <- allele_freq_reference$AN[idx_switched] - allele_freq_reference$AC[idx_switched]
  allele_freq_reference$AF[idx_switched] <- 1 - allele_freq_reference$AF[idx_switched] 
}

#Now match based on SNP
allele_freq_test$SNP <-  paste(allele_freq_test$CHROM, allele_freq_test$POS, allele_freq_test$REF, allele_freq_test$ALT, sep = ":")
allele_freq_reference$SNP <- paste(chr,  allele_freq_reference$POS, allele_freq_reference$REF, allele_freq_reference$ALT, sep = ":")

allele_freq_reference <- allele_freq_reference[allele_freq_reference$SNP %in% allele_freq_test$SNP,]

###########################################################

#get residual variance and final sample size
n_residual_variance <- read.table(argv$testSizeResidualVariance, header = T)
n <- as.numeric(n_residual_variance[1,1])
residual_variance <- as.numeric(n_residual_variance[1,2])

#read in file of gene limits
genes <- read.table(argv$geneCoord, header = F)
names(genes) <- c("Gene", "Start", "End")

#function to convert chr:pos:ref:alt to pos
id_format <- function(s){
  return(as.numeric(unlist(strsplit(s, split=":"))[2]))
}

#Name of output file 
output_cov_files <- argv$output

idx_gene <- which(genes[[1]] == argv$gene)
gene_start <- as.numeric(genes[idx_gene, 2])
gene_end <- as.numeric(genes[idx_gene, 3])

#For InPSYght data chr is actually chrNum instead of Num (ie. chr5 instead of 5)
chr_long <- paste0("chr", chr)
vars_gt_current_positions <- original_vars

current_range <- ScanVcfParam(which = GRanges(seqnames = Rle(chr_long), ranges = IRanges(start = gene_start, end = gene_end)))

gt_current <- tryCatch(expr = 
                         {gt_current <- readGT(file = argv$vcf, param = current_range)
                         gt_current <- t(gt_current)
                         #Keep only single-allelic variants that are in the test set 
                         gt_current <- gt_current[,vars_gt_current_positions %in% allele_freq_reference[[2]], drop = FALSE] # Note there could be no variants shared between the test and reference
                         
                         gt_current},
                       error = function(cond){
                         
                         samples <- read.table(argv$refsamples, header = F)
                         
                         #ad hoc just try to cut into ten chunks
                         chunk_number <- 10
                         sample_chunks <- split(samples[[1]],
                                                cut(seq_along(samples[[1]]),
                                                    chunk_number,
                                                    labels = FALSE))
                         
                         matrix_list <- list(length = length(sample_chunks))
                         for(j in 1:length(sample_chunks)){
                           print(j)
                           current_range <- ScanVcfParam(samples = as.character(sample_chunks[[j]]),
                                                         which = GRanges(seqnames = Rle(chr_long), ranges = IRanges(start = gene_start, end = gene_end)))
                           
                           matrix_list[[j]] <- readGT(file = argv$vcf, param = current_range)
                           matrix_list[[j]] <- matrix_list[[j]][vars_gt_current_positions %in% allele_freq_reference[[2]],, drop = FALSE]
                           matrix_list[[j]] <- t(matrix_list[[j]])
                         }
                         
                         gt_current <- do.call(base::rbind, matrix_list)
                         return(gt_current)
                         
                       })


#Now switch the entries for the alleles that had switched alleles 
idx_gt_switched <- which(allele_freq_reference$ALLELE_SWITCH == 1)

if(length(idx_gt_switched) > 0){
  gt_current[,idx_gt_switched][gt_current[,idx_gt_switched] %in% c("0/0", "0|0")] <- 2
  gt_current[,idx_gt_switched][gt_current[,idx_gt_switched] %in% c("1/1", "1|1")] <- 0
}

#Now the rest should be the same
if(dim(gt_current)[2] > 0){ #Should be >= 1 variant in the gene
  allele_freq_reference_temp <- allele_freq_reference
  allele_freq_test_temp <- allele_freq_test ###
  
  gt_current[gt_current %in% c("0/0", "0|0")] <- 0
  gt_current[gt_current %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
  gt_current[gt_current %in% c("1/1", "1|1")] <- 2
  
  n_row <- nrow(gt_current)
  n_col <- ncol(gt_current)
  
  gt_current <- suppressWarnings(Matrix(as.numeric(gt_current), nrow = n_row, ncol = n_col, sparse = TRUE)) #Note the NAs for missing values; expect warning here so suppress it
  
  gt_current <- sweep(gt_current, 2, 2*allele_freq_reference_temp$AF) #subtract the allele frequencies rowwise
  
  #Now impute the missing values
  is_na <- which(is.na(gt_current), arr.ind = TRUE)
  gt_current[is_na] <- colMeans(gt_current, na.rm = TRUE)[is_na[,"col"]]
  
  #If there are still NAs left, should make them 0 (these correspond to variants with no calls)
  gt_current[is_na] <- 0
  
  #Now add in the variants that weren't in the reference panel at all 
  ## Need to come up with way to add the missing variants to the matrix
  gt <- matrix(0, nrow = nrow(gt_current), ncol = nrow(allele_freq_test))
  
  #Match columns
  idx_in_reference_panel <- match(allele_freq_reference$SNP, allele_freq_test$SNP)
  
  gt[,idx_in_reference_panel] <- Matrix::as.matrix(gt_current)
  gt <- Matrix(gt, sparse = TRUE)
  
}else{
  stop("No variants in this gene!")
}

#Check that the test allele frequency data frame and genotype matrix have the
#same number of variants - shouldn't ever happen
if(nrow(allele_freq_test) != ncol(gt)){
  stop("Allele frequency data frame and genotype matrix have differing number of variants.")
}

C_mat <- crossprod(gt)
diag(C_mat)[diag(C_mat) == 0] <- 1
D_g <- Diagonal(x = sqrt(diag(C_mat)))
D_mat <- Diagonal(x = sqrt(2*n*allele_freq_test$AF*(1-allele_freq_test$AF)))
D_g_inv <- solve(D_g)
output_mat <- (D_mat %*% D_g_inv %*% C_mat %*% D_g_inv %*% D_mat)/(n*residual_variance)

#Now print the matrix in RAREMETAL form
output_mat_list <- asplit(output_mat, 2)
N <- length(output_mat_list)
output_mat_list_final <- lapply(1:N, function(x) paste0(c(output_mat_list[[x]][x:N], ""), collapse = ","))
output_markers_in_window_list <- lapply(1:N, function(x) paste0(c(allele_freq_test$POS[x:N], ""), collapse = ","))

to_write <- vector(length = 3)
to_write[1] <- "##ProgramName=RareMetalWorker"
to_write[2] <- "##Version=4.15.1"
to_write[3] <- "##CHROM\tCURRENT_POS\tMARKERS_IN_WINDOW\tCOV_MATRICES"

to_write <- c(to_write, paste(allele_freq_test$CHROM, allele_freq_test$POS,
                              output_markers_in_window_list, output_mat_list_final,
                              sep = "\t"))

fileConn <- file(output_cov_files)
writeLines(to_write, fileConn)
close(fileConn)
