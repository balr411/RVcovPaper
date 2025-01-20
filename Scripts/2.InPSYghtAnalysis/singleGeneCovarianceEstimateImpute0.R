library("argparser")
library("VariantAnnotation")
library("dplyr")
library("bigmemory")
library("data.table")

p <- arg_parser("Program that computes estimated covariance matrix between variants for a specific gene")
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

#remove variants with 0 AF in test set
allele_freq_test <- allele_freq_test[allele_freq_test[[5]] > 0,]

#Now reduce ref allele freq to only include variants found in test set - note we may want to make sure variants with call rate = 0 in reference are removed
allele_freq_reference <- allele_freq_reference[allele_freq_reference[[2]] %in% allele_freq_test[[2]],]

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

current_range <- ScanVcfParam(which = GRanges(seqnames = Rle(chr), ranges = IRanges(start = gene_start, end = gene_end)))
#Sometimes the dimensions of the genotype matrix are too big to be read into R. So we will use the bigmemory package in R to handle this
variables_list <- tryCatch(expr = {gt_current <- readGT(file = argv$vcf, param = current_range)
gt_current <- t(gt_current)

vars_gt_current <- colnames(gt_current)
vars_gt_current_positions <- as.numeric(sapply(vars_gt_current, id_format))

#Keep only single-allelic variants that are in the test set
gt_current <- gt_current[,vars_gt_current_positions %in% allele_freq_test[[2]]]

if(!is.null(dim(gt_current))){ #This means there are either 0 variants or > 1 variant in the gene
  if(dim(gt_current)[2] > 0){
    allele_freq_reference_temp <- allele_freq_reference[allele_freq_reference[[2]] %in% vars_gt_current_positions,]
    allele_freq_test_temp <- allele_freq_test[allele_freq_test[[2]] %in% vars_gt_current_positions,] #This should be unnecessary
    
    genotype_matrix_current <- matrix(nrow=dim(gt_current)[1], ncol=dim(gt_current)[2])
    
    for(j in 1:dim(gt_current)[2]){
      gt_current[,j][gt_current[,j] %in% c("0/0", "0|0")] <- -2*allele_freq_reference_temp[[7]][j]
      gt_current[,j][gt_current[,j] %in% c("0/1", "1/0", "0|1", "1|0")] <- 1-2*allele_freq_reference_temp[[7]][j]
      gt_current[,j][gt_current[,j] %in% c("1/1", "1|1")] <- 2-2*allele_freq_reference_temp[[7]][j]
      gt_current[,j][gt_current[,j]=="./."] <- mean(as.numeric(gt_current[,j][gt_current[,j]!="./."]))
      genotype_matrix_current[,j] <- as.numeric(gt_current[,j])
    }
  }
}else{#This means there is only 1 variant in the gene
  #Get temporary test and reference allele frequency files with only variants in this gene shared between both sets
  allele_freq_reference_temp <- allele_freq_reference[allele_freq_reference[[2]] %in% vars_gt_current_positions,]
  allele_freq_test_temp <- allele_freq_test[allele_freq_test[[2]] %in% vars_gt_current_positions,]
  
  gt_current[gt_current %in% c("0/0", "0|0")] <- -2*allele_freq_reference_temp[[7]][1]
  gt_current[gt_current %in% c("0/1", "1/0", "0|1", "1|0")] <- 1-2*allele_freq_reference_temp[[7]][1]
  gt_current[gt_current %in% c("1/1", "1|1")] <- 2-2*allele_freq_reference_temp[[7]][1]
  
  gt_current[gt_current == "./."] <- mean(as.numeric(gt_current[gt_current != "./."]))
  
  genotype_matrix_current <- as.numeric(gt_current)
  genotype_matrix_current_square <- genotype_matrix_current^2
  
}

list(vars_gt_current_positions, allele_freq_reference_temp, allele_freq_test_temp, genotype_matrix_current, genotype_matrix_current^2)

},

error = function(cond){
  samples <- read.table(argv$refsamples, header = F)
  num_vars <- sum(original_vars %in% as.numeric(gene_start):as.numeric(gene_end))
  n_samples_vars <- dim(samples)[1]*num_vars
  num_splits <- 0
  while(is.na(n_samples_vars)){
    num_splits <- num_splits + 1
    n_samples_vars <- ceiling((dim(samples)[1]/(2^num_splits)) * num_vars)
  }
  
  chunk_number <- 2^num_splits
  sample_chunks <- split(samples[[1]],
                         cut(seq_along(samples[[1]]),
                             chunk_number,
                             labels = FALSE))
  
  matrix_list <- list(length = length(sample_chunks))
  for(j in 1:length(sample_chunks)){
    current_range <- ScanVcfParam(samples = as.character(sample_chunks[[j]]),
                                  which = GRanges(seqnames = Rle(chr), ranges = IRanges(start = gene_start, end = gene_end)))
    
    matrix_list[[j]] <- t(readGT(file = argv$vcf, param = current_range))
  }
  
  vars_gt_current <- c(colnames(matrix_list[[1]]))
  vars_gt_current_positions <- as.numeric(sapply(vars_gt_current, id_format))
  
  for(j in 1:length(sample_chunks)){
    matrix_list[[j]] <- matrix_list[[j]][,vars_gt_current_positions %in% allele_freq_reference[[2]]]
  }
  
  allele_freq_reference_temp <- allele_freq_reference[allele_freq_reference[[2]] %in% vars_gt_current_positions,]
  allele_freq_test_temp <- allele_freq_test[allele_freq_test[[2]] %in% vars_gt_current_positions,]
  
  genotype_matrix_current <- big.matrix(nrow = dim(samples)[1], ncol =  dim(matrix_list[[1]])[2])
  
  for(j in 1:dim(genotype_matrix_current)[2]){
    current_var <- matrix_list[[1]][,j]
    for(k in 2:length(matrix_list)){
      current_var <- c(current_var, matrix_list[[k]][,j])
    }
    genotype_matrix_current[,j][current_var %in% c("0/0", "0|0")] <- -2*allele_freq_reference_temp[[7]][j]
    genotype_matrix_current[,j][current_var %in% c("0/1", "1/0", "0|1", "1|0")] <- 1-2*allele_freq_reference_temp[[7]][j]
    genotype_matrix_current[,j][current_var %in% c("1/1", "1|1")] <- 2-2*allele_freq_reference_temp[[7]][j]
    genotype_matrix_current[,j][current_var=="./."] <- mean(as.numeric(!is.na(genotype_matrix_current[,j])))
  }
  
  genotype_matrix_current_square <- big.matrix(nrow = dim(genotype_matrix_current)[1], ncol = dim(genotype_matrix_current)[2])
  for(k in 1:dim(genotype_matrix_current_square)[2]){
    genotype_matrix_current_square[,k] <- genotype_matrix_current[,k]^2
  }
  
  return(list(vars_gt_current_positions, allele_freq_reference_temp, allele_freq_test_temp, genotype_matrix_current, genotype_matrix_current_square))
  
})

vars_gt_current_positions <- variables_list[[1]]
allele_freq_reference_temp <- variables_list[[2]]
allele_freq_test_temp <- variables_list[[3]]
genotype_matrix_current <- variables_list[[4]]
genotype_matrix_current_square <- variables_list[[5]]

if(!is.null(dim(genotype_matrix_current))){ #This means there are either 0 variants or > 1 variant in the gene
  if(dim(genotype_matrix_current)[2] > 0){
    
    #get D_W(j) and D_j for each variant j 
    
    if(is.big.matrix(genotype_matrix_current)){
      D_W <- biganalytics::apply(genotype_matrix_current_square, 2, sum)
    }else{
      D_W <- colSums(genotype_matrix_current_square) 
    }
    D <- 2*allele_freq_test_temp[[5]]*(1-allele_freq_test_temp[[5]])*n
    
    N <- dim(allele_freq_test_temp)[1] 
    output_markers_in_window_list <- vector(length=N, mode="list")
    output_covariance_list <- vector(length=N, mode="list")
    
    for (j in 1:N){
      current_var <- allele_freq_test_temp[j, 2]
      output_markers_in_window_list[[j]] <- allele_freq_test_temp[[2]][allele_freq_test_temp[[2]] >= current_var] 
      
      
      #match indices for D_W matrix
      num_in_window <- length(output_markers_in_window_list[[j]])
      idx_markers_in_window <- seq(j, j + num_in_window - 1, 1)
      
      #match indices for D matrix 
      D_start <- j
      idx_D_markers_in_window <- seq(D_start, D_start + num_in_window -1, 1)
      
      mat_i_j <- genotype_matrix_current[,j]*genotype_matrix_current[,idx_markers_in_window] #this works for big matrix but may not if the dimensions are too big so keep an eye on this
      sqrt_cov <- sqrt((D[D_start]*D[idx_D_markers_in_window])/(D_W[j]*D_W[idx_markers_in_window]))
      if(!is.null(dim(mat_i_j))){
        sum_cov <- apply(mat_i_j, 2, sum)
      }else{
        sum_cov <- sum(mat_i_j)
      }
      
      output_covariance_list[[j]] <- (sqrt_cov * sum_cov)/(n*residual_variance)
      
      #Replace all the covariance NA by 0, and the variance NA by the variance estimated from the test set
      if(sum(is.na(output_covariance_list[[j]])) == length(output_covariance_list[[j]])){ #This means current variant is missing
        output_covariance_list[[j]][1] <- (2*allele_freq_test_temp[[5]][j]*(1-allele_freq_test_temp[[5]][j]))/(residual_variance)
        if(length(output_covariance_list[[j]]) > 1){
          output_covariance_list[[j]][2:length(output_covariance_list[[j]])] <- 0
        }
      }else{
        idx_na <- which(is.na(output_covariance_list[[j]]))
        if(length(idx_na) > 0){
          output_covariance_list[[j]][idx_na] <- 0 
        }
      }
    }
    to_write<-vector(length=N+3)
    to_write[1] <- "##ProgramName=RareMetalWorker"
    to_write[2] <- "##Version=4.15.1"
    to_write[3] <- "##CHROM\tCURRENT_POS\tMARKERS_IN_WINDOW\tCOV_MATRICES"
    for (k in 1:N){
      to_write[k+3] <- paste(c(paste(c(paste(c(allele_freq_test_temp[k,1], allele_freq_test_temp[k,2], paste(output_markers_in_window_list[[k]], collapse=",")), collapse="\t"), 
                                       paste(output_covariance_list[[k]], collapse=",")), collapse=",\t"), ","), collapse="")
    }
    fileConn <- file(output_cov_files)
    writeLines(to_write, fileConn)
    close(fileConn)
  }else{ #This means there is only 1 variant in the gene
    #get D_W(j) and D_j for each variant j 
    D_W <- sum(genotype_matrix_current^2)
    D <- 2*allele_freq_test_temp[[5]]*(1-allele_freq_test_temp[[5]])*n
    
    N <- dim(allele_freq_test_temp)[1] 
    output_markers_in_window_list <- vector(length=N, mode="list")
    output_covariance_list <- vector(length=N, mode="list")
    
    for (j in 1:N){
      current_var <- allele_freq_test_temp[j, 2]
      output_markers_in_window_list[[j]] <- allele_freq_test_temp[[2]][allele_freq_test_temp[[2]] >= current_var] 
      
      
      #match indices for D_W matrix
      num_in_window <- length(output_markers_in_window_list[[j]])
      idx_markers_in_window <- seq(j, j + num_in_window - 1, 1)
      
      #match indices for D matrix 
      D_start <- j
      idx_D_markers_in_window <- seq(D_start, D_start + num_in_window -1, 1)
      
      mat_i_j <- genotype_matrix_current[j]*genotype_matrix_current[idx_markers_in_window]
      sqrt_cov <- sqrt((D[D_start]*D[idx_D_markers_in_window])/(D_W[j]*D_W[idx_markers_in_window]))
      if(!is.null(dim(mat_i_j))){
        sum_cov <- apply(mat_i_j, 2, sum)
      }else{
        sum_cov <- sum(mat_i_j)
      }
      
      output_covariance_list[[j]] <- (sqrt_cov * sum_cov)/(n*residual_variance)
      
      if(is.na(output_covariance_list[[j]][1])){
        output_covariance_list[[j]][1] <- (2*allele_freq_test_temp[[5]][j]*(1-allele_freq_test_temp[[5]][j]))/(residual_variance)
      }
      
    }
    
    to_write<-vector(length=N+3)
    to_write[1] <- "##ProgramName=RareMetalWorker"
    to_write[2] <- "##Version=4.15.1"
    to_write[3] <- "##CHROM\tCURRENT_POS\tMARKERS_IN_WINDOW\tCOV_MATRICES"
    for (k in 1:N){
      to_write[k+3] <- paste(c(paste(c(paste(c(allele_freq_test_temp[k,1], allele_freq_test_temp[k,2], paste(output_markers_in_window_list[[k]], collapse=",")), collapse="\t"), 
                                       paste(output_covariance_list[[k]], collapse=",")), collapse=",\t"), ","), collapse="")
    }
    
    fileConn <- file(output_cov_files)
    writeLines(to_write, fileConn)
    close(fileConn)
  }
}else{#this means there is only one variant in the gene
  D_W <- sum(genotype_matrix_current^2)
  D <- 2*allele_freq_test_temp[[5]]*(1-allele_freq_test_temp[[5]])*n
  
  N <- dim(allele_freq_test_temp)[1] 
  output_markers_in_window_list <- vector(length=N, mode="list")
  output_covariance_list <- vector(length=N, mode="list")
  
  for (j in 1:N){
    current_var <- allele_freq_test_temp[j, 2]
    output_markers_in_window_list[[j]] <- allele_freq_test_temp[[2]][allele_freq_test_temp[[2]] >= current_var] 
    
    
    #match indices for D_W matrix
    num_in_window <- length(output_markers_in_window_list[[j]])
    idx_markers_in_window <- seq(j, j + num_in_window - 1, 1)
    
    #match indices for D matrix 
    D_start <- j
    idx_D_markers_in_window <- seq(D_start, D_start + num_in_window -1, 1)
    
    mat_i_j <- genotype_matrix_current*genotype_matrix_current 
    sqrt_cov <- sqrt((D[D_start]*D[idx_D_markers_in_window])/(D_W[j]*D_W[idx_markers_in_window]))
    if(!is.null(dim(mat_i_j))){
      sum_cov <- apply(mat_i_j, 2, sum)
    }else{
      sum_cov <- sum(mat_i_j)
    }
    
    output_covariance_list[[j]] <- (sqrt_cov * sum_cov)/(n*residual_variance)
    
    if(is.na(output_covariance_list[[1]][1])){
      output_covariance_list[[1]][1] <- (2*allele_freq_test_temp[[5]][1]*(1-allele_freq_test_temp[[5]][1]))/(residual_variance)
    }
  }
  
  to_write<-vector(length=N+3)
  to_write[1] <- "##ProgramName=RareMetalWorker"
  to_write[2] <- "##Version=4.15.1"
  to_write[3] <- "##CHROM\tCURRENT_POS\tMARKERS_IN_WINDOW\tCOV_MATRICES"
  for (k in 1:N){
    to_write[k+3] <- paste(c(paste(c(paste(c(allele_freq_test_temp[k,1], allele_freq_test_temp[k,2], paste(output_markers_in_window_list[[k]], collapse=",")), collapse="\t"), 
                                     paste(output_covariance_list[[k]], collapse=",")), collapse=",\t"), ","), collapse="")
  }
  
  fileConn <- file(output_cov_files)
  writeLines(to_write, fileConn)
  close(fileConn)
  
}






