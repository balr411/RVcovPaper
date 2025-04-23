library("argparser")
library("VariantAnnotation")
library("dplyr")
library("bigmemory")
library("stringr")

p <- arg_parser("Program that computes an estimated covariance matrix for a single gene by only outputting variance, with all covariance set to 0")

p <- add_argument(p, "--geneCoord", help = "File containing genomic coordinates of genes")
p <- add_argument(p, "--gene", help = "Name of gene")
p <- add_argument(p, "--chr", help = "Chromosome to be worked on")
p <- add_argument(p, "--testsize", help = "Size of test set used")
p <- add_argument(p, "--trait", help = "Name of trait to use")

argv <- parse_args(p)

genes <- read.table(argv$geneCoord, header = F)
chr_long <- argv$chr
chr <- as.numeric(unlist(strsplit(argv$chr, split = "c"))[2])
test <- as.numeric(argv$testsize)
trait <- argv$trait
gene <- argv$gene


for(l in 1:20){
  study <- paste0("S", l)
  argv$testallelefreq <- str_glue("testAF/{chr_long}/{study}/{test}/{trait}/{chr_long}.{study}.{test}.{trait}.{gene}.test.allele.freq")
  argv$testSizeResidualVariance <- str_glue("residVar/{chr_long}/{study}/{test}/{trait}/{chr_long}.{study}.{test}.{gene}.{trait}.testSize_residualVariance")
  argv$output <- str_glue("estCovAll0IndStudy/{chr_long}/{study}/{test}/{trait}/{chr_long}.{study}.{test}.{trait}.{gene}.all0.estimated.cov.gz")
    
  allele_freq_test <- read.table(argv$testallelefreq, header = F)
  names(allele_freq_test) <- c("CHROM", "POS", "REF", "ALT", "AF", "AC", "PVAL")
  n_residual_variance <- read.table(argv$testSizeResidualVariance, header = T)
  n <- as.numeric(n_residual_variance[1,1])
  residual_variance <- as.numeric(n_residual_variance[1,2])
  
  #Remove multiallelic variants in test set 
  allele_freq_test <- anti_join(allele_freq_test, allele_freq_test[duplicated(allele_freq_test[,2]),], by = "POS")

  allele_freq_test <- allele_freq_test[allele_freq_test[[5]] > 0,]
  
  output_cov_files_split <- unlist(strsplit(argv$output, split = ".", fixed = TRUE))
  output_cov_files <- paste(output_cov_files_split[-length(output_cov_files_split)], collapse = ".")
  
  idx_gene <- which(genes[[1]] == argv$gene)
  gene_start <- as.numeric(genes[idx_gene, 2])
  gene_end <- as.numeric(genes[idx_gene, 3])
  
  vars_current_idx <- which(allele_freq_test[[2]] %in% gene_start:gene_end)
  
  
  if(length(vars_current_idx) > 0){
    allele_freq_test_temp <- allele_freq_test[vars_current_idx,]
    
    N <- dim(allele_freq_test_temp)[1] 
    output_markers_in_window_list <- vector(length=N, mode="list")
    output_covariance_list <- vector(length=N, mode="list")
    
    for(j in 1:N){
      current_var <- allele_freq_test_temp[j, 2]
      output_markers_in_window_list[[j]] <- allele_freq_test_temp[[2]][allele_freq_test_temp[[2]] >= current_var] 
      
      output_covariance_list[[j]] <- rep(0, length(output_markers_in_window_list[[j]]))
      output_covariance_list[[j]][1] <- (2*allele_freq_test_temp[[5]][j]*(1-allele_freq_test_temp[[5]][j]))/(residual_variance)
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
    
    #Now bgzip and tabix the covariance
    all0Cov <- output_cov_files
    cm <- str_glue("bgzip {all0Cov}")
    system(cm)
    
    all0Cov_gzip <- paste0(output_cov_files, ".gz")
    cm <- str_glue("tabix -c \"#\" -s 1 -b 2 -e 2 {all0Cov_gzip}")
    system(cm)
    
  }else{ #Means there are no variants in the gene - output an empty file so Snakemake doesn't fail, also output a .ignore file 
    file.create(paste0(output_cov_files, ".gz"))
    file.create(paste0(output_cov_files, ".gz.tbi"))
    file_ignore <- paste0(output_cov_files, ".ignore")
    file.create(file_ignore)
  }
  
  
  
  
  
}




