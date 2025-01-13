library("argparser")
library("stringr")

p <- arg_parser("Program that takes snakemake wildcards and creates file of covariance and score files to be read into raremetal")

p <- add_argument(p, "--chr", help = "Chromosome")
p <- add_argument(p, "--testsize", help = "Size of test set")
p <- add_argument(p, "--refsize", help = "Size of reference panel")
p <- add_argument(p, "--trait", help = "Name of trait")
p <- add_argument(p, "--gene", help = "Name of gene")
p <- add_argument(p, "--refallelefreq", help = "Reference panel allele frequency file")
p <- add_argument(p, "--testallelefreq", help = "Test set allele frequency file")

p <- add_argument(p, "--summaryFiles", help = "Name of raremetal summary file")
p <- add_argument(p, "--trueCovFiles", help = "Name of raremetal true covariance file")
p <- add_argument(p, "--impCovFiles", help = "Name of raremetal estimated covariance file imputed to 0")
p <- add_argument(p, "--all0CovFiles", help = "Name of raremetal estimated covariance file with all covariances set to 0")

argv <- parse_args(p)

scoreFile <- paste(argv$chr, argv$testsize, argv$gene, argv$trait, "singlevar.score.txt", sep = ".")
scoreFile <- paste(c(paste("rmw", argv$chr, argv$testsize, argv$trait, "", sep = "/"), scoreFile), collapse = "")

trueCov <- paste(argv$chr, argv$testsize, argv$gene, argv$trait, "singlevar.cov.txt", sep = ".")
trueCov <- paste(c(paste("rmw",argv$chr, argv$testsize, argv$trait, "", sep = "/"), trueCov), collapse = "")

impCov <- paste(argv$chr, argv$testsize, argv$refsize, argv$trait, argv$gene, "imp0.estimated.cov", sep = ".")
impCov <- paste(c(paste("estCovImp0New", argv$chr, argv$testsize, argv$refsize, argv$trait, "", sep = "/"), impCov), collapse = "")

all0Cov <- paste(argv$chr, argv$testsize, argv$refsize, argv$trait, argv$gene, "all0.estimated.cov", sep = ".")
all0Cov <- paste(c(paste("estCovAll0New", argv$chr, argv$testsize, argv$refsize, argv$trait, "", sep = "/"), all0Cov), collapse = "")

#Now check to see if all0Cov is empty
lines_all0_cov <- readLines(all0Cov, n = 1)

#Note that if all0Cov is empty then so is impCov since they should be the same size 
#If they are empty, create the same files as for the estimated cov

if(length(lines_all0_cov) == 0){
  file_est <- all0Cov
  name_file <- paste0(file_est, ".ignore")
  file.create(name_file)
  
  file_est <- impCov
  name_file <- paste0(file_est, ".ignore")
  file.create(name_file)
}

#Now create the raremetal readable files 
fileConn <- file(argv$summaryFiles)
writeLines(paste0(scoreFile, ".gz"), fileConn)
close(fileConn)

fileConn <- file(argv$trueCovFiles)
writeLines(paste0(trueCov, ".gz"), fileConn)
close(fileConn)

fileConn <- file(argv$impCovFiles)
writeLines(paste0(impCov, ".gz"), fileConn)
close(fileConn)

fileConn <- file(argv$all0CovFiles)
writeLines(paste0(all0Cov, ".gz"), fileConn)
close(fileConn)



