library("argparser")
library("stringr")

p <- arg_parser("Program that takes snakemake wildcards and creates file of covariance and score files to be read into raremetal")

p <- add_argument(p, "--chr", help = "Chromosome")
p <- add_argument(p, "--testsize", help = "Size of test set")
p <- add_argument(p, "--refsize", help = "Size of reference panel")
p <- add_argument(p, "--trait", help = "Name of trait")
p <- add_argument(p, "--gene", help = "Name of gene")

p <- add_argument(p, "--summaryFiles", help = "Name of raremetal summary file")
p <- add_argument(p, "--trueCovFiles", help = "Name of raremetal true covariance file")
p <- add_argument(p, "--impCovFiles", help = "Name of raremetal estimated covariance file imputed to 0")

argv <- parse_args(p)

scoreFile <- paste(argv$chr, argv$testsize, argv$gene, argv$trait, "singlevar.score.txt", sep = ".")
scoreFile <- paste(c(paste("rmw", argv$chr, argv$testsize, argv$trait, "", sep = "/"), scoreFile), collapse = "")

trueCov <- paste(argv$chr, argv$testsize, argv$gene, argv$trait, "singlevar.cov.txt", sep = ".")
trueCov <- paste(c(paste("rmw",argv$chr, argv$testsize, argv$trait, "", sep = "/"), trueCov), collapse = "")

impCovUKBB <- paste(argv$chr, argv$testsize, argv$refsize, argv$trait, argv$gene, "imp0.estimated.cov", sep = ".")
impCovUKBB <- paste(c(paste("estCovImp0New", argv$chr, argv$testsize, argv$refsize, argv$trait, "", sep = "/"), impCovUKBB), collapse = "")

#Now create the raremetal readable files 
fileConn <- file(argv$summaryFiles)
writeLines(paste0(scoreFile, ".gz"), fileConn)
close(fileConn)

fileConn <- file(argv$trueCovFiles)
writeLines(paste0(trueCov, ".gz"), fileConn)
close(fileConn)

fileConn <- file(argv$impCovFiles)
writeLines(paste0(impCovUKBB, ".gz"), fileConn)
close(fileConn)
