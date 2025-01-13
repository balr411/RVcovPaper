library("argparser")
library("stringr")

p <- arg_parser("Program that takes bgzips and then tabix indexes the raremetal output files")

p <- add_argument(p, "--chr", help = "Chromosome")
p <- add_argument(p, "--testsize", help = "Size of test set")
p <- add_argument(p, "--refsize", help = "Size of reference panel")
p <- add_argument(p, "--trait", help = "Name of trait")
p <- add_argument(p, "--gene", help = "Name of gene")
p <- add_argument(p, "--covAll0", help = "Name of unzipped all0 covariance file")
p <- add_argument(p, "--estCovImp0", help = "Name of unzipped external reference panel estimated covariance file using all variants")


argv <- parse_args(p)


#First bgzip everything
score_file <- paste0("rmw/", argv$chr, "/", argv$testsize, "/", argv$trait, "/", argv$chr, ".", argv$testsize, ".", argv$gene, ".", argv$trait, ".singlevar.score.txt")
cm <- str_glue("bgzip {score_file}")
if(file.exists(score_file)){
  system(cm)
}
  
true_cov <- paste0("rmw/", argv$chr, "/", argv$testsize, "/", argv$trait, "/", argv$chr, ".", argv$testsize, ".", argv$gene, ".", argv$trait, ".singlevar.cov.txt")
cm <- str_glue("bgzip {true_cov}")
if(file.exists(true_cov)){
  system(cm)
}

all0_cov <- argv$covAll0
cm <- str_glue("bgzip {all0_cov}")
system(cm)

est_cov_imp0 <- argv$estCovImp0
cm <- str_glue("bgzip {est_cov_imp0}")
system(cm)

#Now tabix everything 
score_file_gzip <- paste0(score_file, ".gz")
cm <- str_glue("tabix -c \"#\" -s 1 -b 2 -e 2 {score_file_gzip}")
if(file.exists(score_file_gzip)){
  system(cm)
}

true_cov_gzip <- paste0(true_cov, ".gz")
cm <- str_glue("tabix -c \"#\" -s 1 -b 2 -e 2 {true_cov_gzip}")
if(file.exists(true_cov_gzip)){
  system(cm)
}
  
all0_cov_gzip <- paste0(all0_cov, ".gz")
cm <- str_glue("tabix -c \"#\" -s 1 -b 2 -e 2 {all0_cov_gzip}")
system(cm)

est_cov_imp0_gzip <- paste0(est_cov_imp0, ".gz")
cm <- str_glue("tabix -c \"#\" -s 1 -b 2 -e 2 {est_cov_imp0_gzip}")
system(cm)
