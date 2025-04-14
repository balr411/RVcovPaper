library("argparser")
library("stringr")

p <- arg_parser("Program that takes bgzips and then tabix indexes the imputed0 raremetal output files")

p <- add_argument(p, "--estCovImp0UKBB", help = "Name of unzipped UKBB external reference panel estimated covariance file using all variants")

argv <- parse_args(p)

#First bgzip everything
est_cov_imp0_UKBB <- argv$estCovImp0UKBB
cm <- str_glue("bgzip {est_cov_imp0_UKBB}")
system(cm)

#Now tabix index everything
imp0_cov_UKBB_gzip <- paste0(est_cov_imp0_UKBB, ".gz")
cm <- str_glue("tabix -c \"#\" -s 1 -b 2 -e 2 {imp0_cov_UKBB_gzip}")
system(cm)

