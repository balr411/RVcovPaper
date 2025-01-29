library("argparser")
library("stringr")

p <- arg_parser("Program that takes bgzips and then tabix indexes the imputed0 raremetal output files")

p <- add_argument(p, "--estCovAll0", help = "Name of unzipped all0 covariance estimate")
p <- add_argument(p, "--trueCov", help = "Name of unzipped raremetalworker estimated covariance")
p <- add_argument(p, "--scoreFile", help = "Name of raremetalworker score statistic file")

argv <- parse_args(p)

#First bgzip everything
all0Cov <- argv$estCovAll0
cm <- str_glue("bgzip {all0Cov}")
system(cm)

trueCov <- argv$trueCov
cm <- str_glue("bgzip {trueCov}")
system(cm)

scoreFile <- argv$scoreFile
cm <- str_glue("bgzip {scoreFile}")
system(cm)

#Now tabix index everything
all0Cov_gzip <- paste0(all0Cov, ".gz")
cm <- str_glue("tabix -c \"#\" -s 1 -b 2 -e 2 {all0Cov_gzip}")
system(cm)

trueCov_gzip <- paste0(trueCov, ".gz")
cm <- str_glue("tabix -c \"#\" -s 1 -b 2 -e 2 {trueCov_gzip}")
system(cm)

scoreFile_gzip <- paste0(scoreFile, ".gz")
cm <- str_glue("tabix -c \"#\" -s 1 -b 2 -e 2 {scoreFile_gzip}")
system(cm)
