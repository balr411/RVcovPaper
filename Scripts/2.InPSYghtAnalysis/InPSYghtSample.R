library("argparser")

p <- arg_parser("R script that takes as argument the desired reference panel size and file of samples outputs a file of samples of the reference panel size")

p <- add_argument(p, "--samples", help = "Name of sample file")
p <- add_argument(p, "--refsize", help = "Size of reference panel")
p <- add_argument(p, "--output", help = "Name of file to write to")

argv <- parse_args(p)

samps <- read.table(argv$samples, header = F)[[1]]
nref <- as.numeric(argv$refsize)

set.seed(nref)
samps_ref <- sample(samps, nref, replace = FALSE)

fileConn <- file(argv$output)
writeLines(samps_ref, fileConn)
close(fileConn)
