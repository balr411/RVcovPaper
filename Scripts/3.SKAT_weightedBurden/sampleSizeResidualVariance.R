library("argparser")

p <- arg_parser("Program that reads in score statistic file and outputs file with sample size and residual variance")
p <- add_argument(p, "--scoreStat", help = "File containing score statisic information")
p <- add_argument(p, "--output", help = "File to output to")

argv <- parse_args(p)

ss <- readLines(argv$scoreStat)

sampleSize <- strsplit(ss[4], split = "=")[[1]][2]

stop <- FALSE
i <- 1
while(!stop){
  if(strsplit(ss[i], split = "\t")[[1]][1] == "##Sigma_e2_Hat"){
    residual_variance <- strsplit(ss[i], split = "\t")[[1]][2]
    stop <- TRUE
  }
  i <- i + 1
}

to_return <- data.frame(sampleSize = sampleSize, residualVariance = residual_variance)

write.table(to_return, file = argv$output)