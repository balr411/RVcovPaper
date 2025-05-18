#Script that calls the agg_test() function from our package on chromosome 2

library("argparser")
library("RVcov")

p <- arg_parser("Script that calls the agg_test() function from our package on chromosome 2")

p <- add_argument(p, "--twoStage", help = "Perform two-stage method?")
p <- add_argument(p, "--scoreStat", help = "Score statistic file")
p <- add_argument(p, "--vcf", help = "Name of refrence panel VCF file")
p <- add_argument(p, "--groupFile", help = "Name of group file")
p <- add_argument(p, "--gene", help = "Name of gene")
p <- add_argument(p, "--altCovariancePath", help = "Name of covariance path to use")
p <- add_argument(p, "--altRaremetalPath", help = "Name of RAREMETAL path to use")

argv <- parse_args(p)

score_stat_file <- argv$scoreStat
vcf_file <- argv$vcf
chr <- 2
burden <- TRUE
wburden <- TRUE
SKAT <- TRUE
anno_file <- NULL
anno <- NULL
two_stage <- as.logical(argv$twoStage)
two_stage_threshold <- 3
group_file <- c(broadMissense = argv$groupFile)
pLOF <- pLOF_narrowMissense <- pLOF_broadMissense <- FALSE
altGroupFilePath <- NULL
altCovariancePath <- argv$altCovariancePath
altRaremetalPath <- argv$altRaremetalPath
altRaremetalName <- NULL
mafThreshold <- 0.01
gene <- argv$gene
gene_start <- gene_end <- NULL
hwe <- 0.000001

agg_test(score_stat_file, vcf_file, chr = chr_num, burden, wburden,
         SKAT, anno_file,
         anno, two_stage = two_stage, two_stage_threshold,
         group_file, pLOF,
         pLOF_narrowMissense, pLOF_broadMissense,
         altGroupFilePath, altCovariancePath,
         altRaremetalPath, altRaremetalName,
         mafThreshold, gene,
         gene_start, gene_end, hwe)
