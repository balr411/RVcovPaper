library("argparser")
library("VariantAnnotation")
library("dplyr")
library("bigmemory")
library("stringr")

p <- arg_parser("Program that computes estimated covariance matrix between variants for a specific gene using an external reference panel")
p <- add_argument(p, "--refallelefreq", help = "Reference panel allele frequency file for gene in form CHROM,POS,REF,ALT,AN,AC")
p <- add_argument(p, "--vcf", help = "Reference set vcf to read in")
p <- add_argument(p, "--geneCoord", help = "File containing genomic coordinates of genes")
p <- add_argument(p, "--gene", help = "Name of gene")
p <- add_argument(p, "--chr", help = "Chromosome to be worked on")
p <- add_argument(p, "--trait", help = "Name of trait to use")
p <- add_argument(p, "--refsamples", help = "File of reference samples")
p <- add_argument(p, "--refsize", help = "Size of reference panel")
p <- add_argument(p, "--testsize", help = "Size of test set")

argv <- parse_args(p)
chr_long <- argv$chr

test <- as.numeric(argv$testsize)
ref <- as.numeric(argv$refsize)
trait <- argv$trait
gene <- argv$gene

for(l in 1:20){
  study <- paste0("S", l)
  testallelefreq <- str_glue("testAF/{chr_long}/{study}/{test}/{trait}/{chr_long}.{study}.{test}.{trait}.{gene}.test.allele.freq")
  testSizeResidualVariance <- str_glue("residVar/{chr_long}/{study}/{test}/{trait}/{chr_long}.{study}.{test}.{gene}.{trait}.testSize_residualVariance")
  output <- str_glue("estCovImp0NewIndStudy/{chr_long}/{study}/{test}/{ref}/{trait}/{chr_long}.{study}.{test}.{ref}.{trait}.{gene}.imp0.estimated.cov.gz")
  
  pooled_af <- "nothing"
  
  refAlleleFreq <- argv$refallelefreq
  vcf <- argv$vcf
  ref_samples <- argv$refsamples
  genes <- argv$geneCoord
  
  cm <- str_glue("Rscript Scripts/4.metaAnalysis/singleGeneCovarianceEstimateImpute0Zipped.R --testallelefreq {testallelefreq} --refallelefreq {refAlleleFreq} --vcf {vcf} --refsamples {ref_samples} --testSizeResidualVariance {testSizeResidualVariance} --geneCoord {genes} --gene {gene} --chr {chr_long} --refsize {ref} --testsize {test} --output {output} --pooledAF {pooled_af}")
  
  system(cm)
  
  
}

