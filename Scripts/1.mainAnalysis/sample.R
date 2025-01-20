library("argparser")

p <- arg_parser("Script that reads in UKBB test set size, reference panel size, a 
                list of sample IDs, an exclusion list, and a list of related samples
                and then samples and outputs a set of samples to serve as test and 
                reference panel")

p <- add_argument(p, "--testsize", help = "Size of test set")
p <- add_argument(p, "--refsize", help = "Size of reference set")
p <- add_argument(p, "--phenosamples", help = "File containing phenotype sample IDs", default = "Phenotypes/sample.id")
p <- add_argument(p, "--exclusionlist", help = "File containing samples to be excluded from analysis due to withdrawn consent")
p <- add_argument(p, "--relatedsamples", help = "File containing samples to be excluded based on being related to someone")
p <- add_argument(p, "--testsamples", help = "Output file containing test sample IDs")
p <- add_argument(p, "--refsamples", help = "Output file containing reference sample IDs")

argv <- parse_args(p)

phenotype_samples <- read.table(as.character(argv$phenosamples), header=T) #note sometimes header should be FALSE, so check before running
exclusion_list <- read.csv(as.character(argv$exclusionlist), header=F)
related_samples <- read.table(argv$relatedsamples, header = F)

#Note below we would normally intersect with the chromosome samples but I am now doing this before writing the phenotype samples file
common_samples <- phenotype_samples[[1]]

#Check to see if any of these samples overlap with exclusion list
common_int_exclusion <- intersect(common_samples, exclusion_list[[1]])
if(length(common_int_exclusion)>0){
  common_samples <- common_samples[-which(common_samples %in% common_int_exclusion)]
}

#Now check to see if any of these samples are related
common_int_rel<- intersect(common_samples, related_samples[[1]])
if(length(common_int_rel) > 0){
  common_samples <- common_samples[-which(common_samples %in% common_int_rel)]
}

#Note that I am also excluding all related individuals and those on the exclusion list already so this shouldn't delete anyone

n_test <- as.numeric(argv$testsize)
n_ref <- as.numeric(argv$refsize)

common_samples <- as.character(common_samples)

set.seed(n_test)
test_set_idx <- sample.int(length(common_samples), n_test)
test_set <- common_samples[test_set_idx]
common_samples_new <- common_samples[-test_set_idx]
set.seed(n_ref)
reference_set_idx <- sample.int(length(common_samples_new), n_ref)
reference_set <- common_samples_new[reference_set_idx]

to_write_test <- vector(length=n_test)
to_write_ref <- vector(length=n_ref)

for(i in 1:n_test){
  to_write_test[i] <- as.character(test_set[i])
}

for(i in 1:n_ref){
  to_write_ref[i] <- as.character(reference_set[i])
}

fileConn <- file(as.character(argv$testsamples))
writeLines(to_write_test, fileConn)
close(fileConn)

fileConn <- file(as.character(argv$refsamples))
writeLines(to_write_ref, fileConn)
close(fileConn)



