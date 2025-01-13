library("argparser")
library("stringr")

p <- arg_parser("Program that creates prefix and calls raremetal for a single gene and multiple different group files")

p <- add_argument(p, "--summaryfiles", help = "Name of summary files")

p <- add_argument(p, "--chr", help = "Chromosome")
p <- add_argument(p, "--testsize", help = "Size of test set")
p <- add_argument(p, "--refsize", help = "Size of reference panel")
p <- add_argument(p, "--traitname", help = "Name of trait")
p <- add_argument(p, "--gene", help = "Name of gene")

p <- add_argument(p, "--trueCov", help = "Run raremetal using true covariance? Options are T/F. If shared = T,
                  raremetal is also run on the group file using only shared variants")
p <- add_argument(p, "--estCov", help = "Run raremetal using estimated covariance? Options are T/F. If shared = T,
                  raremetal is also run using covariance file with no missing variants. By default raremetal is called
                  on full group files where missing covariance in estimated covariance file is imputed to 0")
p <- add_argument(p, "--all0Cov", help = "Run raremetal using null covariance? Options are T/F. Shared = T has no 
                  effect here")

p <- add_argument(p, "--covFilesAll0", help = "Location of all0 covariance files for raremetal", default = NULL)
p <- add_argument(p, "--trueCovFiles", help = "Location of true covariance files for raremetal")
p <- add_argument(p, "--estCovFilesImp0", help = "Location of imp0 covariance files for raremetal")
p <- add_argument(p, "--summaryFiles", help = "Location of summary files for raremetal")

p <- add_argument(p, "--done", help = "File name to mark end of pipeline")

argv <- parse_args(p)

group_file_prefix <- c("BROADMISSENSE") 

#Note that broad missense group file must be in the format "groupFilesNew/BROADMISSENSE/{chr}/{test}/{ref}/{trait}/{chr}.{test}.{ref}.{trait}.{gene}.BROADMISSENSE.group.file"

#These variables used in all raremetal calls
chr <- argv$chr #note in other scripts we use chr as the numeric chromosome but here it has the form c#
test <- argv$testsize
ref <- argv$refsize
traitname <- argv$traitname
gene <- argv$gene 

summary_files <- argv$summaryFiles

#For each desired group file, read in first line to see if it actually has any variants
#No need for the for loop here but this was copied from generatePLOF_Missense_Nonsyn.R so it is staying
empty_group_file_prefix <- c()
non_empty_group_file_prefix <- c()
for(prefix in group_file_prefix){
  group_file_path <- paste("groupFilesNew", prefix, chr, test, ref, traitname, "", sep = "/")
  group_file <- paste(c(group_file_path, paste(chr, test, ref, traitname, gene, prefix, "group.file", sep = ".")), collapse = "")
  file_current <- readLines(group_file, n = 1)
  l <- unlist(strsplit(file_current, split = "\t"))
  
  if(length(l) == 1){
    empty_group_file_prefix <- c(empty_group_file_prefix, prefix)
  }else{
    non_empty_group_file_prefix <- c(non_empty_group_file_prefix, prefix)
  }
}

if(as.logical(argv$trueCov)){
  #First check to see that there are any variants in gene and if not create empty file  (if no variants, all0 cov file is empty)
  file_true_prefix <- paste("estCovAll0New", argv$chr, argv$testsize, argv$refsize, argv$traitname, "", sep = "/")
  file_true <- paste(c(file_true_prefix, paste(argv$chr, argv$testsize, argv$refsize, argv$traitname, argv$gene, ".all0.estimated.cov.ignore", sep = ".")), collapse = "")
  if(file.exists(file_true)){
    for(prefix in group_file_prefix){
      full_prefix <- paste0("UNRESTRICTED_", prefix)
      burden_path <- paste("raremetalNew/burden", prefix, chr, test, ref, traitname, gene, "", sep = "/")
      file_burden <- paste(c(burden_path, paste(argv$chr, argv$testsize, argv$refsize, argv$traitname, argv$gene, full_prefix, "meta.burden.results", sep = ".")), collapse = "")
      file.create(file_burden)
    }
  }else{
    cov_file <- argv$trueCovFiles
    for(prefix in group_file_prefix){
      if(prefix %in% non_empty_group_file_prefix){
        group_file_path <- paste("groupFilesNew", prefix, chr, test, ref, traitname, "", sep = "/")
        group_file <- paste(c(group_file_path, paste(chr, test, ref, traitname, gene, prefix, "group.file", sep = ".")), collapse = "")
        full_prefix <-  paste(chr, test, ref, traitname, gene, paste0("UNRESTRICTED_", prefix), sep = ".")
        burden_path <- paste("raremetalNew/burden", prefix, chr, test, ref, traitname, gene, "", sep = "/")
        full_prefix <- paste0(burden_path, full_prefix)
        
        cm <- str_glue("/net/snowwhite/home/welchr/projects/raremetal/build/release/bin/raremetal --summaryFiles {summary_files} --covFiles {cov_file} --groupFile {group_file} --burden --maf 1 --hwe 0.000001 --prefix {full_prefix}")
        system(cm)
      }else{
        full_prefix <- paste0("UNRESTRICTED_", prefix)
        burden_path <- paste("raremetalNew/burden", prefix, chr, test, ref, traitname, gene, "", sep = "/")
        file_burden <- paste(c(burden_path, paste(argv$chr, argv$testsize, argv$refsize, argv$traitname, argv$gene, full_prefix, "meta.burden.results", sep = ".")), collapse = "")
        file.create(file_burden)
      }
    }
  }
}

#Call raremetal on estimated covFiles

if(as.logical(argv$estCov)){
  #First check to see that there are any variants in gene and if not create empty file  (if no variants, all0 cov file is empty)
  file_true_prefix <- paste("estCovAll0New", argv$chr, argv$testsize, argv$refsize, argv$traitname, "", sep = "/")
  file_true <- paste(c(file_true_prefix, paste(argv$chr, argv$testsize, argv$refsize, argv$traitname, argv$gene, ".all0.estimated.cov.ignore", sep = ".")), collapse = "")
  if(file.exists(file_true)){
    for(prefix in group_file_prefix){
      full_prefix <- paste0("IMPUTED0_", prefix)
      burden_path <- paste("raremetalNew/burden", prefix, chr, test, ref, traitname, gene, "", sep = "/")
      file_burden <- paste(c(burden_path, paste(argv$chr, argv$testsize, argv$refsize, argv$traitname, argv$gene, full_prefix, "meta.burden.results", sep = ".")), collapse = "")
      file.create(file_burden)
    }
  }else{
    cov_file <- argv$estCovFilesImp0
    for(prefix in group_file_prefix){
      if(prefix %in% non_empty_group_file_prefix){
        group_file_path <- paste("groupFilesNew", prefix, chr, test, ref, traitname, "", sep = "/")
        group_file <- paste(c(group_file_path, paste(chr, test, ref, traitname, gene, prefix, "group.file", sep = ".")), collapse = "")
        full_prefix <-  paste(chr, test, ref, traitname, gene, paste0("IMPUTED0_", prefix), sep = ".")
        burden_path <- paste("raremetalNew/burden", prefix, chr, test, ref, traitname, gene, "", sep = "/")
        full_prefix <- paste0(burden_path, full_prefix)
        cm <- str_glue("/net/snowwhite/home/welchr/projects/raremetal/build/release/bin/raremetal --summaryFiles {summary_files} --covFiles {cov_file} --groupFile {group_file} --burden --maf 1 --hwe 0.000001 --prefix {full_prefix}")
        system(cm)
      }else{
        full_prefix <- paste0("IMPUTED0_", prefix)
        burden_path <- paste("raremetalNew/burden", prefix, chr, test, ref, traitname, gene, "", sep = "/")
        file_burden <- paste(c(burden_path, paste(argv$chr, argv$testsize, argv$refsize, argv$traitname, argv$gene, full_prefix, "meta.burden.results", sep = ".")), collapse = "")
        file.create(file_burden)
      }
    }
  }
}

#Call raremetal on null covariance 
if(as.logical(argv$all0Cov)){
  #First check to see that there are any variants in gene and if not create empty file  (if no variants, all0 cov file is empty)
  file_true_prefix <- paste("estCovAll0New", argv$chr, argv$testsize, argv$refsize, argv$traitname, "", sep = "/")
  file_true <- paste(c(file_true_prefix, paste(argv$chr, argv$testsize, argv$refsize, argv$traitname, argv$gene, ".all0.estimated.cov.ignore", sep = ".")), collapse = "")
  if(file.exists(file_true)){
    for(prefix in group_file_prefix){
      full_prefix <- paste0("ALL0_", prefix)
      burden_path <- paste("raremetalNew/burden", prefix, chr, test, ref, traitname, gene, "", sep = "/")
      file_burden <- paste(c(burden_path, paste(argv$chr, argv$testsize, argv$refsize, argv$traitname, argv$gene, full_prefix, "meta.burden.results", sep = ".")), collapse = "")
      file.create(file_burden)
    }
  }else{
    cov_file <- argv$covFilesAll0
    for(prefix in group_file_prefix){
      if(prefix %in% non_empty_group_file_prefix){
        group_file_path <- paste("groupFilesNew", prefix, chr, test, ref, traitname, "", sep = "/")
        group_file <- paste(c(group_file_path, paste(chr, test, ref, traitname, gene, prefix, "group.file", sep = ".")), collapse = "")
        full_prefix <-  paste(chr, test, ref, traitname, gene, paste0("ALL0_", prefix), sep = ".")
        burden_path <- paste("raremetalNew/burden", prefix, chr, test, ref, traitname, gene, "", sep = "/")
        full_prefix <- paste0(burden_path, full_prefix)
        cm <- str_glue("/net/snowwhite/home/welchr/projects/raremetal/build/release/bin/raremetal --summaryFiles {summary_files} --covFiles {cov_file} --groupFile {group_file} --burden --maf 1 --hwe 0.000001 --prefix {full_prefix}")
        system(cm)
      }else{
        full_prefix <- paste0("ALL0_", prefix)
        burden_path <- paste("raremetalNew/burden", prefix, chr, test, ref, traitname, gene, "", sep = "/")
        file_burden <- paste(c(burden_path, paste(argv$chr, argv$testsize, argv$refsize, argv$traitname, argv$gene, full_prefix, "meta.burden.results", sep = ".")), collapse = "")
        file.create(file_burden)
      }
    }
  }
}

file.create(argv$done)












