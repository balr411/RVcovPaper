library("argparser")
library("stringr")

p <- arg_parser("Program that creates prefix and calls raremetal for a single gene and multiple different group files")

p <- add_argument(p, "--PTV", help = "Include protein truncating variant group file? Options are T/F, name of file must be 
                  in format {chr}.{test}.{trait}.{gene}.PTV.group.file")
p <- add_argument(p, "--narrowMissense", help = "Include pLOF + missense(narrow) variant group file? Options are T/F, name of file must be 
                  in format {chr}.{test}.{trait}.{gene}.MISSENSE.group.file")
p <- add_argument(p, "--broadMissense", help = "Include pLOF + missense(broad) variant group file? Options are T/F, name of file must be 
                  in format {chr}.{test}.{trait}.{gene}.BROADMISSENSE.group.file")

p <- add_argument(p, "--chr", help = "Chromosome")
p <- add_argument(p, "--testsize", help = "Size of test set")
p <- add_argument(p, "--refsize", help = "Size of reference panel. Only needed when estCov is T.")
p <- add_argument(p, "--traitname", help = "Name of trait")
p <- add_argument(p, "--gene", help = "Name of gene")

p <- add_argument(p, "--trueCov", help = "Run raremetal using true covariance? Options are T/F.")
p <- add_argument(p, "--estCov", help = "Run raremetal using estimated covariance? Options are T/F.
                  By default raremetal is called on full group files where missing covariance in 
                  estimated covariance file is imputed to 0")
p <- add_argument(p, "--all0Cov", help = "Run raremetal using null covariance? Options are T/F.")

p <- add_argument(p, "--covFilesAll0", help = "Location of all0 covariance files for raremetal", default = NULL)
p <- add_argument(p, "--trueCovFiles", help = "Location of true covariance files for raremetal", default = NULL)
p <- add_argument(p, "--estCovFilesImp0", help = "Location of imp0 covariance files for raremetal", default = NULL)
p <- add_argument(p, "--summaryFiles", help = "Location of summary files for raremetal")

p <- add_argument(p, "--done", help = "File name to mark end of pipeline")

argv <- parse_args(p)


#First find which group files to use
group_file_prefix <- c() 

if(as.logical(argv$PTV)){
  group_file_prefix <- c(group_file_prefix, "PTV")
}

if(as.logical(argv$narrowMissense)){
  group_file_prefix <- c(group_file_prefix, "MISSENSE")
}

if(as.logical(argv$broadMissense)){
  group_file_prefix <- c(group_file_prefix, "BROADMISSENSE")
}

#These variables used in all raremetal calls
chr <- argv$chr #note in other scripts we use chr as the numeric chromosome but here it has the form c#
test <- argv$testsize
traitname <- argv$traitname
gene <- argv$gene 

summary_files <- argv$summaryFiles

#For each desired group file, read in first line to see if it actually has any variants
empty_group_file_prefix <- c()
non_empty_group_file_prefix <- c()
for(prefix in group_file_prefix){
  group_file_path <- paste("groupFilesNew", prefix, chr, test, traitname, "", sep = "/")
  group_file <- paste(c(group_file_path, paste(chr, test, traitname, gene, prefix, "group.file", sep = ".")), collapse = "")
  file_current <- readLines(group_file, n = 1)
  l <- unlist(strsplit(file_current, split = "\t"))
  
  if(length(l) == 1){
    empty_group_file_prefix <- c(empty_group_file_prefix, prefix)
  }else{
    non_empty_group_file_prefix <- c(non_empty_group_file_prefix, prefix)
  }
}


#Call raremetal on true covariance files 

if(as.logical(argv$trueCov)){
  cov_file <- argv$trueCovFiles
  for(prefix in group_file_prefix){
    if(prefix %in% non_empty_group_file_prefix){
      group_file_path <- paste("groupFilesNew", prefix, chr, test, traitname, "", sep = "/")
      group_file <- paste(c(group_file_path, paste(chr, test,  traitname, gene, prefix, "group.file", sep = ".")), collapse = "")
      full_prefix <-  paste(chr, test,  traitname, gene, paste0("UNRESTRICTED_", prefix), sep = ".")
      burden_path <- paste("raremetalNewRerun/burden", prefix, chr, test, traitname, gene, "", sep = "/")
      full_prefix <- paste0(burden_path, full_prefix)
      
      cm <- str_glue("/net/dumbo/home/welchr/projects/nullmetal/cmake-build-release/bin/raremetal --summaryFiles {summary_files} --covFiles {cov_file} --groupFile {group_file} --BBeta --SKAT --maf 1 --hwe 0.000001 --prefix {full_prefix}")
      system(cm)
    }else{
      full_prefix <- paste0("UNRESTRICTED_", prefix)
      burden_path <- paste("raremetalNewRerun/burden", prefix, chr, test, traitname, gene, "", sep = "/")
      file_burden <- paste(c(burden_path, paste(argv$chr, argv$testsize, argv$traitname, argv$gene, full_prefix, "meta.BBeta.results", sep = ".")), collapse = "")
      file_SKAT <- paste(c(burden_path, paste(argv$chr, argv$testsize, argv$traitname, argv$gene, full_prefix, "meta.SKAT_.results", sep = ".")), collapse = "")
      file.create(file_burden)
      file.create(file_SKAT)
    }
  }
}

if(as.logical(argv$estCov)){
  ref <- argv$refsize
  for(prefix in group_file_prefix){
    if(prefix %in% non_empty_group_file_prefix){
      cov_file <- argv$estCovFilesImp0
      group_file_path <- paste("groupFilesNew", prefix, chr, test, traitname, "", sep = "/")
      group_file <- paste(c(group_file_path, paste(chr, test, traitname, gene, prefix, "group.file", sep = ".")), collapse = "")
      full_prefix <-  paste(chr, test, ref, traitname, gene, paste0("IMPUTED0_", prefix), sep = ".")
      burden_path <- paste("raremetalNewRerun/burden", prefix, chr, test, ref, traitname, gene, "", sep = "/")
      full_prefix <- paste0(burden_path, full_prefix)
      cm <- str_glue("/net/dumbo/home/welchr/projects/nullmetal/cmake-build-release/bin/raremetal --summaryFiles {summary_files} --covFiles {cov_file} --groupFile {group_file} --BBeta --SKAT --maf 1 --hwe 0.000001 --prefix {full_prefix}")
      system(cm)
    }else{
      full_prefix <- paste0("IMPUTED0_", prefix)
      burden_path <- paste("raremetalNewRerun/burden", prefix, chr, test, ref, traitname, gene, "", sep = "/")
      file_burden <- paste(c(burden_path, paste(argv$chr, argv$testsize, argv$refsize, argv$traitname, argv$gene, full_prefix, "meta.BBeta.results", sep = ".")), collapse = "")
      file_SKAT <- paste(c(burden_path, paste(argv$chr, argv$testsize, argv$refsize, argv$traitname, argv$gene, full_prefix, "meta.SKAT_.results", sep = ".")), collapse = "")
      file.create(file_burden)
      file.create(file_SKAT)
    }
  }
}

#Call raremetal on null covariance 
if(as.logical(argv$all0Cov)){
  cov_file <- argv$covFilesAll0
  for(prefix in group_file_prefix){
    if(prefix %in% non_empty_group_file_prefix){
      group_file_path <- paste("groupFilesNew", prefix, chr, test, traitname, "", sep = "/")
      group_file <- paste(c(group_file_path, paste(chr, test, traitname, gene, prefix, "group.file", sep = ".")), collapse = "")
      full_prefix <-  paste(chr, test, traitname, gene, paste0("ALL0_", prefix), sep = ".")
      burden_path <- paste("raremetalNewRerun/burden", prefix, chr, test, traitname, gene, "", sep = "/")
      full_prefix <- paste0(burden_path, full_prefix)
      cm <- str_glue("/net/dumbo/home/welchr/projects/nullmetal/cmake-build-release/bin/raremetal --summaryFiles {summary_files} --covFiles {cov_file} --groupFile {group_file} --BBeta --SKAT --maf 1 --hwe 0.000001 --prefix {full_prefix}")
      system(cm)
    }else{
      full_prefix <- paste0("ALL0_", prefix)
      burden_path <- paste("raremetalNewRerun/burden", prefix, chr, test, traitname, gene, "", sep = "/")
      file_burden <- paste(c(burden_path, paste(argv$chr, argv$testsize, argv$traitname, argv$gene, full_prefix, "meta.BBeta.results", sep = ".")), collapse = "")
      file_SKAT <- paste(c(burden_path, paste(argv$chr, argv$testsize, argv$traitname, argv$gene, full_prefix, "meta.SKAT_.results", sep = ".")), collapse = "")
      file.create(file_burden)
      file.create(file_SKAT)
    }
  }
}

file.create(argv$done)
