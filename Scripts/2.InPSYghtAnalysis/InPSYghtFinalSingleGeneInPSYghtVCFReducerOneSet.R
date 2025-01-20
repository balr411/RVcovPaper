library("argparser")
library("stringr")

p <- arg_parser("Program that takes samples, gene, and chromosome, reduces bcf to include only the samples
                and only variants included in the gene boundaries and tabix indexes it")

p <- add_argument(p, "--samples", help = "Name of file containing desired samples, one per line")
p <- add_argument(p, "--gene", help = "Name of desired gene")
p <- add_argument(p, "--chr", help = "Chromosome that gene is on")
p <- add_argument(p, "--output", help = "Name of desired output vcf file")

argv <- parse_args(p)

gene_table <- read.table(paste0("../Genes/", argv$chr, ".genes"), header = FALSE)

idx <- which(gene_table[[1]] == argv$gene)
range_1 <- gene_table[[2]][idx]
range_2 <- gene_table[[3]][idx]

chr <- paste0("chr", as.numeric(unlist(strsplit(argv$chr, split = "c"))[2]))
bcf <- str_glue("../InPSYght/bcfFiles/freeze9.inpsyght-topmed-aa.{chr}.filtered.gtonly.minDP0.minAC1.bcf")
chr_num <- as.numeric(unlist(strsplit(argv$chr, split = "c"))[2])
range <- paste0(chr, ":", as.numeric(range_1), "-", as.numeric(range_2)) #Note the as.numeric() could be important o/w space gets added

#Test vcf
test_samples <- argv$samples
test_output <- argv$output


cm <- str_glue("bcftools view -Oz --samples-file {test_samples} -r {range} {bcf} > {test_output}")
system(cm)

cm <- str_glue("tabix -p vcf {test_output}")
system(cm)