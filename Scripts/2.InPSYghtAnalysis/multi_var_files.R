#Script that reads in InPSYght allele frequency files and outputs a list of the multi-allelic
#positions

library("stringr")
library("data.table")
library("dplyr")

for(i in 1:22){
  print(i)
  chr <- paste0("c", i)
  inpsyght_file <- str_glue("InPSYght_multi_allelic_investigation/InPSYght_full_allele_freq_noTopMED/{chr}.noTopMED.allele.freq.gz")
  df_inpsyght <- fread(cmd = paste0("zcat ", inpsyght_file), data.table = F)
  names(df_inpsyght) <- c("CHR", "POS", "REF", "ALT", "AN", "AC")
  
  df_dup <- df_inpsyght[df_inpsyght$POS %in% df_inpsyght$POS[duplicated(df_inpsyght$POS)],]
  
  file_to_write <- str_glue("InPSYght_duplicated/{chr}.InPSYght.duplicated.vars")
  fwrite(df_dup, file = file_to_write, quote = F, row.names = F, col.names = F, sep="\t")
}