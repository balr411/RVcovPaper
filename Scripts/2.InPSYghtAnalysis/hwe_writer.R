#Script that reads in raremetalworker allele frequency files and writes out
#data frame of those not passing hwe p-value threshold

library("stringr")
library("data.table")
library("dplyr")

for(i in 1:22){
  print(i)
  chr_curr <- paste0("c", i)
  file_curr <- str_glue("InPSYght_hwe/{chr_curr}.randomCov.singlevar.score.txt.gz")
  df_inpsyght_curr <- fread(cmd = str_glue("zgrep -v ^## {file_curr}"), header = TRUE, data.table = FALSE)
  df_inpsyght_curr <- filter(df_inpsyght_curr, HWE_PVALUE < 0.000001)
  
  file_to_write <- str_glue("InPSYght_hwe/{chr_curr}.InPSYght.failHWE.vars")
  fwrite(df_inpsyght_curr, file = file_to_write, quote = F, row.names = F, col.names = F, sep="\t")
}
