#This script reads in the results of the PLINK HWE experiment and the file containing
#the variant information for InPSYght from RAREMETALWORKER and writes it to files
#InPSYght_duplicated/{chr}.InPSYght.duplicated.vars and InPSYght_hwe/{chr}.InPSYght.failHWE.vars

#Note that the files {chr_curr}.positions.gz were obtained by running bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' freeze9.inpsyght-topmed-aa.{chr_curr}.filtered.gtonly.minDP0.minAC1.bcf | bgzip > {chr_curr}.positions.gz
#And the files {chr_curr}_plink_hwe.hwe.gz were obtained by running plink19 --bcf freeze9.inpsyght-topmed-aa.{chr_curr}.filtered.gtonly.minDP0.minAC1.bcf --hardy gz --out {chr_curr}_plink_hwe --allow-extra-chr


library("stringr")
library("data.table")
library("dplyr")

chrs <- paste0("chr", 1:22)

for(i in 1:length(chrs)){
  print(i)
  chr_curr <- chrs[i]
  df_dup <- fread(str_glue("{chr_curr}.positions.gz"), header = FALSE, data.table = FALSE)
  df_hwe <- fread(str_glue("{chr_curr}_plink_hwe.hwe.gz"), header = TRUE, data.table = FALSE)
  
  #Get variants failing hwe:
  idx_failHwe <- which(df_hwe$P < 10^-6)
  hwe_to_output <- df_dup[idx_failHwe,]
  
  #Now get duplicated positions:
  df_subset <- df_dup[duplicated(df_dup$V2) | duplicated(df_dup$V2, fromLast = TRUE), ]
  
  #Write out the files
  hwe_out <- str_glue("InPSYght_hwe/{chr_curr}.InPSYght.failHWE.vars")
  fwrite(hwe_to_output, file = hwe_out, sep = "\t", row.names = FALSE, col.names = FALSE, 
         quote = FALSE)
  
  dup_out <- str_glue("InPSYght_duplicated/{chr_curr}.InPSYght.duplicated.vars")
  fwrite(df_subset, file = dup_out, sep = "\t", row.names = FALSE, col.names = FALSE, 
         quote = FALSE)
  
}

#Note that both InPSYght_duplicated/{chr}.InPSYght.duplicated.vars and InPSYght_hwe/{chr}.InPSYght.failHWE.vars
#should have the variant position in the second column. The other columns are not used
