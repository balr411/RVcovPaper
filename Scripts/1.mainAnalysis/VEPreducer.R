library("argparser")
library("data.table")
library("dplyr")
library("stringr")

p <- arg_parser("Program that takes VEP annoation file and reduces to only variant name, consequence, and group")

p <- add_argument(p, "--chr", help = "Chromosome")

argv <- parse_args(p)


anno_name <- paste0(argv$chr, ".vcf.vep.out.gz")
anno <- fread(cmd = paste0("zgrep -v ^## ", anno_name), h = T, data.table = F)
anno <- anno[anno$PICK == 1,]

#Find variants with deleterious consequences
cons <- unique(anno[[7]])
imp_vars <- c()
bad_cons <- c("stop_gained", "splice_acceptor_variant", "stop_lost", "splice_donor_variant", "frameshift_variant", "start_lost", "missense_variant")
for(i in 1:length(cons)){
  temp <- strsplit(cons[i], split = ",")
  if(sum(bad_cons %in% unlist(temp)) > 0){
    imp_vars <- c(imp_vars, cons[i])
  }
}

anno <- anno[anno[[7]] %in% imp_vars,]
anno_names <- c("#Uploaded_variation", "Consequence", "SYMBOL", "MutationTaster_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "SIFT4G_pred", "LRT_pred")
vep_out_sub <- anno[,anno_names]

#Now add predictions
vep_out_sub$LRT_pred <- str_count(vep_out_sub$LRT_pred, "D")
vep_out_sub$MutationTaster_pred <- str_count(vep_out_sub$MutationTaster_pred, "D") + str_count(vep_out_sub$MutationTaster_pred, "A")
vep_out_sub$Polyphen2_HDIV_pred <- str_count(vep_out_sub$Polyphen2_HDIV_pred, "D") + str_count(vep_out_sub$Polyphen2_HDIV_pred, "P")
vep_out_sub$Polyphen2_HVAR_pred <- str_count(vep_out_sub$Polyphen2_HVAR_pred, "D") + str_count(vep_out_sub$Polyphen2_HVAR_pred, "P")
vep_out_sub$SIFT4G_pred <- str_count(vep_out_sub$SIFT4G_pred, "D")
vep_out_sub$LRT_pred <- ifelse(vep_out_sub$LRT_pred > 0, 1, 0)
vep_out_sub$MutationTaster_pred <- ifelse(vep_out_sub$MutationTaster_pred > 0, 1, 0)
vep_out_sub$Polyphen2_HDIV_pred <- ifelse(vep_out_sub$Polyphen2_HDIV_pred > 0, 1, 0)
vep_out_sub$Polyphen2_HVAR_pred <- ifelse(vep_out_sub$Polyphen2_HVAR_pred > 0, 1, 0)
vep_out_sub$SIFT4G_pred <- ifelse(vep_out_sub$SIFT4G_pred > 0, 1, 0)
vep_out_sub$Sum_Algs <- rowSums(vep_out_sub[,c(4:8)])
vep_out_sub$Sum_Algs <- ifelse(vep_out_sub$Sum_Algs == 5, 1, 0)

table_name <- paste0(argv$chr, ".vcf.vep.reduced")
write.table(vep_out_sub, table_name)