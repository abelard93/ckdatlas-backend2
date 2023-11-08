#################################################
##  Script to format GTEx v8 for Atlas
##
##  HowTo: only need to adjust base directory
##
##  Author: Maria Wörheide
##  Date: 2020-06-16
#################################################


zap()

library(tidyverse)
library(data.table)


## SET BASE + RESULTS DIRECTORY --------------------------------------
base_dir <- "/storage/metabostore/project/ADatlas/repo/code/GTEx/"
result_dir <- "/storage/metabostore/project/ADatlas/repo/data/QTLs/eQTL/"

setwd(base_dir)

# make & go to folder
system("mkdir data")
setwd(paste0(base_dir,"data"))

## DOWNLOAD DATA ------------------------------------------

# download cis eQTLs
download.file("https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar",destfile = "GTEx_v8_cis.tar")
system("tar -xvf GTEx_v8_cis.tar")
system("rm GTEx_v8_cis.tar")
system("rm GTEx_Analysis_v8_eQTL/*v8.egenes.txt.gz")

# download trans eQTLs
download.file("https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_trans_eGenes_fdr05.txt",destfile = "GTEx_v8_trans.txt")

# download lookup table variant id to rsID
download.file("https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz",destfile = "GTEx_v8_lookup.txt.gz")



## LOOKUP TABLES ------------------------------------------
lookup <- data.table::fread("GTEx_v8_lookup.txt.gz")

# variant ids
lookup_vec <- setNames(object = lookup$rs_id_dbSNP151_GRCh38p7,nm = lookup$variant_id )



## CIS ------------------------------------------


# 
# for(file in list.files("GTEx_Analysis_v8_eQTL/",full.names = T)){
#   name <- basename(file) %>% str_replace(".gz","")
#   eqtl <- data.table::fread(file)
#   eqtl$rs_id_dbSNP151_GRCh38p7 <-  lookup_vec[eqtl$variant_id]
#   fwrite(eqtl,file = paste0(result_dir,"cis/",name),sep = ",",col.names = T,row.names = F,na="")
#   message(paste0("Finished processing file ", name))
# }

for(file in list.files("GTEx_Analysis_v8_eQTL/",full.names = T)){
  ## get tissue name and sid from file
  # TISSUE
  tissue <- basename(file) %>% sub(.,pattern = ".v8.signif_variant_gene_pairs.txt.gz",replacement = "") %>%
    sub(.,pattern = "_",replacement = " ")
  # SID
  id <- basename(file) %>% sub(.,pattern = ".v8.signif_variant_gene_pairs.txt.gz",replacement = " GTEx v8 ") %>% 
    sub(.,pattern = "_",replacement = " ")
  # OUTPUT
  output <- basename(file) %>% 
    str_replace_all(.,"[.]", "_") %>% 
    sub(.,pattern = "_txt_gz",replacement = ".csv")
  ## READ GTEx FILE
  gtex.file <- fread(file,stringsAsFactors = F) %>% 
    mutate(gene_id_original=gene_id,
           gene_id=map_chr(gene_id,~strsplit(.,"\\.")[[1]][1])) %>% 
    rename_all(~sub(x=.,"pval","pvalue")) 
  ## ADD rsID
  gtex.file$rs_id_dbSNP151_GRCh38p7 <-  lookup_vec[gtex.file$variant_id]
  gtex.file <- gtex.file %>% filter(rs_id_dbSNP151_GRCh38p7!=".")
  ## GENERATE INFO 
  gtex.info <- data.frame(
    c("sid","file_name","organism","sample_type","version","publication","project","status","study_type"), 
    c(id,output,"Homo Sapiens",tissue,"v8","GTExConsortium_2019","Genotype-Tissue Expression (GTEx)","healthy","population-based"))
  ## WRITE
  write.table(gtex.file,file = paste0(result_dir,"cis/",output),sep = ",",col.names = T,row.names = F)
  write.table(gtex.info,file = paste0(result_dir,"cis/",str_replace(output,".csv","_info.csv")),sep = ",",col.names = F,row.names = F)
  message(paste0("Finished processing file ", tissue))
}



## TRANS ------------------------------------------

# file <- "/home/icb/maria.woerheide/Documents/GTEx_v8/GTEx_Analysis_v8_trans_eGenes_fdr05.txt"
# name <- basename(file) 
# eqtl <- data.table::fread(file)
# eqtl$rs_id_dbSNP151_GRCh38p7 <-  lookup_vec[eqtl$variant_id]
# fwrite(eqtl,file = paste0(result_dir,name),sep = ",",col.names = T,row.names = F,na="")
# message(paste0("Finished processing file ", name))
# 

  
## READ GTEx FILE
gtex.file <- fread("GTEx_v8_trans.txt",stringsAsFactors = F) %>% 
  # filter(qval<=0.05) %>% 
  mutate(gene_id_original=gene_id,gene_id=map_chr(gene_id,~strsplit(.,"\\.")[[1]][1])) %>% 
  rename_all(~sub(x=.,"pval","pvalue")) 
## ADD rsID
gtex.file$rs_id_dbSNP151_GRCh38p7 <-  lookup_vec[gtex.file$variant_id]
gtex.file <- gtex.file %>% filter(rs_id_dbSNP151_GRCh38p7!=".")

for(tissue in unique(gtex.file$tissue_id)){
  
  gtex.tissue <- gtex.file %>% filter(tissue_id==tissue)
  output <- basename("GTEx_v8_trans.txt") %>% 
    str_replace_all(.,"[.]", "_") %>% 
    sub(.,pattern = "_txt",replacement = ".csv") %>% 
    paste0(tissue,"_",.)
  ## get tissue name and sid from file
  # TISSUE
  tissue <- tissue %>%
    sub(.,pattern = "_",replacement = " ")
  # SID
  id <- paste0(tissue," GTEx v8 trans")
  # OUTPUT
  
  ## GENERATE INFO 
  gtex.info <- data.frame(
    c("sid","file_name","organism","sample_type","version","publication","project","status","study_type"), 
    c(id,output,"Homo Sapiens",tissue,"v8","GTExConsortium_2019","Genotype-Tissue Expression (GTEx)","healthy","population-based"))
  ## WRITE
  write.table(gtex.tissue,file = paste0(result_dir,"trans/",output),sep = ",",col.names = T,row.names = F)
  write.table(gtex.info,file = paste0(result_dir,"trans/",str_replace(output,".csv","_info.csv")),sep = ",",col.names = F,row.names = F)
  
}
  




########### LEGACY ----------------------------------------
########### Version 7

file.info <- "~/Documents/PHD/krumsieklab/neo4jserver/code/gTEX/data_v7/raw"
file.output <- "~/Documents/PHD/krumsieklab/neo4jserver/code/gTEX/data_v7/processed/"

for(file in list.files(file.info,full.names = T,recursive = T)){
  ## get tissue name and sid from file
  # TISSUE
  tissue <- basename(file) %>% sub(tissue,pattern = ".v7.egenes.txt.gz",replacement = "") %>%
    sub(.,pattern = "_",replacement = " ")
  # SID
  id <- basename(file) %>% sub(.,pattern = ".v7.egenes.txt.gz",replacement = " GTEx v7 ") %>% 
    sub(.,pattern = "_",replacement = " ")
  # OUTPUT
  output <- basename(file) %>% 
    str_replace_all(.,"[.]", "_") %>% 
    sub(.,pattern = "_txt_gz",replacement = ".csv")
  ## READ GTEx FILE
  gtex.file <- read.delim(file,stringsAsFactors = F) %>% 
    filter(qval<=0.05) %>% 
    mutate(gene_id_original=gene_id,gene_id=map_chr(gene_id,~strsplit(.,"\\.")[[1]][1])) %>% 
    rename_all(~sub(x=.,"pval","pvalue")) 
  ## GENERATE INFO 
  gtex.info <- data.frame(
    c("sid","file_name","organism","sample_type","version","publication","project","status","study_type"), 
    c(id,output,"Homo Sapiens",tissue,"v7","lonsdale_2013_23715323","Genotype-Tissue Expression (GTEx)","healthy","population-based"))
  ## WRITE
  write.table(gtex.file,file = paste0(file.output,output),sep = ",",col.names = T,row.names = F)
  write.table(gtex.info,file = paste0(file.output,str_replace(output,".csv","_info.csv")),sep = ",",col.names = F,row.names = F)
}

#################

file.info <- "~/Documents/PHD/krumsieklab/neo4jserver/code/gTEX/data/processed/"
file.output <- "~/Documents/PHD/krumsieklab/neo4jserver/code/gTEX/data/"
file.mapping <- fread("~/Desktop/snp_2_gene_map.csv.gz")
vector.mapping <- setNames(file.mapping$gene,file.mapping$sid)

for(file in list.files(file.info,full.names = T,recursive = T,pattern = "*gene_pairs.csv.gz")){
  
  ## get tissue name and sid from file
  # TISSUE
  tissue <- basename(file) %>% sub(.,pattern = ".v8.signif_variant_gene_pairs.csv.gz",replacement = "") %>%
    sub(.,pattern = "_",replacement = " ")
  # SID
  id <- basename(file) %>% sub(.,pattern = ".v8.signif_variant_gene_pairs.csv.gz",replacement = " GTEx v8 ") %>% 
    sub(.,pattern = "_",replacement = " ")
  # OUTPUT
  output <- basename(file) %>% 
    str_replace_all(.,"[.]", "_") %>% 
    sub(.,pattern = "signif_variant_gene_pairs_csv_gz",replacement = "summary.csv")
  ## READ GTEx FILE
  gtex.file <- fread(file,stringsAsFactors = F) %>% transmute(gene=gene_id,pvalue=pvalue_nominal,rsid=rs_id_dbSNP151_GRCh38p7)
  
  gtex.file$gene_snp <- vector.mapping[gtex.file$rsid]
  gtex.file$tissue <- tissue
  gtex.file.summary <- gtex.file %>% filter(!is.na(gene_snp)) %>% group_by(gene,gene_snp) %>% summarize(pvalue_summary=paste0("</style>
                                                                                                                              <table class='snipa-plot-tooltip'><thead><tr><th colspan='3'>SIGNIFICANT_COREG</th></tr></thead>
                                                                                                                              <tbody><tr><td>min_pvalue:&nbsp;</td><td>",min(pvalue),"</td></tr>
                                                                                                                              <tr><td>snps:&nbsp;</td><td>",n(),"</td></tr>
                                                                                                                              </tbody></table>
                                                                                                                              <table class='snipa-plot-tooltip'><tr><th colspan='3'>",unique(tissue),"</th></tr>",
                                                                                                                              paste0("<td>",rsid,"</td><td colspan='2'>",pvalue,"</td></tr>",collapse = " "),"</table>"), min_pvalue=min(pvalue))
  names(gtex.file.summary) <- c("gene","gene_snp",paste0("label_",tissue %>% tolower %>% str_replace_all(.," ","_")),paste0("min_pvalue_",tissue %>% tolower %>% str_replace_all(.," ","_")))
  
  
  
  ## WRITE
  write.table(gtex.file.summary,file = paste0(file.output,output),sep = ",",col.names = T,row.names = F)
  #write.table(gtex.info,file = paste0(file.output,str_replace(output,".csv","_info.csv")),sep = ",",col.names = F,row.names = F)
}


