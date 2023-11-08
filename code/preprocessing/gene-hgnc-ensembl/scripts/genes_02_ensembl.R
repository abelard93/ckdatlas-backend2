require(tidyverse)
require(biomaRt)
##
##
## PARAMS
args = commandArgs(trailingOnly = T)
file.out <- args[1]
version <- args[2]
##
## file.in <- "../raw/genes/hgnc_complete_set.txt"
##
################################################################################
## READ
################################################################################

library(biomaRt)
if(is.na(version)){
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
}else{
  archive <- listEnsemblArchives() %>% filter(version == !!version) %>% .$url
  ensembl = useMart(host= archive,biomart= "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
}

#attributes = listAttributes(ensembl)

gene2transcript2protein=getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","ensembl_peptide_id" ), mart=ensembl) %>% 
  mutate_all(~na_if(.,""))
ensembl2entrez=getBM(attributes=c("ensembl_gene_id","entrezgene" ), mart=ensembl) %>% 
  mutate_all(~na_if(.,""))
gene=getBM(attributes=c("ensembl_gene_id","description","band","chromosome_name","start_position","end_position","strand","gene_biotype", "hgnc_symbol","external_gene_name"), mart=ensembl) %>% 
  mutate_all(~na_if(.,""))
transcript=getBM(attributes=c("ensembl_transcript_id","transcription_start_site","transcript_length","transcript_biotype","transcript_start","transcript_end","band", "chromosome_name","external_transcript_name"), mart=ensembl) %>% 
  mutate_all(~na_if(.,""))
protein=getBM(attributes=c("ensembl_peptide_id","uniprotswissprot","uniprotsptrembl"), mart=ensembl) %>% 
  mutate_all(~na_if(.,""))

################################################################################
## WRITE
################################################################################
write.csv(x = gene2transcript2protein, file = paste0(file.out,substitute(gene2transcript2protein),".csv"), quote = T,na="",row.names = F)
write.csv(x = ensembl2entrez, file = paste0(file.out,substitute(ensembl2entrez),".csv"), quote = T,na="",row.names = F)
write.csv(x = gene, file = paste0(file.out,substitute(gene),".csv"), quote = T,na="",row.names = F)
write.csv(x = transcript, file = paste0(file.out,substitute(transcript),".csv"), quote = T,na="",row.names = F)
write.csv(x = protein, file = paste0(file.out,substitute(protein),".csv"), quote = T,na="",row.names = F)



