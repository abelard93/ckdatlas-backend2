library(biomaRt)
tmp_path="/Users/mwoerheide/Documents/ICB/neo4j/data/genes/"

##get gene list
ensembl=useMart("ensembl")  # using ensembl database data
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)   # from ensembl using homosapien gene data
genemap=getBM(attributes=c("entrezgene","hgnc_symbol","ensembl_gene_id" ), mart=ensembl) # fuction to get  gene id's and gene name from data base

write.csv(genemap,file=paste(tmp_path,"biomaRt_hsa.csv",sep = ""),na = "",row.names = F)

##mouse
ensembl=useMart("ensembl")  # using ensembl database data
ensembl=useDataset("mmusculus_gene_ensembl",mart=ensembl)   # from ensembl using mmusculus gene data
genemap=getBM(attributes=c("entrezgene","mgi_symbol","ensembl_gene_id" ), mart=ensembl) # fuction to get  gene id's and gene name from data base

write.csv(genemap,file=paste(tmp_path,"biomaRt_mmu.csv",sep = ""),na = "",row.names = F)
