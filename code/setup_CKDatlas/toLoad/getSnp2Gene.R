library(doParallel)
library(foreach)
library(reticulate)

getSnp2gene <- function(cores=10){
  
  tic("Extracting SNP to gene mapping from db")
  
  if(!dir.exists("CKDatlas/tmp/snp_to_gene_tmp")){
    dir.create("CKDatlas/tmp/snp_to_gene_tmp")
  }
  
  detectCores()
  myCluster <- makeCluster(10, # number of cores to use
                           type = "FORK") # type of cluster
  
  clusterEvalQ(myCluster, {
    #uri = "neo4j://metabo-interactive:7687"
    uri = "bolt://10.0.0.114:7687" # address of neo4j instance
    neo4j <- import(module = "neo4j",as = "neo4j")
    token <- neo4j$basic_auth("neo4j","sysdiab")
    driver <- neo4j$GraphDatabase$driver(uri, auth=token)
  })
  
  registerDoParallel(myCluster)
  
  genes <- "MATCH (g:proteinCoding) where exists((:SNP)<-[:MAP]-(g)) RETURN g.ensembl_id as gene" %>% runQuery(read = TRUE) %>% .$gene
  
  total <-c()
  r<-c()
  counter=0
  for(i in split(1:length(genes), 1:length(genes)%/%70)){
    genes_test <- genes[i]
    counter=counter+1
    r <- foreach(gene = genes_test, .combine='rbind', .noexport = "driver") %dopar% { 
      
      "MATCH (g:proteinCoding) where g.ensembl_id=$properties with g MATCH (g)-[:MAP]->(s:SNP) return distinct g.ensembl_id as gene, s.sid as sid" %>% runQuery(read = TRUE, param = TRUE, param_list = gene)
    }
    fwrite(r,file = paste0("CKDatlas/tmp/snp_to_gene_tmp/snp_to_gene_",counter,".csv"))
    
  }
  
  stopCluster(myCluster)
  
  ## READ file
  snp2gene = lapply(getFiles(dir="CKDatlas/tmp/snp_to_gene_tmp/", pat="snp_to_gene_*"), fread, header=T)
  snp2gene = do.call(rbind,snp2gene)
  fwrite(snp2gene,file = paste0("CKDatlas/tmp/snp_to_gene_mapping.csv"))
  file.remove(getFiles(dir="CKDatlas/tmp/snp_to_gene_tmp/", pat="snp_to_gene_[0-9]*.csv"))
  
  info(log,capture.output(toc()))
  return (snp2gene)
}