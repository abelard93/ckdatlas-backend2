#################
## SNP - dummy (toDO)
#################
library(data.table)

addSNP <- function(){

      ## TO DO:: ADD SOURCE (meta file)!!
      file <- "SNP/snp_2_gene_map.csv.gz" 
      source <- basename(file) # temporary
      
      tic(paste0("Added SNPs ","QTL files"))
      
      ## all snps from SNiPA mapping
      df <- data.table::fread("SNP/snp_ids.csv")
      snps <- c()
      
      ## get all snps that are in QTL files 
      for (file in getFiles("QTLs/")) {
        snps <- unique(c(snps, data.table::fread(file) %>% select(one_of("rs_id_dbSNP151_GRCh38p7","rsID", "SNP")) %>% unlist %>% as.vector %>% unique))
        message(paste0("Reading done: ",file))
      }
      
      
      snps <- data.frame(sid=snps)
      snps <- snps$sid[!snps$sid%in%df$sid]
      rm(df)
      
      
      ## ADD SNPs without info from SNiPA
      query <- paste("CALL apoc.periodic.iterate(
                    \"UNWIND $properties AS prop return prop as snp\",
                     \"CREATE (:SNP {sid:snp, source:'QTL'})\", {batchSize:1000, parallel:true, iterateList:true, concurrency:8, params: {properties:$properties}})
                     ",sep="")
      
      runQuery(query, param = TRUE, periodic = TRUE, param_list = snps)

      
      info(log,capture.output(toc()))
      
      tic(paste0("Added SNPs ",source))
      
      file <- "SNP/snp_2_gene_map.csv.gz" 
      
      # add SNPs from SNiPA
      query <- paste("CALL apoc.periodic.iterate(
                     \"LOAD CSV WITH HEADERS FROM 'file:///SNP/snp_ids.csv' AS props with distinct props.sid as snp RETURN snp\",
                     \"CREATE (n:SNP {sid:snp,source:'SNiPA'})\", {batchSize:1000, parallel:true, iterateList:true, concurrency:8})
                     ",sep="") 
      
      runQuery(query, periodic = TRUE)
      
      # make sure genes from SNiPA are all in DB
      query <- paste("CALL apoc.periodic.iterate(
                     \"LOAD CSV WITH HEADERS FROM 'file:///'+$properties AS props with distinct props.gene as gene RETURN gene\",
                     \"MERGE (g:gene:ensembl {sid:gene}) ON CREATE SET g.source='SNiPA', g.ensembl_id=gene\", {batchSize:1000, parallel:true, iterateList:true, concurrency:8,params: {properties:$properties}})
                     ",sep="")
      
      runQuery(query, param = TRUE, periodic = TRUE, param_list = file)
      
      # add snp to gene mapping from SNiPA
      query <- paste("USING PERIODIC COMMIT 1000
                      LOAD CSV WITH HEADERS FROM 'file:///'+$properties AS props with props.sid as snp, props.gene as gene
                      MATCH (a:gene:ensembl {sid:gene}) 
                      MATCH (b:SNP {sid:snp}) 
                      CREATE (a)-[:MAP]->(b)
                     ",sep="")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = file)

      info(log,capture.output(toc()))

}