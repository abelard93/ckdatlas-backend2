## loop through files and make sure all SNPS are in DB 

gatherSNPs <- function() {
  
    snps <- c()
    for (file in getFiles("QTLs/")) {
      
      snps <- unique(c(snps, data.table::fread(file) %>% select(one_of("rs_id_dbSNP151_GRCh38p7","rsID", "SNP")) %>% unlist %>% as.vector %>% unique))
      message(paste0("Reading done: ",file))
    }
  
    query <- paste("CALL apoc.periodic.iterate(
                   \"UNWIND {properties} AS prop return prop as snp\",
                   \"MERGE (:SNP {sid:snp})\", {batchSize:500, parallel:true, iterateList:true, concurrency:8, params: {properties:{properties}}})
                   ",sep="")
    cypher(db,query, properties=snps)

}