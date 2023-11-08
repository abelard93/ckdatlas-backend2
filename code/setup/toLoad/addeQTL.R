## expression quatitative trait loci (mQTL)
addeQTL <- function() {
  
  ## CIS
  for (file in getFiles("QTLs/eQTL/cis")) {
    #read info on dataset
    source <- getInfo(file, c("sid", "file_name"), err)
    #skip file if no/incomplete info
    if (!all(is.na(source))) {
      tic(paste0("Added eQTL file ", source["file_name"]))
      #check/add source node for dataset
      checkSource(source, log)
      
      ## ADD edgeDEG NODE
      # sid is source_id
      
      query <- paste0("CALL apoc.periodic.iterate(
                       \"LOAD CSV WITH HEADERS FROM 'file:///'+$properties AS props return distinct props.gene_id as g\",
                      \"MERGE (ge:gene:ensembl {sid:g, ensembl_id:g}) ON CREATE SET ge.source = 'eQTL'\", {batchSize:50, parallel:true, retries:5,iterateList:true, concurrency:8, params: {properties:$properties}})
                      ")
      runQuery(query, param = TRUE, param_list = file)
     
      
      query <- paste("
                     USING PERIODIC COMMIT 200
                     LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS props 
                     MATCH (s:source {sid:$properties.source})
                     MATCH (snp:SNP {sid:props.rs_id_dbSNP151_GRCh38p7})
                     MATCH (g:gene {sid:props.gene_id}) 
                     CREATE (n:eQTL:cis {
                                       pvalue:toFloat(props.pvalue_nominal),
                                       pvalue_beta:COALESCE(toFloat(props.pvalue_beta),null),
                                       beta:COALESCE(toFloat(props.beta),null),
                                       pvalue_fdr:COALESCE(toFloat(props.fdr),null),
                                       variant_id:COALESCE(props.variant_id,null), 
                                       gene_id_version:COALESCE(props.gene_id_original,null),
                                       pvalue_cutoff:COALESCE(props.pvalue_nominal_threshold,null),
                                       maf:COALESCE(props.maf,null)}) 
                     CREATE (n)-[:FROM]->(s)
                     CREATE (snp)-[:HIT]->(n)-[:HIT]->(g)
                     MERGE (g)-[:MAP]->(snp)
                     ",sep="")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=file,source=source["sid"]$sid))

      info(log, capture.output(toc()))
      
    }
  }
  
  ## TRANS
  for (file in getFiles("QTLs/eQTL/trans")) {
    #read info on dataset
    source <- getInfo(file, c("sid", "file_name"), err)
    #skip file if no/incomplete info
    if (!all(is.na(source))) {
      tic(paste0("Added eQTL file ", source["file_name"]))
      #check/add source node for dataset
      checkSource(source, log)
      
      ## ADD edgeDEG NODE
      # sid is source_id
      
      query <- paste0("CALL apoc.periodic.iterate(
                      \"LOAD CSV WITH HEADERS FROM 'file:///'+$properties AS props return distinct props.gene_id as g\",
                      \"MERGE (ge:gene:ensembl {sid:g, ensembl_id:g}) ON CREATE SET ge.source = 'eQTL'\", {batchSize:500, parallel:true, retries:5,iterateList:true, concurrency:4, params: {properties:$properties}})
                      ")
      runQuery(query, param = TRUE, param_list = file)
      
      query <- paste("
                     USING PERIODIC COMMIT 200
                     LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file AS props 
                     MATCH (s:source {sid:$properties.source})
                     MATCH (snp:SNP {sid:props.rs_id_dbSNP151_GRCh38p7})
                     MATCH (g:gene {sid:props.gene_id}) 
                     CREATE (n:eQTL:trans {
                                       pvalue:toFloat(props.pvalue_nominal),
                                       pvalue_fdr:toFloat(props.fdr),
                                       variant_id:COALESCE(props.variant_id,null), 
                                       gene_id_version:COALESCE(props.gene_id_original,null),
                                       chromosome:COALESCE(props.gene_chr,null),
                                       tissue_af:COALESCE(props.tissue_af,null)}) 
                     CREATE (n)-[:FROM]->(s)
                     CREATE (snp)-[:HIT]->(n)-[:HIT]->(g)
                     MERGE (g)-[:MAP]->(snp)
                     ",sep="")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=file,source=source["sid"]$sid))
      
      info(log, capture.output(toc()))
    
    }
  }
}
