#################
## genes
#################

addGeneSymbol <- function(){
  
  tic("Added gene symbols info")
  
  #####################
  ## ADD GENE SYMBOLS
  #####################
  query <- "
  USING PERIODIC COMMIT 1000 
  LOAD CSV WITH HEADERS FROM 'file:///geneSymbol/hgnc_ensembl.csv' AS line 
  MATCH (g:gene:ensembl {sid:line.sid}) 
  MERGE (g)-[:SYNONYM]->(n:geneSymbol {sid:line.symbol, symbol_id:line.symbol})
  ON CREATE SET n.symbol_source='hgnc'
  "
  runQuery(query, periodic = TRUE)
  
  info(log,capture.output(toc()))
}
