################
## merge gene nodes
#################

mergeGenes <- function(){
  
  tic("Merging gene info")
  
  #####################
  ## ADD GENE SYMBOLS
  #####################
  
  query <- "
  match (n:geneSymbol) set n:metagene
  "
  runQuery(query)
  
  query <- "
  match (n:gene) set n:metagene
  "
  runQuery(query)
  
  query <- "
  CREATE INDEX ON :metagene(partition)
  "
  runQuery(query)
  
  query <- "
  CALL gds.wcc.write({nodeQuery: 'MATCH (n:metagene) RETURN id(n) as id', relationshipQuery: 'MATCH (n:metagene)-[r:SYNONYM]-(m:metagene) RETURN id(n) AS source, id(m) AS target, type(r) as type', writeProperty: 'partition'})
  YIELD componentCount, nodePropertiesWritten, computeMillis, writeMillis"
  runQuery(query)

  
  query <- "
  MATCH (a:metagene) WITH distinct a.partition AS part, collect(a) AS nodes
  CALL apoc.refactor.mergeNodes(nodes, {properties:'combine', mergeRels: true}) YIELD node return count(node)
  "
  runQuery(query)
  
  query <- "
  match (c:metagene)-[r:SYNONYM]-(c:metagene) delete r
  "
  runQuery(query)  
  
  
  query <- "
  MATCH (a:metagene) remove a:metagene
  "
  runQuery(query) 

  query <- "
  DROP INDEX ON :metagene(partition)
  "
  runQuery(query) 
  
  info(log,capture.output(toc()))
}
