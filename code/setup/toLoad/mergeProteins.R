#################
## merge gene nodes
#################

mergeProteins <- function(){
  
  tic("Merging protein info")
  
  #####################
  ## ADD GENE SYMBOLS
  #####################
  
  query <- "
  match (n:protein) set n:metaprotein
  "
  runQuery(query)
  
  query <- "
  match (n:uniprot) set n:metaprotein
  "
  runQuery(query)
  
  query <- "
  CREATE INDEX ON :metaprotein(partition)
  "
  runQuery(query)
  
  query <- "
  CALL gds.wcc.write({nodeQuery: 'MATCH (n:metaprotein) RETURN id(n) as id', relationshipQuery: 'MATCH (n:metaprotein)-[r:MAP]-(m:metaprotein) RETURN id(n) AS source, id(m) AS target, type(r) as type', writeProperty: 'partition'})
  YIELD componentCount, nodePropertiesWritten, computeMillis, writeMillis"
  runQuery(query)

  
  query <- "
  MATCH (a:metaprotein) WITH distinct a.partition AS part, collect(a) AS nodes
  CALL apoc.refactor.mergeNodes(nodes, {properties:'combine', mergeRels: true}) YIELD node return count(node)
  "
  runQuery(query)
  
  query <- "
  match (c:metaprotein)-[r:MAP]-(c:metaprotein) delete r
  "
  runQuery(query) 
  
  query <- "
  MATCH (a:metaprotein) remove a:metaprotein
  "
  runQuery(query)
  
  query <- "
  DROP INDEX ON :metaprotein(partition)
  "
  runQuery(query) 
  
  info(log,capture.output(toc()))
}
