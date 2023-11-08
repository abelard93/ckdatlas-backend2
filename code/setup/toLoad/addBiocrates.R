#' Add metabolites measured on a Biocrates platform 
#'
#' @description
#' Reads annotations of Biocrates metabolites.
#' 1 - Add :metBiocrates:metMeasured nodes
#' 2 - Create pseudo ontology layer (using chemical id)
#' CAUTION: only implemented to read data from meta file for now - different from Metabolon where results files can be added.

addBiocrates <- function(){

  for(file in getFiles("metBiocrates","_meta")){
   
    ## add biocrates nodes
    query <- paste0("
    LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS line with line 
    MERGE (n:metBiocrates:metMeasured {sid:apoc.coll.toSet(split(line.sid,'|'))})
    with n, line 
    where exists(line.platform)
    call apoc.create.addLabels(id(n), split(line.platform,'|')) yield node 
    return count(*)
    ")
    runQuery(query)
    
    ## add Ontology
    query <- paste0("
    LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS line with line 
    MATCH (n:metBiocrates {sid:apoc.coll.toSet(split(line.sid,'|'))})
    MERGE (m:metOntology {sid:n.sid[0]})
    MERGE (n)-[:MAP]->(m) with m, line 
    where exists(line.metHMDB) 
    with m, split(line.metHMDB, '|') as ids unwind ids as id
    MATCH (h1:accHMDB)-[]-(h:metHMDB) where h1.sid=id 
    MERGE (m)-[:MAP]->(h)
    ")
    runQuery(query)
    
    ## add biochemical names (as nodes)
    query <- paste0("
    LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS line with line 
    where exists(line.metName)
    MATCH (n:metBiocrates {sid:apoc.coll.toSet(split(line.sid,'|'))})
    with n, split(line.metName, '|') as names unwind names as name
    MERGE (m:metName {sid:name})
    MERGE (n)-[:HAS_NAME]->(m)
    ")
    runQuery(query)
    
    ## set name property of biocrates node to sid if no name given
    query <- paste0("
                    LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS line with line 
                    where exists(line.sid)
                    MATCH (n:metBiocrates {sid:apoc.coll.toSet(split(line.sid,'|'))})
                    with n
                    where not exists(n.name)
                    SET n.name=n.sid[0]
                    ")
    runQuery(query)

    ## set name property of biocrates node
    query <- paste0("
                    LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS line with line 
                    MATCH (n:metBiocrates {sid:apoc.coll.toSet(split(line.sid,'|'))})
                    with n, line
                    where exists(line.name)
                    SET n.name=line.name
                    with n, line.name as name
                    MERGE (m:metName {sid:name})
                    MERGE (n)-[:HAS_NAME]->(m)
                    ")
    runQuery(query)
    
    ## IDS
    for(i in c("metKEGG","metCAS","metIUPAC","metSmiles")){
      query <- paste0("
      LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS line with line 
      where exists(line.",i,")
      MATCH (n:metBiocrates {sid:apoc.coll.toSet(split(line.sid,'|'))})
      with n, split(line.",i,", '|') as ids unwind ids as id
      MERGE (m:",i," {sid:id})
      MERGE (n)-[:HAS_ID]->(m)
      ")
      runQuery(query)
    }
    
    ## classification
    query <- paste0("
    LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS line with line 
    where exists(line.metClassification)
    MATCH (n:metBiocrates {sid:apoc.coll.toSet(split(line.sid,'|'))})
    with n, split(line.metClassification, '|') as names unwind names as name
    MERGE (m:metClassification {sid:name})
    MERGE (n)-[:CLASSIFIED_AS]->(m)
    ")
    runQuery(query)
    
    ## superPathway
    query <- paste0("
                    LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS line with line 
                    where exists(line.metSuperPathway)
                    MATCH (n:metSuperPathway {sid:line.metSuperPathway})
                    with n, line.metSuperPathway as name
                    MERGE (m:metSuperPathway {sid:name})
                    MERGE (n)-[:IN_PATHWAY]->(m)
                    ")
    runQuery(query)
    
    ## subPathway
    query <- paste0("
                    LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS line with line 
                    where exists(line.metSubPathway)
                    MATCH (n:metSubPathway {sid:line.metSubPathway})
                    with n, line.metSubPathway as name
                    MERGE (m:metSubPathway {sid:name})
                    MERGE (n)-[:IN_PATHWAY]->(m)
                    ")
    runQuery(query)
    
    ## weight
    query <- paste0("
    LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS line with line 
    where exists(line.metWeight)
    MATCH (n:metBiocrates {sid:apoc.coll.toSet(split(line.sid,'|'))})
    with n, split(line.metWeight, '|') as ids unwind ids as id
    MERGE (m:metWeight {sid:toFloat(id)})
    MERGE (n)-[:HAS_WEIGHT]->(m)
    ")
    runQuery(query)
    
    ## formula
    query <- paste0("
    LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS line with line 
    where exists(line.metFormula)
    MATCH (n:metBiocrates {sid:apoc.coll.toSet(split(line.sid,'|'))})
    MERGE (m:metFormula {sid:line.metFormula})
    MERGE (n)-[:HAS_FORMULA]->(m)
    ")
    runQuery(query)
    
    ## HMDB
    query <- paste0("
    LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS line with line 
    where exists(line.metHMDB)
    MATCH (n:metBiocrates {sid:apoc.coll.toSet(split(line.sid,'|'))})
    with n, split(line.metHMDB, '|') as ids unwind ids as id
    MATCH (h1:accHMDB)-[]-(h:metHMDB) where h1.sid=id 
    MERGE (n)-[:HAS_ID]->(h) 
    ")
    runQuery(query)
    
    ## ChEBI
    query <- paste0("
    LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS line with line 
    where exists(line.metChEBI)
    MATCH (n:metBiocrates {sid:apoc.coll.toSet(split(line.sid,'|'))})
    with n, split(line.metChEBI, '|') as ids unwind ids as id
    MATCH (h1:accChEBI)-[]-(h:metChEBI) where h1.sid=id 
    MERGE (n)-[:HAS_ID]->(h) 
    ")
    runQuery(query)
    
    ## add mapping
    query <- paste0("
    LOAD CSV WITH HEADERS FROM 'file:///",file,"' AS line with line 
    where exists(line.metMetabolon)   
    MATCH (n:metBiocrates {sid:apoc.coll.toSet(split(line.sid,'|'))})
    with n, split(line.metMetabolon, '|') as ids unwind ids as id
    MATCH (m:metMetabolon {COMP_ID:id})
    MERGE (n)-[:MAP {info:'manual mapping'}]->(m)
    ")
    runQuery(query)
  }
}
