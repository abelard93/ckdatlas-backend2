#################
## BIOCRATES
#################

## only data from meta file for now

addMetaMetabolite <- function(){

  
  tic("Created metMeta")
  
 

  # merge/map metabolites that have the same COMP_ID (even if they were measured on diferent platforms)
  ## COMP_ID (metabolon)
  "Call apoc.periodic.iterate(\"match (a:metMetabolon),(b:metMetabolon) where id(a)>id(b) and a.COMP_ID=b.COMP_ID return a,b\",\"MERGE (a)-[:MAP]->(b) return count(*)\",{parallel:false,batchSize:500,iterateList:true})" %>% runQuery(periodic = TRUE)
  
  "Call apoc.periodic.iterate(\"match (a:metMetabolon)-[:MAP]-(c:metOntology)-[:MAP]-(b:metMetabolon) where id(a)>id(b) and a.name=b.name and not exists((a)-[:MAP]-(b)) return a,b\",\"MERGE (a)-[:MAP]->(b) return count(*)\",{parallel:false,batchSize:500,iterateList:true})" %>% runQuery(periodic = TRUE)
  
  
  # need to drop constraits to duplicate nodes
  "DROP CONSTRAINT ON (m:metMetabolon) ASSERT m.sid IS UNIQUE" %>% runQuery()
  "CREATE INDEX ON :metMetabolon(sid)" %>% runQuery()
  "CREATE INDEX ON :metMeasured(sid)" %>% runQuery()
  "CREATE INDEX ON :metBiocrates(sid)" %>% runQuery()
 # "DROP CONSTRAINT ON (m:metMeasured) ASSERT m.sid IS UNIQUE" %>% cypher(db,.)
 
  ##
  # duplicate all metabolite nodes
  # give duplicated new label metMeta
  # remove other labels on duplicated
  

  # find nodes that should be merged later on
  "CALL gds.wcc.write({nodeQuery: 'MATCH (n:metMeasured) RETURN id(n) as id', relationshipQuery: 'MATCH (n:metMeasured)-[r:MAP]-(m:metMeasured) RETURN id(n) AS source, id(m) AS target, type(r) as type', writeProperty: 'partition'})
  YIELD componentCount, nodePropertiesWritten, computeMillis, writeMillis" %>% runQuery()
  
  #mark nodes that need to be copied
  "match (n:metMeasured) WHERE (n)-[:HIT]-() or (n)-[:GGM]-() set n.copied=\"false\"" %>% runQuery()

  # copy nodes, copy relationships, add label metMeta, remove other labels 
  "Call apoc.periodic.iterate(
    \"match (n:metMeasured)  WHERE exists( n.copied ) return distinct n\",
    \"REMOVE n.copied with n 
    CALL apoc.refactor.cloneNodesWithRelationships([n]) YIELD output
    CALL apoc.create.addLabels(output, ['metMeta']) YIELD node as no
    CALL apoc.create.removeLabels(no,['metMeasured','metMetabolon','metBiocrates','p150','p180','bileAcidKit']) YIELD node
    RETURN count(*)\", {parallel:false,batchSize:10,iterateList:true})" %>% runQuery(periodic = TRUE)

  # add constraint again
  "DROP INDEX ON :metMetabolon(sid)" %>% runQuery()
  "CREATE CONSTRAINT ON (m:metMetabolon) ASSERT m.sid IS UNIQUE" %>% runQuery()
  
  # add index on new metMeta
  "CREATE INDEX ON :metMeta(sid)" %>% runQuery()
  "CREATE INDEX ON :metMeta(partition)" %>% runQuery()

  # merge nodes connected by map
  # https://github.com/neo4j-contrib/neo4j-apoc-procedures/issues/1408
  "Call apoc.periodic.iterate(
    \"MATCH (a:metMeta) WITH distinct a.partition AS part, collect(a) AS nodes RETURN nodes\",
    \"CALL apoc.nodes.get(nodes) YIELD node WITH collect(node) AS nodes2 CALL apoc.refactor.mergeNodes(nodes2, {properties:'combine', mergeRels: true}) YIELD node return count(*)\", {batchSize:1})" %>% runQuery(periodic = TRUE)
  
  

  # remove self loops
  "MATCH (a:metMeta)-[r:MAP]-(a:metMeta) delete r" %>% runQuery()
  
  
  # all names and sid properties to list (to avoid errors when doing list operations later on)
  "MATCH (n:metMeta) where apoc.meta.type(n.name)='STRING' set n.name=[n.name]" %>% runQuery()

  "MATCH (n:metMeta) where apoc.meta.type(n.sid)='STRING' set n.sid=[n.sid]" %>% runQuery()

  # merge nodes that are connected to the same unknown metabolite (metabolon X-**)
  # reason: when identified metabolon ids get new COMP id => same metabolite with dif COMP_IDS
  # https://github.com/neo4j-contrib/neo4j-apoc-procedures/issues/1408
  "Call apoc.periodic.iterate(\"match (a:metMeta)-[:HAS_NAME]-(c:metName) where replace(c.sid,' ','') STARTS WITH 'X-' with distinct replace(c.sid,' ','') as unkown,  collect(distinct a) as metabo where size(metabo) >1 return metabo\" ,
  \"CALL apoc.nodes.get(metabo) YIELD node WITH collect(node) AS nodes2 CALL apoc.refactor.mergeNodes(nodes2, {properties:'combine', mergeRels: true}) YIELD node return count(*)\",
  {batchSize:1})" %>% runQuery(periodic = TRUE)
  
  ## order name of metabolite such that shortest name first 
  "match (m:metMeta) 
   WITH m.name AS name , m 
   UNWIND name AS n 
   WITH DISTINCT m, n 
   ORDER BY size(n) 
   WITH m, collect(n) AS list 
   SET m.name=list" %>% runQuery()
  
  # give thos metabolites which have been identified names 
  "match (a:metMeta)-[:HAS_NAME]-(c:metName) where replace(c.sid,' ','') STARTS WITH 'X-' and size(a.name)>1 with distinct a set a.name=[last(a.name)]" %>% runQuery()
  
  # add constraint on metMeta
  "DROP INDEX ON :metMeta(sid)" %>% runQuery()
  "CREATE CONSTRAINT ON (m:metMeta) ASSERT m.sid IS UNIQUE" %>% runQuery()
  

  info(log,capture.output(toc()))

}
