## Sample code
# - shows the differences between neo4r and RNeo4j
# - how to pass properties in RNeo4j
# - how to visualize a graph with neo4r

############
## RNeo4j ##
############

library(RNeo4j)

#connect to database
db = startGraph("http://vmetaneo4j:7474/db/data/", username = "neo4j", password = "sysdiab")

#############################
## PROPERTIES AS VECTOR #####
#############################

## query that assumes that the gene symbol is annotaed in geBiGG node
x=c("TYR","COMT","AOC2","AOC3")

query <- "UNWIND {properties} AS prop 
          MATCH (g:geBiGG:Recon3D {symbol:prop})-[:CATALYZES]-(:reBiGG:Recon3D)-[:PARTICIPATES]-(m:metBiGG:Recon3D) 
          RETURN distinct g.symbol as gene, collect(distinct m.sid) as metabolite"
## return df
df <- cypher(db,query,properties=x)

## query that assumes that the gene symbol is not annotaed in geBiGG node
query <- "UNWIND {properties} AS prop 
          MATCH (g:gene)-[:MAP]-(:geBiGG:Recon3D)-[:CATALYZES]-(:reBiGG:Recon3D)-[:PARTICIPATES]-(m:metBiGG:Recon3D) 
          where prop in g.sid
          RETURN distinct g.symbol as gene, collect(distinct m.sid) as metabolite"
## return df
df <- cypher(db,query,properties=x)

#########################
## PROPERTIES AS DF #####
#########################

## each row holds info for one "query 
x=data.frame(sid=c("TYR","COMT","AOC2","AOC3"),additional_info=c("SIGNIFICANT","SIGNIFICANT","SIGNIFICANT","NOT SIGNIFICANT"))

query <- "UNWIND {properties} AS prop 
          MATCH (g:geBiGG:Recon3D {symbol:prop.sid})-[:CATALYZES]-(:reBiGG:Recon3D)-[:PARTICIPATES]-(m:metBiGG:Recon3D) 
          RETURN distinct g.symbol as gene, collect(distinct m.sid) as metabolite, prop.additional_info as info"
## return df
df <- cypher(db,query,properties=x)



#####################
## ALGORITHMS   #####
#####################

df <- "CALL algo.degree.stream('MATCH (g:gene) match (m:metBiocrates) with g,m limit 30 return collect(distinct id(g))+collect(distinct id(m))','MATCH (g:gene) match (m:metBiocrates) with g,m limit 30 with collect(distinct id(g))+collect(distinct id(m)) as nodes CALL apoc.algo.cover(nodes) YIELD rel return rel',{direction: \"incoming\",graph: \"cypher\", params: {new:new} })
YIELD nodeId, score
RETURN algo.asNode(nodeId).sid AS name, score AS followers
ORDER BY followers DESC" %>% 
  cypher(.,query)


# MATCH (g:proteinCoding) where g.ensembl_id in ['ENSG00000244474','ENSG00000189357','ENSG00000106346'] WITH collect(distinct g) as genes unwind genes as n
# OPTIONAL MATCH p=(n)-[:SIGNIFICANT_COEXPR]-(g1) where all(x in relationships(p) where any( y in ['DLPFC','FP'] where y IN split(x.tissue, '|')))
# with collect(distinct g1)+genes as genes2 unwind genes2 as n2
# OPTIONAL MATCH (n2)-[rel:SIGNIFICANT_QTL]-(m)
# OPTIONAL MATCH (n2)-[:SIGNIFICANT_traitQTL]-(t)
# with collect(distinct m)+collect(distinct t)+genes2 as node_subset unwind node_subset as ni with collect(distinct ID(ni)) as new, collect(distinct ni) as nodi
# CALL algo.degree.stream('with $new','CALL apoc.algo.cover($new) YIELD rel',{direction: "incoming",graph: "cypher", params: {new:new} })
# YIELD nodeId, score
# RETURN algo.asNode(nodeId).sid AS name, score AS followers
# ORDER BY followers DESC


############
## neo4r ##
############

library(neo4r)

## connect to db
con <- neo4j_api$new(
  url = "http://vmetaneo4j:7474/", 
  user = "neo4j", 
  password = "sysdiab"
)
con$ping()

## each row holds info for one "query 
x=data.frame(sid=c("TYR","COMT","AOC2","AOC3"),additional_info=c("SIGNIFICANT","SIGNIFICANT","SIGNIFICANT","NOT SIGNIFICANT"))


#### HACK ####
# pass properties as map within cypher query
properties_map <- paste0("{sid:'",x$sid,"',additional_info:'",x$additional_info,"'}",collapse ="," ) 

## QUERY
G <-paste0("UNWIND [",properties_map,"] AS prop 
     MATCH p=(:geBiGG:Recon3D {symbol:prop.sid})-[:CATALYZES]-(:reBiGG:Recon3D)-[:PARTICIPATES]-(m:metBiGG:Recon3D) 
     RETURN p") %>% 
  call_neo4j(con,type = "graph") 

## NODES
# only extract first label for each node 
G$nodes$label <- lapply( G$nodes$label, function(x) x[[1]][1])
# unnest 
G$nodes <- G$nodes %>%
  unnest_nodes(what = "label")
nodes <- data.frame(id=G$nodes$id,group=G$nodes$value)
head(G$nodes)  

## EDGES
relationships <- G$relationships %>%
  unnest_relationships() %>%
  dplyr::select(from = startNode, to = endNode)
head(relationships)

## GRAPH
visNetwork::visNetwork(nodes, relationships)
