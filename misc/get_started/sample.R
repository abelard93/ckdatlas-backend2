#load packages
library(RNeo4j)

#connect to database
db_addr <- "http://dzdcon1:6464/db/data/"
db = startGraph(db_addr)

# important functions: cypher() and cypherToList()
# cypher() - used when returning a dataframe/ resulty in table form
# cypherToList() - used when returning whole node "objects", provides list of nodes with all properties



###################
# Example queries #
###################

#---------------------------------------
# [1] get all measured metabolites from:
#---------------------------------------

# a) specific file
#------------------

query="
      MATCH (m:metMetabolon)-[:FROM]-(s:source)
      WHERE s.sid='Karpas-20422-20met_alt'
      RETURN m.sid as metabolites"

cypher(db,query)

#alternative

query="
      MATCH (m:metMetabolon)-[:FROM]-(:source {sid:'Karpas-20422-20met_alt'})
      RETURN m.sid as metabolite, m.BIOCHEMICAL as name"

cypher(db,query)


# b) all files starting with 'HELM'
#---------------------------------

#need to use DISTINCT to return unique metabolite ids

query="
      MATCH (m:metMetabolon)-[:FROM]-(s:source)
      WHERE s.sid STARTS WITH 'HELM'
      RETURN DISTINCT m.sid as metabolites "

cypher(db,query)


# c) liver samples
#-----------------

query="
      MATCH (m:metMetabolon)-[:FROM]-(s:source)
      WHERE s.sample_type='liver'
      RETURN DISTINCT m.sid as metabolites"

cypher(db,query)


# d) liver samples with list of measured metabolites per file/study
#------------------------------------------------------------------

query="
      MATCH (m:metMetabolon)-[:FROM]-(s:source {sample_type:'liver'})
      RETURN s.sid as source, COLLECT(m.sid) as metabolites"

res=cypher(db,query)


#------------------------------------------------------------------------------------------------
# [2] provide data table containing all metabolite-gene pairings that are "connected" by reaction
#------------------------------------------------------------------------------------------------

# a) for mouse, also including reactions
#---------------------------------------

# returning table with metabolites participating/reaction/genes catalyzing reaction

query="
      MATCH (g:gene)-[:CATALYZES]-(r:reBiGG)-[:PARTICIPATES]-(n:metBiGG)-[:MAP*]-(m:metMetabolon)
      WHERE (n)-[:FROM]-(:source {sid:'iMM1415'}) 
      return collect(distinct m.BIOCHEMICAL) as Metabolites,  r.sid as reaction, collect(distinct g.Symbol) as Genes"

#pass list of metabolites to query
cypher(db,query)


# b) for specific list of metabolites, only reporting metabolite/genes
#---------------------------------------------------------------------

#metabolites that should be inspected
metabos <- c("35631", "35305", "33961", "34416", "35628", "45968", "36602", "33955", "35186", "34214", "34419", "32635", "20675", "39270", "15681", "3155", "56","18289", "1419", "1494", "32484")

#list of metabolites resides in 'properties'
#by unwinding we go through every entry/metabolite

query="
      UNWIND {properties} as metabo
      MATCH (g:gene)-[:CATALYZES]-(:reBiGG)-[:PARTICIPATES]-(n:metBiGG)-[:MAP*]-(m:metMetabolon)
      WHERE (n)-[:FROM]-(:source {sid:'iMM1415'}) and m.sid = metabo 
      return m.sid as Metabolite, collect(distinct g.sid) as Genes"

#pass list of metabolites to query
cypher(db,query,properties = metabos)

#---------------------------------------------------------------------------------------------
# [3] visualize (GUI - http://ibisdb02.scidom.de:7475/browser/) or extract all nodes of a path 
#---------------------------------------------------------------------------------------------

query="
       MATCH p = (g:gene {sid:'29953'})-[:CATALYZES]-(:reBiGG)-[:PARTICIPATES]-(n:metBiGG)-[:MAP*]-(m:metMetabolon {sid:'1494'}) 
       return p"

#returns object that contains path object in [[1]]$p and length of path in [[1]]$p$length
path = cypherToList(db,query)

#extract nodes in path object
nodes(path[[1]]$p)
