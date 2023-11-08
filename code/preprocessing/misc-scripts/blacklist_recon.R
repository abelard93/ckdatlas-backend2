zap()
library(RNeo4j)
library(ggplot2)

## CONNECT TO DB 
url = "http://dzdcon1:6464/db/data/"
db = startGraph(url,username = "neo4j",password = "sysdiab")

## GET BIGG METABOLITES, KEGG_ID AND # OF REACTIONS
query=("MATCH (n:metBiGG)-[r:PARTICIPATES]-(:reBiGG)  RETURN n.id
       as metabolite,n.metNames as name, n.metKeggID as kegg_id, count(r) as reactions
       order by reactions DESC")

blacklist <- cypher(db,query)

## PLOT TO SEE WHICH METABO PARTICIPATES IN THE MOST REACTIONS
ggplot(blacklist[1:50,], aes(x=reorder(metabolite,order(blacklist[1:50,]$reactions,decreasing = T)), y=reactions,group=1)) +
  geom_path() + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("top 50 BiGG metabolites")

## READ IN BLACKLIST FROM PARVIZ
blacklist_parviz <- read.csv2("~/Documents/ICB/sharedcodes/packages/neo4jqueries/test_data/blackList_parviz_180522.csv")
## ONLY KEEP SUPPOSED COFACTORS
blacklist_parviz <- blacklist_parviz[blacklist_parviz$Modified=="TRUE",]

#sum(blacklist_parviz$kegg_id%in%blacklist$kegg_id)

## MERGE ON KEGG_ID
blacklist <- merge(blacklist,blacklist_parviz, all=T)
## SORT
blacklist <- blacklist[order(blacklist$reactions,decreasing = T),]
rownames(blacklist) <- NULL
## SET TOP 45 MOST FRQUENT BIGG METABO TO TRUE
blacklist[1:45,]$Modified <- TRUE
## MANUALL INSPECTION OF TOP 45
blacklist$Modified[blacklist$metabolite%in%c("crn","gly","gthrd")] <- FALSE
blacklist$Modified[blacklist$metabolite%in%c("fe2","fe3")] <- TRUE

## FORMAT
blacklist <-blacklist[blacklist$Modified=="TRUE" & !is.na(blacklist$Modified),]
rownames(blacklist) <- NULL


query=("UNWIND {ids} as id
        MATCH (n:metBiGG {id:id}) RETURN n.sid as sid,n.id as metabolite")

res <- cypher(db,query,ids=na.omit(blacklist$metabolite))

blacklist <- merge(res,blacklist,all=T)
## WRITE
write.csv2(blacklist,"~/Documents/ICB/sharedcodes/packages/neo4jqueries/test_data/blacklist.csv",row.names = F)
#recon3d <- load("Desktop/Recon3D_301_SBML_FixIDs/recon_graph.RData")



##########################
## REACTIONS
##########################

## GET BIGG REACTIONS AND # OF PARTICIPATING METABO
query=("MATCH (n:metBiGG)-[:PARTICIPATES]-(r:reBiGG)  RETURN r.sid
as reaction, r.rxnNames as names, count(n) as metabo
order by metabo DESC")

blacklist <- cypher(db,query)

## PLOT TO SEE WHICH METABO PARTICIPATES IN THE MOST REACTIONS
ggplot(blacklist[1:50,], aes(x=reorder(reaction,order(blacklist[1:50,]$metabo,decreasing = T)), y=metabo,group=1)) +
  geom_path() + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("top 50 BiGG reactions")

## CHOOSE TOP 49 REACTIONS: BIOMASS AND PROTEIN ASSEMBLY/DEGRADATION
blacklist <- blacklist[1:49,]

## WRITE
write.csv2(blacklist,"~/Documents/ICB/sharedcodes/packages/neo4jqueries/test_data/blacklist_reactions.csv",row.names = F)
#recon3d <- load("Desktop/Recon3D_301_SBML_FixIDs/recon_graph.RData")
