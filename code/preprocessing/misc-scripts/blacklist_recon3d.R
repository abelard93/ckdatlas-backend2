## BLACKLIST FOR RECON3D
# using Regex List from Jonas

zap()

require(readxl)
require(tidyverse)
require(RNeo4j)
# Connect to db:
url = "http://localhost:7474/db/data/"
db = startGraph(url, username = "neo4j",password = "neo4j")

## read in excel sheets
metabo <- read_excel("~/Desktop/blacklist_wo_glucose.xlsx",sheet = "metabolites") %>% 
  select( field, regex) %>% 
  mutate(regex = map_chr(regex,str_replace_all,pattern="\\\\\\[",replacement="_")) %>%
  mutate(field = map_chr(field,str_replace_all,pattern="id",replacement="sid"))
react <- read_excel("~/Desktop/blacklist_wo_glucose.xlsx", sheet = "reactions") %>% select( field, regex) %>%
  mutate(field = map_chr(field,str_replace_all,pattern="id",replacement="sid"))


query <- paste0("
                    MATCH (m:metBiGG:recon3D)
                    RETURN distinct m.sid as sid
                    ")

m <- cypher(db,query)

query <- paste0("
                    MATCH (r:reBiGG:recon3D)
                    OPTIONAL MATCH (m:metBiGG:recon3D)-[:PARTICIPATES]-(r)
                    RETURN distinct r.sid as sid, r.subsystem as subsystem, r.name as name, m.sid as metab
                    ")

r <- cypher(db,query)


get_ids<- function(field, regex, data){
  
   qf <- rlang::sym(field)
   ids <- data %>%
     filter(str_detect(!!qf, regex(regex, ignore_case = T))) %>%
     .$sid
  
   ids
}


bl_metab <- metabo %>%
  mutate(ids = map2(field, regex, get_ids, data = m)) %>%
  unnest(ids) 

bl_react <- react %>%
  mutate(ids = map2(field, regex, get_ids, data = r)) %>%
  unnest(ids) 

blacklist.recon3D <- as.data.frame(unique(c(bl_metab$ids, bl_react$ids))) 
names(blacklist.recon3D) <- "sid"

write.csv2(blacklist.recon3D,"~/Documents/ICB/neo4j/data/blacklist/recon3D_maria.csv",row.names = F)
