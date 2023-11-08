library(tidyverse)
zap()


######################################
## metabolite mapping 
## matze: biocrates on 25-07-2019
######################################

# BIOCRATES meta files
# READ AND FORMAT
readxl::excel_sheets("~/Cross-platform_metabolite_mapping_v1.xlsx")
biocrates_meta <- readxl::read_excel(
  "~/Cross-platform_metabolite_mapping_v1.xlsx", 
  sheet=1,
  col_names = T,
  na = "NA") %>% 
  transmute(
    name = `Curated metabolite name (from Biocrates)`,
    p180 = `ADNI Biocrates P180`,
    bileAcidKit = `ADNI Biocrates bile acid kit`,
    p150 = `Draisma Biocrates P150`,
    metabolon_id = strsplit(`Shin Metabolon Id`,","),
    metClassification = ifelse(is.na(`Biocrates classification`),`Bile acid classification`,`Biocrates classification`),
    metabolon_name = `Shin Metabolon name`
  ) %>% 
  unnest(metabolon_id) %>% 
  mutate(metabolon_id = str_remove(str_remove_all(metabolon_id, "M"), "^0+"))

# BIOCRATES mapping 
# READ AND FORMAT
readxl::excel_sheets("/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/metabolite_mapping/p180_annotation.xlsx")
p180_meta <- readxl::read_excel(
  "~/p180_annotation.xlsx", 
  sheet=1,
  col_names = T,
  na = "NA") %>% 
  mutate(
    CHEBI = ifelse(is.na(CHEBI),NA,paste0("CHEBI:",CHEBI))
  ) %>% 
  group_by(SHORT_NAME) %>% 
  summarise_all(.funs = ~paste(unique(na.omit(.)), collapse = '|')) %>% 
  ungroup() %>% 
  transmute(
    p180 = SHORT_NAME,
    name = NAME,
    other_name = MOLECULE,
    metKEGG = KEGG,
    metCAS = CAS,
    metChEBI = CHEBI,
    metHMDB = HMDB,
    metWeight = MOLECULEWEIGHT,
    metFormula = MOLECULARFORMULA,
    metMESH = MESH
  )
# EMPTY STRING TO NA
p180_meta[p180_meta == ""] <- NA

# JOIN META AND MAPPING INFO 
biocrates_meta <- biocrates_meta %>% full_join(., p180_meta) %>%
  ## GROUP BY EXACT SAME NAME
  group_by(name) %>% # one entry for each metabolite
  summarise_all(.funs = ~ paste(unique(na.omit(.)), collapse = '|')) %>%
  ungroup() %>%
  ## COLLECT NAME * SYNONYM
  mutate(
    #name = tidyr::unite("name",name,other_name,remove = TRUE,sep = "|")
    name = map2_chr(name, other_name, function (x,y){ paste(x,y,sep="|")
  }) 
   
  ) %>%
  dplyr::select(-other_name) %>%
  mutate_all(na_if, "") %>% # empty string to na
  ## WIDE --> LONG
  gather(key = "platform", value = "value", p180, bileAcidKit, p150) %>%
  filter(!is.na(value)) %>% # delete metabolites not measured on platform
  ## GROUP BY SAME ID
  group_by(value) %>%
  summarise_all(.funs = ~ paste(unique(na.omit(.)), collapse = '|')) %>% # group together metabolites with same id
  ungroup() %>%
  mutate_all(na_if, "") %>% # empty string to na
  ## UNNEST IDS AND JOIN BACK TOGETHER 
  # problem with multiple names on one platform eg. C12.1|C12:1 (p180) and C12.1 (p150)
  mutate(value = strsplit(value, "\\|")) %>%
  unnest(value) %>%
  group_by_at(vars(-value, -platform)) %>% 
  summarise_all(.funs = ~ paste(unique(na.omit(.)), collapse = '|')) %>% ungroup %>% 
  dplyr::rename(metName=name, sid=value, metMetabolon=metabolon_id) %>% 
  mutate(
    metName = map2_chr(metName, sid, function(x,y){
      
      tmp <- c(strsplit(x, "\\|") %>% unlist,strsplit(y, "\\|") %>% unlist) %>% unique
      tmp_processed <- tmp %>% unlist %>% tolower %>% stringr::str_replace_all("[^a-z0-9]","") %>% unique
      c(tmp,tmp_processed) %>% unlist %>% paste0(.,collapse = "|")
  
    #  paste(x, y ,sep = "|")
    })#,
  #  sid = map_chr(sid, ~ stringi::stri_split_fixed(., pattern = "|", n = 2)[[1]][1])
    )

biocrates_meta <- biocrates_meta %>% select(sid, everything())

write.csv(biocrates_meta,"/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/metBiocrates/biocrates_meta.csv",row.names = F,quote = T,na = "")



## PATRICKS META FILE ############


library(RNeo4j)

library(tidyverse)

biocrates_meta <- read.csv("metBiocrates/biocrates_meta.csv")
## with mapping

anno <- readxl::read_excel("~/biocrates_anno.xlsx",trim_ws = T) %>% 
 # mutate(`biocrates_synonymes`= str_replace_all(`biocrates_synonymes`, ", ","|")) %>% 
 # unite(.,metName,`biocrates_commonly_used_names`,`biocrates_synonymes`,sep = "\n",remove = T,na.rm=T) %>% 
  unite(.,metName,`biocrates_commonly_used_names`,sep = "\n",remove = T,na.rm=T) %>% 
  select(-`biocrates_synonymes`) %>%
  unite(.,metClassification,biocrates_super_class,biocrates_super_path,biocrates_sub_path,class,sep = "\n",remove = T,na.rm=T) %>% 
  mutate_all(list(~map_chr(.,function(x) str_split(x, "[\r\n][\r\n]|[\r\n]") %>% unlist %>% unique %>% unlist %>% paste0(.,collapse = "|") ))) %>% 
  mutate(kegg = map_chr(pathways, function(x)
    regmatches(x, gregexpr('map[0-9]{5}', x)) %>% unlist %>% unique %>% unlist %>% paste0(collapse = "|")
    
  ),
  smpdb =  map_chr(pathways, function(x)
    regmatches(x, gregexpr('SMP[0-9]{5}', x)) %>% unlist %>% unique %>% unlist %>% paste0(collapse = "|")
  )) %>% 
  mutate_all(list(~map_chr(.,function(x) str_replace(x,"\\|$","")))) %>% select(-pathways,-descriptor,-levels) %>% 
  mutate(
    metName = map2_chr(metName, biocrates_name, function(x,y){
      tmp <- c(strsplit(x, "\\|") %>% unlist, strsplit(y, "\\|") %>% unlist) %>% unique
      tmp_processed <- tmp %>% unlist %>% tolower %>% stringr::str_replace_all("[^a-z0-9]","") %>% unique
      c(tmp,tmp_processed) %>% unlist %>% paste0(.,collapse = "|")
      
      #  paste(x, y ,sep = "|")
    })#,
    #  sid = map_chr(sid, ~ stringi::stri_split_fixed(., pattern = "|", n = 2)[[1]][1])
  )

anno[anno == ""] <- NA
anno[anno == "NA"] <- NA

#### get names
db_addr <- "http://metabo-interactive:7474/db/data/"
## CONNECT TO NEO4J
db = startGraph(db_addr, username = "neo4j", password = "sysdiab")


metabo_name <- "UNWIND {properties} AS prop with prop return [(m:metName {sid:prop})-[:HAS_NAME]-(:metMetabolon)|m.sid][0] as name" %>%
  cypher(db,.,properties=anno$metabolon_name) 

metabo_name_matched <- "UNWIND {properties} AS prop with prop return [(m:metName {sid:prop})-[:HAS_NAME]-(:metMetabolon)|m.sid][0] as name" %>%
  cypher(db,.,properties=anno$metabolon_name_matched)


anno$name <- purrr::map2_chr(metabo_name$name,metabo_name_matched$name, function(x,y){
  ifelse(!is.na(y),y,x)
}) %>% purrr::map2_chr(.,anno$metabolon_name_matched,function(x,y){
  ifelse(!is.na(x),x,y)
})

anno <- anno %>% dplyr::select(-metabolon_name_matched,-metabolon_name)
names(anno) <- c("id", "metName","metHMDB","metKEGG","metFormula","metIUPAC","metCAS","metSmiles","metClassification","metINCHI","metSuperPathway","metSubPathway", "metMetabolon", "mapKEGG","SMPDB","name")
# add corresponding sid from meta 

anno$sid <- sapply(anno$id,function (y) {
  index <- sapply(biocrates_meta$sid,function (x) {
    y%in%(str_split(x,"\\|") %>% unlist)
  }) %>% unlist %>% which 
  ifelse(is_empty(biocrates_meta$sid[index]),y,biocrates_meta$sid[index])
} 
) %>% unlist %>% as.vector

anno <- anno %>% select(sid, name, everything()) %>% mutate(metMetabolon=str_remove(metMetabolon,"M|^0*"))

write.csv(anno,"/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/metBiocrates/biocrates_annotation_meta.csv",row.names = F,quote = T,na = "")




## without mapping 


anno <- readxl::read_excel("~/biocrates_anno.xlsx",trim_ws = T) %>%
 # mutate(`Attached Synonyms`= str_replace_all(`Attached Synonyms`, ", ","|")) %>% 
 # unite(.,metName,`Commonly used names`,`Attached Synonyms`,sep = "\r\n",remove = T,na.rm=T) %>% 
  unite(.,metName,`Commonly used names`,sep = "\r\n",remove = T,na.rm=T) %>% 
  select(-`Attached Synonyms`) %>%
  unite(.,metClassification,Super_Class,Superpath,Subpath,Class,sep = "\n",remove = T,na.rm=T) %>% 
  mutate_all(list(~map_chr(.,function(x) str_split(x, "[\r\n][\r\n]|[\r\n]") %>% unlist %>% unique %>% unlist %>% paste0(.,collapse = "|") ))) %>% 
  mutate(kegg = map_chr(Pathways, function(x)
    regmatches(x, gregexpr('map[0-9]{5}', x)) %>% unlist %>% unique %>% unlist %>% paste0(collapse = "|")
    
  ),
  smpdb =  map_chr(Pathways, function(x)
    regmatches(x, gregexpr('SMP[0-9]{5}', x)) %>% unlist %>% unique %>% unlist %>% paste0(collapse = "|")
  )) %>% 
  mutate_all(list(~map_chr(.,function(x) str_replace(x,"\\|$","")))) %>% select(-Pathways,-Descriptor,-Levels) %>% 
  mutate(
    metName = map2_chr(metName, Metabolite, function(x,y){
      tmp <- c(strsplit(x, "\\|") %>% unlist, strsplit(y, "\\|") %>% unlist) %>% unique
      tmp_processed <- tmp %>% unlist %>% tolower %>% stringr::str_replace_all("[^a-z0-9]","") %>% unique
      c(tmp,tmp_processed) %>% unlist %>% paste0(.,collapse = "|")
      
      #  paste(x, y ,sep = "|")
    })#,
    #  sid = map_chr(sid, ~ stringi::stri_split_fixed(., pattern = "|", n = 2)[[1]][1])
  )

anno[anno == ""] <- NA
anno[anno == "NA"] <- NA

names(anno) <- c("id", "metName","metHMDB","metKEGG","metFormula","metIUPAC","metCAS","metSmiles","metClassification","metINCHI","mapKEGG","SMPDB")

# add corresponding sid from meta 

anno$sid <- sapply(anno$id,function (y) {
      index <- sapply(biocrates_meta$sid,function (x) {
        y%in%(str_split(x,"\\|") %>% unlist)
        }) %>% unlist %>% which 
      ifelse(is_empty(biocrates_meta$sid[index]),y,biocrates_meta$sid[index])
      } 
  ) %>% unlist %>% as.vector

anno <- anno %>% select(sid, everything())

write.csv(anno,"metBiocrates/biocrates_annotation_meta.csv",row.names = F,quote = T,na = "")


