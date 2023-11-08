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

## fix
biocrates_meta$p150[startsWith(biocrates_meta$name,"PC ae C")] <- biocrates_meta$p180[startsWith(biocrates_meta$name,"PC ae C")]
biocrates_meta$name[biocrates_meta$name=="Glutarylcarnitine* (Hydroxyhexanoylcarnitine)"] <- "Glutarylcarnitine (Hydroxyhexanoylcarnitine)"

# BIOCRATES mapping 
# READ AND FORMAT
readxl::excel_sheets("~/p180_annotation.xlsx")
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

p180_meta$other_name[startsWith(p180_meta$p180,"lysoPC")] <-paste(p180_meta$name[startsWith(p180_meta$p180,"lysoPC")],p180_meta$other_name[startsWith(p180_meta$p180,"lysoPC")],sep="|")
p180_meta$name[startsWith(p180_meta$p180,"lysoPC")] <- p180_meta$p180[startsWith(p180_meta$p180,"lysoPC")] 

# JOIN META AND MAPPING INFO 
biocrates_meta <- biocrates_meta %>% full_join(., p180_meta) %>%
  ## GROUP BY EXACT SAME NAME
  group_by(name) %>% # one entry for each metabolite
  summarise_all(.funs = ~ paste(unique(na.omit(.)), collapse = '|')) %>%
  ungroup() %>% 
  mutate(metabolon_name = str_replace(metabolon_name,",","|")) %>%
  mutate(metabolon_name = stringr::str_replace_all(metabolon_name, "\\(([a-z\\s]*)\\:excluded\\)", "\\1")) %>%
  tidyr::unite("metName",name,other_name,metabolon_name, remove = F ,sep = "|") %>%
  select(-other_name,-metabolon_name) %>%
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
  dplyr::rename(sid=value, metMetabolon=metabolon_id) %>% 
  mutate(
    metName = map2_chr(metName, sid, function(x,y){
      
      tmp <- c(strsplit(x, "\\|") %>% unlist,strsplit(y, "\\|") %>% unlist) %>% unique
      tmp_processed <- tmp %>% unlist %>% tolower %>% stringr::str_replace_all("[^a-z0-9]","") %>% unique
      c(tmp,tmp_processed) %>% unlist %>% paste0(.,collapse = "|")

    }),
    metMetabolon = stringr::str_replace_all(metMetabolon, "\\(([0-9\\s]*)\\:ex[c]*luded\\)", "\\1")
  )

biocrates_meta$join_id = c(1:nrow(biocrates_meta))

biocrates_meta <- biocrates_meta %>% select(sid, everything())

biocrates_meta <- biocrates_meta %>% mutate(sid = strsplit(sid, "\\|"), metName = str_replace_all(metName,"\\|\\|","|")) %>%
  unnest(sid) 



## PATRICKS META FILE ############


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
    }),
    sid = Metabolite
  ) %>% select(-Metabolite)

anno[anno == ""] <- NA
anno[anno == "NA"] <- NA

names(anno) <- c("metName","metHMDB","metKEGG","metFormula","metIUPAC","metCAS","metSmiles","metClassification","metINCHI","mapKEGG","SMPDB","sid")
anno <- anno %>% select(sid, everything()) %>% mutate(platform="p180",  metName = str_replace_all(metName,"\\|\\|","|"))

##### CONCATINATE BOTH FILES 

meta <- rbind.fill(biocrates_meta,anno) %>% group_by(sid) %>%
  summarise_all(.funs = ~ paste(unique(na.omit(.)), collapse = '|')) %>% # group together metabolites with same id
  ungroup() %>%
  mutate_all(na_if, "") 

meta$join_id[is.na(meta$join_id)] <- 999

meta <- meta %>% mutate(join_id = strsplit(join_id, "\\|")) %>%
  unnest(join_id) %>% group_by(join_id) %>%
  summarise_all(.funs = ~ paste(unique(na.omit(.)), collapse = '|')) %>% # group together metabolites with same id
  ungroup() %>%
  mutate_all(na_if, "")  %>% select(-join_id)


### FORMATING ISSUES 

meta <- meta %>% mutate(metName = str_replace(metName,"\\|\\|","|"))

### NAMING ISSUES 
meta$name[meta$sid=="PC.ae.C32.1"] <- "PC ae C32:1"
meta$name[meta$sid=="PC.ae.C32.2"] <- "PC ae C32:2"
meta$name[is.na(meta$name)] <- meta$sid[is.na(meta$name)]

meta <- meta %>% mutate(name = map_chr(name, function(x) { str_split(x,"\\|")[[1]][1]}),
                        platform = map_chr(platform, function(x){ str_split(x,"\\|")[[1]] %>% unique %>% paste0(collapse = "|")}))

write.csv(meta,"metBiocrates/biocrates_meta.csv",row.names = F,quote = T,na = "")

