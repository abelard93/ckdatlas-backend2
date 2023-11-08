library(tidyverse)
zap()


######################################
## metabolite mapping 
## parviz: metabolon on 19-07-2019
######################################

# METABOLON meta file
# READ
readxl::excel_sheets("/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/metabolite_mapping/metabolon_lookup_table_20190508.xlsx")
metabolon_meta <- readxl::read_excel(
  "/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/metabolite_mapping/metabolon_lookup_table_20190508.xlsx", 
  sheet=1,
  col_names = T)
# FORMAT
metabolon_meta <- metabolon_meta %>% 
  mutate(
    metabolonID = str_remove(str_remove(metabolonID, "M"), "^0+"), # remove leading 0 and M
    externalID_type = purrr::map_chr(externalID_type, function (x) # adapt external identifier name
      paste0("met", ifelse(x == "Recon", "BiGG", ifelse(x == "PUBCHEM", "PubChem", x)))
    )
  ) %>% 
  dplyr::filter(!(externalID_type == "metHMDB" & !grepl("^HMDB", externalID))) %>% #delete HMDB ids in wrong format
  mutate(row_id = row_number()) %>%  # add row id (needed for spread)
  spread(externalID_type, externalID) %>% # long -> wide
  dplyr::select(-row_id) %>%
  group_by_at(vars(metabolonID)) %>% # one entry for every comp_id
  summarise_all(.funs = ~ paste(unique(na.omit(.)), collapse = '|')) %>% # paste multiple entries together
  ungroup() %>% 
  mutate(metBiGG = map_chr(metBiGG,  ~ unlist(strsplit(.x, "-"))[1])) # fix format of BiGG id
# EMPTY STRING TO NA
metabolon_meta[metabolon_meta == ""] <- NA

# WRITE
write.csv(
  metabolon_meta,
  "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/metMetabolon/metMetabolon_meta.csv",
  row.names = F,
  quote = T,
  na = ""
)