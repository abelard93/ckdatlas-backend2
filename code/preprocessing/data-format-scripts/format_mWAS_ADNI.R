library(tidyverse)
zap()

########
## mWAS
########
zap()

# READ DATA
mapping <- read.csv("/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/MWAS outcomes & GGMs/mWAS_phenotype_mapping.csv") %>%
  dplyr::transmute(sid = name,
                   id = MWAS_trait_name,
                   metaTrait = phenotype_class)
ADNI_bile <- data.table::fread(
  "/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/MWAS outcomes & GGMs/ADNI-1GO2-bile_acids_trait_associations.csv"
)
ADNI_p180 <- data.table::fread(
  "/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/MWAS outcomes & GGMs/ADNI-1GO2-p180_trait_associations.csv"
)

# LOOP THROUGH TRAITS (1)
traits_bile <- ADNI_bile %>% .$trait %>% unique
for (i in traits_bile) {
  tmp <- ADNI_bile %>% 
    filter(trait == i) %>%
    transmute(
      met,
      trait = mapping %>% dplyr::filter(id == i) %>% .$sid, # convert to trait name (use mapping table)
      beta = cross.b, # only crossectional for now
      SE = cross.se,
      pvalue = cross.p,
      qvalue = cross.q
    )
  name <- paste0("ADNI_bile_", str_replace_all(str_replace_all(str_replace_all(i, "[^A-Za-z0-9]|' '", "_"), "__", "_"), "_$", "")) # fix name
  
  write.table(
    tmp,
    file = paste0("/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/QTLs/mWAS/ADNI_bile/",name,".csv"),
    sep = ",",
    col.names = T,
    row.names = F
  )
  # FILE INFO
  data = c(
    "file_name", paste0(name, ".csv"),
    "sid", name,
    "organism", "Homo Sapiens",
    "node_type", "mWAS",
    "metaboliteID", "bileAcidKit",
    "cohort", "ADNI"
  )
  # MAKE DF
  info <- data.frame(matrix(data = data, ncol = 2, byrow = T))
  # WRITE
  write.table(
    info,
    file = file.path(
      "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/QTLs/mWAS/ADNI_bile",
      paste0(name, "_info.csv")
    ),
    sep = ",",
    col.names = F,
    row.names = F
  )
  
}

# LOOP THROUGH TRAITS (2)
traits_p180 <- ADNI_p180 %>% .$trait %>% unique
for (i in traits_p180) {
  tmp <- ADNI_p180 %>% 
    filter(trait == i) %>%
    transmute(
      met,
      trait = mapping %>% dplyr::filter(id == i) %>% .$sid, # convert to trait name (use mapping table)
      beta = cross.b, # only crossectional for now
      SE = cross.se,
      pvalue = cross.p
    )
  name <- paste0("ADNI_p180_", str_replace_all(str_replace_all(str_replace_all(i, "[^A-Za-z0-9]|' '", "_"), "__", "_"), "_$", ""))
  # WRITE
  write.table(
    tmp,
    file = paste0("/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/QTLs/mWAS/ADNI_p180/",name,".csv"),
    sep = ",",
    col.names = T,
    row.names = F
  )
  # FILE INFO 
  data = c(
    "file_name", paste0(name, ".csv"),
    "sid", name,
    "organism", "Homo Sapiens",
    "node_type", "mWAS",
    "metaboliteID", "p180",
    "cohort", "ADNI"
  )
  # MAKE DF
  info <- data.frame(matrix(data = data,ncol = 2,byrow = T))
  # WRITE
  write.table(
    info,
    file = file.path(
      "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/QTLs/mWAS/ADNI_p180",
      paste0(name, "_info.csv")
    ),
    sep = ",",
    col.names = F,
    row.names = F
  )
}
# WRITE TRAITS
write.table(
  mapping,
  file = paste0(
    "/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/trait/",
    "ADNI_bile_p180_traits.csv"
  ),
  sep = ",",
  col.names = T,
  row.names = F
)

#####################