library(tidyverse)
zap()

########
## mGWAS
########


# READ DATA
mapping <- read_excel("vmetaneo4j/adni_mgwas/Global_Trait_Map.xlsx") %>%
  filter(subset=="ADNI_MET") %>% 
  dplyr::transmute(sid = name,
                   id = variable,
                   source = filename)

# LOOP THROUGH FILES (1)

for(i in 1:nrow(mapping)){
  
  load(file=paste0("vmetaneo4j/adni_mgwas/",mapping[i,"source"]))
  df <- get(mapping[[i,"id"]]) %>% 
    transmute(rsID = SNP,
           chromosome = CHR,
           bp = BP,
           effect_allele = EA,
           effect_direction = EFF.DIR,
           pvalue = PVALUE,
           metabolite = tolower(str_replace_all(mapping[[i,"sid"]],"[^0-9A-Za-z]",""))
    )
  
  rm(list = eval(mapping[[i, "id"]]))
  
  name <- str_replace_all(str_replace_all(tolower(str_replace_all(mapping[[i,"id"]],"\\.","_")),"__","_"),"_$","")
  
  write.table(
    df,
    file = paste0("vmetaneo4j/mgwas_neo4j/",name,".csv"),
    sep = ",",
    col.names = T,
    row.names = F
  )
  
  
  # FILE INFO
  data = c(
    "file_name", paste0(name, ".csv"),
    "sid", name,
    "organism", "Homo Sapiens",
    "node_type", "mQTL",
    "metaboliteID", "BIOCHEMICAL",
    "platform", "biocrates",
    "cohort", "ADNI"
  )
  # MAKE DF
  info <- data.frame(matrix(data = data, ncol = 2, byrow = T))
  # WRITE
  write.table(
    info,
    file = paste0("vmetaneo4j/mgwas_neo4j/",name,"_info.csv"),
    sep = ",",
    col.names = F,
    row.names = F
  )
}


