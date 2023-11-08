library(tidyverse)
zap()

####################
## traitQTLs
####################
zap()
mapping <- readxl::read_excel("~/vmetaneo4j/Phenotype_mapping.xlsx")

# PREPARE FILE INFO
ADNI <- c("cohort", "ADNI")
BEECHAM <-
  c("cohort",
    paste0(
      c(
        'ACT',
        'ADCs',
        'MAYO',
        'MBB',
        'NIA-LOAD',
        'OHSU',
        'ROSMAP',
        'TGEN2',
        'UM',
        'VU',
        'MSSM',
        'UP',
        'ADGC'
      ),
      collapse = "|"
    ),
    "publication",
    "beecham_2014_25188341")
DEMING <-
  c("cohort",
    paste0(
      c(
        'Knight ADRC',
        'ADNI',
        'BIOCARD',
        'HB',
        'MAYO',
        'SWEDEN',
        'UPENN',
        'UW'
      ),
      collapse = "|"
    ),
    "publication",
    "deming_2017_28247064")
IGAP <-
  c("cohort",
    paste0(
      c(
        'I-GAP',
        'ADGC',
        'CHARGE',
        'EADI',
        'GERAD',
        'ACT',
        'ADC1',
        'ADC2',
        'ADC3',
        'ADNI',
        'GSK',
        'LOAD',
        'MAYO',
        'MIRAGE',
        'OHSU',
        'ROSMAP',
        'TGEN2',
        'UMVUMSS',
        'UPITT',
        'WASHU',
        'AGES',
        'CHS.inc',
        'CHS.prev',
        'FHS.inc',
        'FHS.prev',
        'RS.inc',
        'RS.prev'
      ),
      collapse = "|"
    ),
    "publication",
    "lambert_2013_24162737")

# LOOP THROUGH STUDIES
for (sub in unique(mapping$subset)) {
  # READ MAPPING
  mapping <- readxl::read_excel("~/vmetaneo4j/Phenotype_mapping.xlsx") %>% 
    dplyr::filter(subset == sub)
  traits <- mapping %>% 
    dplyr::transmute(
      sid = name, 
      id = variable,
      metaTrait = `phenotype class`)
  folder <- strsplit(sub, "_")[[1]][1]
  # LOOP THROUGH TRAITS
  for (file in mapping[, "filename"][[1]]) {
    print(file)
    # READ FILE
    x = load(paste0("~/vmetaneo4j/data/", file))
    y = get(x)
    rm(x)
    names(y) = c("SNP",
                 "chromosome",
                 "bp",
                 "effect_allele",
                 "effect_direction",
                 "pvalue")
    y <- y %>% 
      dplyr::mutate(
        trait = mapping %>% dplyr::filter(filename == file) %>% .$name, # convert to trait name (use mapping table)
        effect_allele = toupper(effect_allele)
      )
    file <- tools::file_path_sans_ext(file) %>% str_replace_all(., "[.]", "_") %>% paste0(., ".csv")
    # WRITE
    write.table(
      y,
      file = file.path("/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/QTLs/traitQTL", folder, file),
      sep = ",",
      col.names = T,
      row.names = F
    )
    # COMPRESS
    R.utils::gzip(file.path("/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/QTLs/traitQTL", folder, file))
    rm(y)
    
    ## INFO
    data = c(
      "file_name", file,
      "sid", gsub(".csv", "", file),
      "organism", "Homo Sapiens",
      "node_type", "traitQTL",
      get(folder) # cohort info
    )
    # MAKE DF
    info <- data.frame(matrix(data = data, ncol = 2, byrow = T))
    # WRITE
    write.table(
      info,
      file = file.path("/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/QTLs/traitQTL", folder, gsub(pattern = ".csv", replacement = "_info.csv", x = file)),
      sep = ",",
      col.names = F,
      row.names = F
    )
  }
}

# WRITE TRAITS
write.table(
  traits,
  file = paste0("~/vmetaneo4j/neo4j/", "ADNI_traits.csv"),
  sep = ",",
  col.names = T,
  row.names = F
)

#####################

############################
## mQTL: metabo ID mapping
############################
# READ
zap()


id_mapping <- read.csv("/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/Population-based metabolite QTLs/id_map/shin_suhre.id.map.csv", stringsAsFactors = F) %>% 
  transmute(
    metabolite = metabolonDescription,
    COMP_ID = map_chr(metabolonID, function(x) str_remove(x, "^M0*")))

# SET WD
setwd("/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/Population-based metabolite QTLs/")

# LOOP THROUGH SPECIFIC FILE
for (file in list.files(".",full.names = T,pattern = "*.txt")){
  
  setwd("/Volumes/HMGU_maria/DukeBox/Box Sync/Molecular Atlas/Population-based metabolite QTLs/")
  old <- read.csv(file, stringsAsFactors = F,sep = "\t") %>% select(-pubmedid,-pubmedlink)
  names(old) <- c("rsID","metabolite","pvalue","source","fluid")
  
  if (basename(file) %in% c("shin_et_al.txt", "suhre_et_al_kora.txt", "suhre_et_al_twinsuk.txt")) {
    
    # ADD ID
    old <- old %>% left_join(., id_mapping)
  }
  # SET WD
  setwd("/Volumes/HMGU_maria/DukeBox/Box Sync/Neo4J Data Repository/neo4jserver/data/QTLs/mQTL/population_based/")
  # WRITE
  write.table(
    old,
    file = gsub(pattern = ".txt", replacement = ".csv", x = file),
    sep = ",",
    col.names = T,
    row.names = F
  )
  
}

#####################

############################
## mQTL: tvs to csv
############################
source(codes.makepath("neo4jserver/code/setup/toLoad/setupFunctions.R"))
for (file in getFiles("QTLs/mQTL")) {
  tsv <- read.delim(file)
  write.csv(
    tsv,
    sub(".txt", ".csv", file),
    row.names = F,
    quote = T,
    na = ""
  )
}

#####################