
#
# Paper: A Consensus Proteomic Analysis of Alzheimer’s Disease Brain 
#       and Cerebrospinal Fluid Reveals Early Changes in Energy Metabolism 
#       Associated with Microglia and Astrocyte Activation
# DOI: https://doi.org/10.1101/802959
#
# Maria Wörheide
# 10-03-2020
#

zap()
library(dplyr)
library(RNeo4j)
library(tidyverse)
db = startGraph("http://metabo-interactive:7474/db/data/", username = "neo4j", password = "sysdiab")

setwd("/Users/mwoerheide/Documents/PHD/projects/AD_atlas_proteomics/johnson_et_al/")



### -----------------------  data preprocessing  ----------------------- 

### ----------------------- ID mapping ----------------------------

proteins <- read.csv("data/3.cleanDat.csv",stringsAsFactors = F) %>%
  .$X %>% 
  # split name|uniprotID
  stringr::str_split_fixed(.,"\\|",2) %>% 
  as.data.frame() %>% 
  `colnames<-`(c("name","uniprot")) %>% mutate_all(as.character) %>%
  mutate(id = read.csv("data/3.cleanDat.csv",stringsAsFactors = F) %>%
           .$X)

id_mapping <- setNames(proteins$uniprot,proteins$id)
# split any info from uniprotID after "-"
# proteins[,c(2,3)] <- stringr::str_split_fixed(proteins$uniprot,"-",2)
# names(proteins)[3] <- "additional_info"
samples <- read.csv("data/3.cleanDat.csv") %>% dplyr::select(-1) %>% colnames

### ----------------------- meta df----------------------------
meta_df <- read.csv("data/0.Traits.csv") %>% 
  dplyr::select(X,Age,Sex,PMI, Group) %>% 
  # filter samples
  filter(X%in%samples) %>% 
  tibble::column_to_rownames(var = "X")

### ----------------------- protein df----------------------------
protein_df <- read.csv("data/2b.Minimally_regressed_Batch_and_Site-corrected_cleanDat.csv",stringsAsFactors = F) %>% 
  # filter proteins
  filter(X%in%proteins$id) %>% 
  # filter samples
  dplyr::select(one_of(c("X",samples))) %>% 
  # id mapping
  mutate(X=purrr::map_chr(X, function(x) id_mapping[x])) %>% tibble::column_to_rownames(var = "X")

### ----------------------- DEP ----------------------------
dep <- readxl::read_excel("data/media-3.xlsx",
                          sheet = 2,skip = 1) 
dep <- dep[,c(1,443:450)]
names(dep) <- c("protein","fvalue","pr(>f)","AsymAD vs. AD","Control vs. AD", "Control vs. AsymAD", "log2(FC)_AsymADAD","log2(FC)_ControlAD","log2(FC)_ControlAsymAD")

dep <- dep %>% separate(protein,into = c("name","uniprot"),sep = "\\|") %>% na.omit() %>% mutate(info="avg log2(group1) - avg log2(group2)")

### ----------------------- save ----------------------------
save(meta_df,protein_df,dep, file = "data/proteomics.Rdata")

### ----------------------- LOAD DATA ----------------------------------
setwd("/Users/mwoerheide/Documents/PHD/projects/AD_atlas_proteomics/johnson_et_al/")

load("data/proteomics.Rdata")

### ----------------------- FORMAT -------------------------------------

names(dep)[2] <- c("protein_id")

### -------------- LOOP THROUGH COMPARISONS ----------------------------

 
for(i in c(5,6,7)){
  
  df <- dep[,c(2,3,4,i,i+3,11)]
  names(df)[4:5] <- c("pvalue","log2(FC)")
  df <- df %>% transmute(protein_id, 
                         fvalue, 
                         pvalue_fdr=pvalue, 
                         `pr(>f)`,
                         `log2(FC)`, 
                         ## unchanged : not significant 
                         ## up : CTRL-CONDITION < 1 
                         direction = ifelse(pvalue_fdr>0.05,"unchanged",ifelse(`log2(FC)`>0,"down","up")),
                         trait=names(dep)[i],
                         info
                      ) %>%
    mutate(uniprot_id=str_replace(protein_id,"[-].*",""))
  
  ### ----------------------- WRITE -------------------------------------
  write.csv(df,paste0("data/dep_emory_",str_remove_all(names(dep)[i],pattern = "[.]|[ ]|vs"),".csv"),row.names = F,quote = T,na = "")
  
  ## FILE INFO
  info <- data.frame(matrix(ncol = 2,nrow = 11))
  info[1, ] <- c("sid",paste0("dep_emory_",str_remove_all(names(dep)[i],pattern = "[.]|[ ]|vs")))
  info[2, ] <- c("node_type","edgeDEP")
  info[3, ] <- c("file_name",paste0("dep_emory_",str_remove_all(names(dep)[i],pattern = "[.]|[ ]|vs"),".csv"))
  info[4, ] <- c("sample_type","brain|DLPFC")
  info[5, ] <- c("cohort","BLSA|ACT|MSSB|Banner")
  info[6, ] <- c("tax","hsa")
  info[7, ] <- c("taxID","9606")
  info[8, ] <- c("organism","Homo sapiens")
  info[9, ] <- c("publication","johnson_2020")
  info[10, ] <- c("significance","fdr<=0.05")
  info[11, ] <- c("status","disease")
  # WRITE
  write.table(
    info,
    file = paste0("data/dep_emory_",str_remove_all(names(dep)[i],pattern = "[.]|[ ]|vs"),"_info.csv"),
    sep = ",",
    col.names = F,
    row.names = F
  )
    
}
