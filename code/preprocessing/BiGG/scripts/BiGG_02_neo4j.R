####
# Script from Jonas Zierer to parse recon3D
# adapted by Maria Woerheide for integration into Neo4j
# last modified: 2019-04-16
####
zap()
require(tidyverse)
require(xml2)
require(stringdist)
## GRAPHS
require(igraph)
library(zoo)
## PARAMS
args = commandArgs(trailingOnly = T)
file.rdata          <- args[1]
model.suffix        <- args[2]
file.ids.gene       <- args[3]
file.ids.reaction   <- args[4]
file.ids.metabo     <- args[5]

##
# file.recon.sbml     <- "~/Desktop/recon3D/Recon3D_301.xml"
# file.hgnc           <- "~/Documents/PHD/neo4j/neo4jserver/data/recon3d/genes_hgnc.rds"
# file.out            <- "~/Desktop/recon3D/recon3d.Rdata"
# file.ids.gene       <- "~/Documents/PHD/neo4j/neo4jserver/code/BiGG/data/20190411_recon_genes.tsv"
# file.ids.reaction   <- "~/Documents/PHD/neo4j/neo4jserver/code/BiGG/data/20190411_recon_reactions.tsv"
# file.ids.metabo     <- "~/Documents/PHD/neo4j/neo4jserver/code/BiGG/data/20190411_metabolites.tsv"
# model.suffix        <- "recon3D"

load(file.rdata)

## READ FILE WITH IDS
ext.ids <- all(!is.na(c(file.ids.gene,file.ids.metabo,file.ids.reaction)))

if(ext.ids){
  ids.metabo <- read_tsv(file.ids.metabo)
  ids.gene <- read_tsv(file.ids.gene)
  ids.reaction <- read_tsv(file.ids.reaction)  
}

################################################################################
## FORMAT for neo4j integration
################################################################################

#######################
## HELPER FUNCTIONS
#######################

## OLD HMDB TO NEW HMDB
old.to.new <- function(hmdb) {
  if (!is.na(hmdb) & str_length(hmdb) == 9) {
    sub("HMDB", "HMDB00", hmdb)
  } else{
    hmdb
  }
}

## HANDLE DUPLICATE COLUMNS AFTER JOIN
process_double <- function(df){
  # which column names are double (id not new)
  double <- grep("\\.x", names(df), value = T) %>% sub("\\.x", "", .) 
  for(i in double){
    # colnames
    x <- paste0(i,".x")
    y <- paste0(i,".y")
    ## if HMDB transform to new id 
    if(i=="hmdb"){
      df[,x] <- with(df,sapply(get(x),old.to.new))
      df[,y] <- with(df,sapply(get(y),old.to.new))
    }
    
    # new column named after property
    # messy but works
    df[,i] <- apply(df,1, function(z){
      
      if(is.na(z[x])&is.na(z[y])){ # both na
        NA
      }else if (is.na(z[x])&!is.na(z[y])){ # one na other contains id
        z[y]
      }else if(!is.na(z[x])&is.na(z[y])){ # one na other contains id
        z[x]
      }else if(z[x]!=z[y]){ # both contain id but not same
        paste(z[x],"|",z[y])
      }else{ # both contain same
        z[x]
      }
    })
    # remove  duplicate cols
    df <- df %>% select(-x,-y) 
  }
  return(df)
}




## NOTE: external id information downloaded from https://vmh.uni.lu/#human/all

#######################
## GENES
#######################

## JOIN WITH EXTERNAL IDS
if(ext.ids){
  neo4j.genes <- left_join(data.genes,
                           ids.gene %>% mutate(gene_number = map_chr(gene_number, as.character)),
                           by = c("label" = "gene_number")) %>% 
    mutate_all(trimws) %>% 
    process_double() %>% 
    mutate(sid=label) %>% 
    mutate(ensembl_gene = replace(ensembl_gene, !startsWith(ensembl_gene,"ENSG"), NA)) %>% unique() ## remove free text concrning pseduogenes 
  
}else{
  neo4j.genes <- data.genes %>% 
    mutate_all(trimws) %>% 
    mutate(sid=label) %>% unique() ## remove free text concrning pseduogenes 
}

## MANUAL FIXES
# remove gene_id 0 => where does this come from?
neo4j.genes <- neo4j.genes %>% filter(label!=0) %>% filter(label!="Unknown")
# microbe geneID with peg is not an entrez id .. set to NA
neo4j.genes$entrez_id[grepl(x=neo4j.genes$entrez_id,pattern = "*.peg.*")] <- NA

if(ext.ids){
  # same entrez_id => same gene information
  names.na <- neo4j.genes %>% filter(is.na(ensembl_gene)) %>% .$entrez_id
  for(i in names.na){
    neo4j.genes[neo4j.genes$entrez_id==i,c(3:ncol(neo4j.genes))] <- lapply(neo4j.genes[neo4j.genes$entrez_id==i,c(3:ncol(neo4j.genes))], function(x) ifelse(length(x)>1,na.locf(x),x))
  }
  # manual curation of genes with missing ensemble id
  ## 5446.1
  tmp <- list(chromosome="7", description="paraoxonase 3 [Source:EntrezGene;Acc:5446]",band="q21.3",gene_type="protein_coding",ensembl_gene="ENSG00000105852")
  neo4j.genes[neo4j.genes$label=="5446.1",names(tmp)] <- tmp
  ## 201288.1
  tmp <- list(chromosome="17", description="nitric oxide synthase 2 pseudogene 2 [Source:EntrezGene;Acc:5446]",band="p11.2",gene_type="pseudogene",ensembl_gene="ENSG00000167494")
  neo4j.genes[neo4j.genes$label=="201288.1",names(tmp)] <- tmp
  ## 645740.1
  tmp <- list(chromosome="17", description="nitric oxide synthase 2 pseudogene 1 [Source:EntrezGene;Acc:645740]",band="q11.2",gene_type="pseudogene",ensembl_gene="ENSG00000265788")
  neo4j.genes[neo4j.genes$label=="645740.1",names(tmp)] <- tmp
  ## 8781.1
  tmp <- list(chromosome="7", description="phosphoserine phosphatase pseudogene 1 [Source:EntrezGene;Acc:8781]",band="p11.2",gene_type="pseudogene",ensembl_gene="ENSG00000226278")
  neo4j.genes[neo4j.genes$label=="8781.1",names(tmp)] <- tmp
  ## 8041.1
  tmp <- list(chromosome="11", description="T-cell activation antigen p250 [DISCONTINUED; Source:EntrezGene;Acc:8041]",band="11pter-p11.2",gene_type="unknown", symbol="TP250")
  neo4j.genes[neo4j.genes$label=="8041.1",names(tmp)] <- tmp  
}


#######################
## METABOLITES
####################### 

## MAKE IDS
neo4j.metabolites <- data.metabolites %>% 
  #    mutate(sid = paste0(id,"_",model.suffix)) %>% 
  separate(id, into = c("mid", "no"), sep = -2, remove = F) %>% 
  select(-no) %>% mutate_all(as.character)

if(ext.ids){
  ## RENAME COLUMNS
  ids.metabo <- ids.metabo %>%
    rename(
      mid = abbreviation,
      name = fullName,
      formula = chargedFormula,
      kegg = keggId,
      pubchem = pubChemId,
      chebi = cheBlId
    ) %>% mutate_all(as.character)
  
  ## JOIN WITH EXTERNAL IDS
  neo4j.metabolites <- left_join(neo4j.metabolites, 
                                 ids.metabo,
                                 by = c("mid" = "mid")) %>% 
    process_double() %>% rename(sid=id) %>% mutate_all(trimws) %>% unique()
}else{
  neo4j.metabolites <- neo4j.metabolites %>% rename(sid=id) %>% mutate_all(trimws) %>% unique()
}


#######################
## REACTIONS
####################### 

if(ext.ids){
  ## RENAME COLUMNS
  ids.reaction<- ids.reaction%>% 
    rename(
      id = abbreviation,
      name = description,
      ec.code = ecnumber,
      kegg = keggId
    ) %>% 
    mutate_all(as.character) %>%  
    # messy
    mutate_at(vars(ec.code), ~map(., function(x) unlist(strsplit(x,split = ", ")))) %>% 
    mutate_at(vars(keggorthology,cog), ~map(., function(x) unlist(strsplit(x,split = " ")))) %>% 
    mutate_all(~map_chr(.,function(x) ifelse(all(!is.na(x)),paste(x,collapse = "|"),x)))
  
  ## JOIN WITH EXTERNAL IDS
  neo4j.reactions <- data.reactions %>%
    select(-products,-reactants,-reactants_stoi,-products_stoi, -genes) %>% 
    mutate_all(~map_chr(.,function(x) ifelse(length(unlist(x)!=0),ifelse(unlist(x)!="",paste(x,collapse = "|"),NA),NA))) %>% 
    mutate(rev=map_dbl(reversible, function(x)ifelse((x),1,0))) %>% left_join(.,ids.reaction,by="id") %>% 
    process_double() %>%# left_join(.,data.reactions %>% select(id,products,reactants,reactants_stoi,products_stoi)) %>% 
    mutate(sid = paste0(id,"_",model.suffix))
}else{
  neo4j.reactions <- data.reactions %>%
    select(-products,-reactants,-reactants_stoi,-products_stoi, -genes) %>% 
    mutate_all(~map_chr(.,function(x) ifelse(length(unlist(x)!=0),ifelse(unlist(x)!="",paste(x,collapse = "|"),NA),NA))) %>% 
    mutate(rev=map_dbl(reversible, function(x)ifelse((x),1,0))) %>%# left_join(.,data.reactions %>% select(id,products,reactants,reactants_stoi,products_stoi)) %>% 
    mutate(sid = paste0(id,"_",model.suffix))
}



## MAKE EDGES
# REACTIONS - METABO
neo4j.edges.r <- data.reactions %>%
  transmute(metBiGG = reactants, reBiGG = id, stoichiometry=reactants_stoi) %>%
  unnest(metBiGG,stoichiometry) %>% mutate(stoichiometry=stoichiometry*-1) %>%
  mutate(reBiGG = paste0(reBiGG,"_",model.suffix)) %>%
  as.data.frame() %>%
  select(metBiGG, reBiGG, everything()) 

neo4j.edges.p <- data.reactions %>%
  transmute(metBiGG = products, reBiGG = id, stoichiometry=products_stoi) %>%
  unnest(metBiGG,stoichiometry) %>%
  mutate(reBiGG = paste0(reBiGG,"_",model.suffix)) %>%
  as.data.frame() %>%
  select(metBiGG, reBiGG, everything()) 

neo4j.met.re <- bind_rows(neo4j.edges.r, neo4j.edges.p) #%>% 
#mutate(metBiGG = paste0(metBiGG,"_",model.suffix))

# REACTIONS - GENE
neo4j.re.gene <- data.reactions %>%
  transmute(gene=genes.orig,reBiGG = id) %>%
  unnest(gene) %>%
  as.data.frame() %>%
  select(gene, reBiGG, everything()) %>% 
  mutate(reBiGG = paste0(reBiGG,"_",model.suffix))


################################################################################
## WRITE 
################################################################################
out.name <- strsplit(basename(file.rdata),"\\.")[[1]][1]

write.csv(file = paste0(dirname(file.rdata),"/",out.name,"_metabolites.csv"),neo4j.metabolites,row.names = F,quote = T,na="")

write.csv(file = paste0(dirname(file.rdata),"/",out.name,"_reactions.csv"),neo4j.reactions,row.names = F,quote = T,na="")

write.csv(file = paste0(dirname(file.rdata),"/",out.name,"_genes.csv"),neo4j.genes,row.names = F,quote = T,na="")

write.csv(file = paste0(dirname(file.rdata),"/",out.name,"_met_re.csv"),neo4j.met.re,row.names = F,quote = T,na="")

write.csv(file = paste0(dirname(file.rdata),"/",out.name,"_re_gene.csv"),neo4j.re.gene,row.names = F,quote = T,na="")

#write.csv(file = "~/Documents/ICB/neo4j/data/fixes/uniformize_ids.csv", mid_replacements,row.names = F,quote = T,na="")


# ### GET DUPLICATE METABO IDS AND FIX NAMES (adapted script from parviz)
# 
# # Get metabolite names that are essentially the same,
# # minus the punctiations
# names_similar <- neo4j.metabolites %>% 
#   distinct(name) %>% 
#   mutate(name_simple = str_replace_all(name, "[[:punct:]+ ]", ""),
#          name_simple = tolower(name_simple)) %>% 
#   group_by(name_simple) %>% 
#   mutate(n = n()) %>% 
#   filter(n>1) %>% 
#   group_by(name_simple) %>% 
#   nest(name) %>% 
#   mutate(data = unlist(data, recursive = FALSE)) %>% 
#   rowwise() %>% 
#   transmute(name1 = data[1],
#             name2 = data[2])
# # Fix neo4j.metabolites
# # Uniformize metabolite names that are essentially the same
# neo4j.metabolites <- neo4j.metabolites %>% 
#   left_join(names_similar, by = c("name" = "name1")) %>% 
#   mutate(name = if_else(!is.na(name2), name2, name)) %>% 
#   dplyr::select(-name2)
# 
# 
# # Find metabolites that have more than one recon mid
# same_names <- neo4j.metabolites %>% 
#   distinct(mid, name, formula) %>% 
#   dplyr::count(name,formula) %>% 
#   filter(n>1) %>% .$name 
# 
# # Select a representative mid for metabolites with multiple mids.
# # Do this based on semantic similarity
# mid_replacements <- neo4j.metabolites %>% 
#   filter(name %in% same_names) %>% 
#   distinct(mid, name) %>% 
#   mutate(sim = stringdist(mid, name),
#          # artificially increase distance for metabolites
#          # with M0 at the beginning, since I know they are
#          # assigned names for presumably unknown metabolites
#          sim = if_else(grepl("^M0", mid), 1000, sim)) %>% 
#   group_by(name) %>%
#   mutate(ranking  = rank(sim, ties.method = "first"),
#          replacement = mid[which(ranking == min(ranking))]) %>% 
#   filter(mid != replacement) %>% 
#   ungroup() %>% 
#   dplyr::select(mid, replacement)
# uniformize the metabolite mids that actually point to the same
# metabolite
# data_metabolites3 <- neo4j.metabolites %>% 
#   left_join(mid_replacements, by = "mid") %>% 
#   mutate(mid = if_else(!is.na(replacement), replacement, mid)) %>% 
#   dplyr::select(-replacement) %>% 
#   mutate(id = glue::glue("{mid}_{compartment.id}"))


################
##maybe add references??