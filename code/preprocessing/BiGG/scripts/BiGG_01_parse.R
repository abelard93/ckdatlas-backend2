####
# Script from Jonas Zierer to parse recon3D
# adapted by Maria Woerheide for integration into Neo4j
# last modified: 2019-04-16
####
zap()
require(tidyverse)
require(xml2)
require(stringdist)
library(zoo)
## PARAMS
args = commandArgs(trailingOnly = T)
file.hgnc       <- args[1]
file.recon.sbml <- args[2]
file.out        <- args[3]
##
#file.recon.sbml     <- "~/Documents/PHD/neo4j/neo4jserver/code/BiGG/data/raw/Recon3D/Recon3D_301.xml"
#file.hgnc           <- "~/Documents/PHD/neo4j/neo4jserver/code/BiGG/data/genes_hgnc.rds"
#file.out            <- "~/Documents/PHD/neo4j/neo4jserver/code/BiGG/data/raw/Recon3D_301.Rdata" 



# ## PARAMS
# args = commandArgs(trailingOnly = T)
# file.recon.sbml <- args[1]
# file.hgnc       <- args[2]
# file.out        <- args[3]
##
## file.recon.sbml <- "../data/raw/Recon3D_301.xml"
## file.hgnc       <- "../data/processed/genes_hgnc.rds"
##
################################################################################
## READ SBML DATA
################################################################################
## HGNC
data.hgnc    <- readRDS(file.hgnc)
gene.symbols <- setNames(data.hgnc$symbol, data.hgnc$hgnc_id)

## ## READ MODEL
model.xml <- read_xml(file.recon.sbml)
## first lebel is SBML
model.xml <- xml_child(model.xml)
## as list
a         <- as_list(model.xml)

################################################################################
## PARSE COMPARTMENTS
################################################################################
data.compartments <- map_dfr(a$listOfCompartments, function(x)attributes(x))


################################################################################
## PARSE METABOLITES
################################################################################
data.metabolites <- map_dfr(a$listOfSpecies     , function(x)attributes(x)[names(attributes(x)) != "names"]) %>%
  rename(compartment.id = compartment) %>%
  left_join(data.compartments %>% transmute(compartment.id = id, compartment = name))

## ADD NOTES & ANNOTATION
## unique(unlist(map(map(a$listOfSpecies, function(x)attributes(x)), "names")))
data.metabolites <- data.metabolites %>%
  mutate(notes      = map(a$listOfSpecies, "notes")) %>%
  mutate(notes      = map_chr(notes, function(n)unname(paste(unlist(n), collapse = "###", sep = "###")))) %>%
  mutate(annotation = map(a$listOfSpecies, ~.x$annotation$RDF$Description$is$Bag)) %>%
  mutate(annotation = map(annotation     , function(l)map_chr(l, attr, "resource"))) %>%
  mutate(annotation = map_chr(annotation, paste, collapse = "###", sep = "###")) %>%
  ## ADD IDS - changed to include hmdb id (Maria)
  mutate(kegg    = map_chr(annotation,
                           function(x)ifelse(grepl("kegg"   , x, ignore.case = T), gsub('.*kegg.compound/(C[0-9]+).*', "\\1", x), NA)),
         hmdb    = map_chr(annotation,
                           function(x)ifelse(grepl("hmdb"   , x, ignore.case = T), gsub('.*hmdb/(HMDB[0-9]+).*', "\\1", x), NA)),
         pubchem = map_chr(annotation,
                           function(x)ifelse(grepl("pubchem", x, ignore.case = T), gsub('.*pubchem.compound/([0-9]+).*', "\\1", x), NA)),
         inchi1  = map_chr(notes,
                           function(x)ifelse(grepl("InChI"  , x, ignore.case = T), gsub('.*(InChI[^!"#]+).*', "\\1", x), NA)),
         inchi2  = map_chr(annotation,
                           function(x)ifelse(grepl("InChI"  , x, ignore.case = T), gsub('.*(InChI[^!"#]+).*', "\\1", x), NA))) %>%
  mutate(inchi = ifelse(is.na(inchi1), inchi2, inchi1)) %>%
  select(-inchi1, -inchi2) %>%
  select(-notes, -annotation)

## FORMAT
data.metabolites <- data.metabolites %>%
  select(-metaid) %>%
  mutate(charge = as.numeric(charge)) %>%
  rename(formula = chemicalFormula)  %>%
  mutate(id = str_replace(id, "__91__([a-z])__93__", "_\\1") %>% str_remove("^M_"))


################################################################################
## PARSE GENES
################################################################################
data.genes <- map_dfr(a$listOfGeneProducts, function(x)attributes(x)[names(attributes(x)) != "names"])
data.genes <- data.genes %>%
  mutate(entrez_id = map_chr(label,function (x) ifelse(grepl(x,pattern = "peg"),x,str_replace(x, "\\..*", "")))) %>%
  left_join(data.hgnc %>% transmute(entrez_id=as.character(entrez_id), symbol) %>% na.omit)


################################################################################
## PARSE SUBSYSTEMS
################################################################################
data.subsystems <-  map_dfr(a$listOfGroups, function(x)attributes(x)[names(attributes(x)) != "names"])
data.subsystems <- data.subsystems %>%
  mutate(members = map(a$listOfGroups, "listOfMembers")) %>%
  mutate(members = map(members, function(x)map_chr(x, function(x)attr(x, "idRef")))) %>%
  select(-sboTerm) %>%
  unnest(members)


################################################################################
## PARSE REACTIONS
################################################################################
unique(unlist(map(map(a$listOfReactions, function(x)attributes(x)), "names")))
data.reactions    <- map_dfr(a$listOfReactions   , function(x)attributes(x)[names(attributes(x)) != "names"]) %>%
  mutate(reversible = reversible == "true") %>%
  mutate(fast       = fast == "true") 
##
## REACTANTS & PRODUCTS
data.reactions <- data.reactions %>%
  mutate(reactants = map(a$listOfReactions, function(x)map_dfr(x$listOfReactants, function(x)attributes(x))),
         products  = map(a$listOfReactions, function(x)map_dfr(x$listOfProducts, function(x)attributes(x)))) %>%
  ## STOI
  mutate(reactants_stoi = map(reactants, "stoichiometry"),
         products_stoi  = map(products , "stoichiometry")) %>%
  mutate(reactants_stoi = map(reactants_stoi, as.integer),
         products_stoi  = map(products_stoi , as.integer)) %>%
  ## NAME ONLY (FOR NOW)
  mutate(reactants = map(reactants, "species"),
         products  = map(products , "species")) %>%
  ## IDS
  mutate(reactants = map(reactants, ~ .x %>% str_replace("__91__([a-z])__93__", "_\\1") %>% str_remove("^M_")),
         products  = map(products , ~ .x %>% str_replace("__91__([a-z])__93__", "_\\1") %>% str_remove("^M_")))

## GENE PRODUCTS
parseGene <- function(n, x){
  ## IF NULL --> STOP
  if(is.null(x)){
    return(list())
  }
  ## IF LEAVE --> RETURN
  if(n == "geneProductRef"){
    return(attr(x, "geneProduct"))
  }
  ## OR NODE --> NOTHING TO DO
  if(n == "or"){
    xx <- map(seq_along(x), function(n)parseGene(names(x)[n], x[[n]]))
    if(any(map_chr(xx, class) == "list"))
      xx <- unlist(xx, recursive = F)
    return(xx)
  }
  ## AND NODE --> DISTRIBUTE (IF ANY)
  if(n == "and"){
    ## DISTRIBUTE IF ORS
    if("or" %in% names(x)){
      ## PARSE CHILDREN
      x.parse <- map(seq_along(x), function(n)parseGene(names(x)[n], x[[n]]))
      ## ALL COMBINATION OF ORS
      x.ors   <- x.parse[ names(x) == "or" ]
      lengths <- map(x.ors, function(i)if(length(i) == 0){0} else {c(1:length(i))})
      ## hack
      x.ors   <- cross(lengths) %>%
        map(function(i){
          cc <- unlist(map(1:length(i), function(j){
            x.ors[[j]][[ i[[j]] ]]
          }))
        })
      ## ADD OTHER
      x.others <- unlist(x.parse[ names(x) != "or" ])
      x.ors    <- map(x.ors, function(a)c(a, x.others))
      ## RETURN
      return(x.ors)
    }else{
      return(map_chr(seq_along(x), function(n)parseGene(names(x)[n], x[[n]])))
    }
  }
  ## OTHERS
  warning("UNKNOWN NODE TYPE")
}

my.replicate <- function (n, expr, simplify = 'array') {sapply(n, eval.parent(substitute(function(...) expr)),simplify = simplify)}
unlockBinding("replicate", as.environment("package:base"))
assign("replicate", my.replicate, "package:base")

data.reactions <- data.reactions %>%
  mutate(genes.orig    = map(a$listOfReactions, function(x)x$geneProductAssociation)) %>%
  mutate(genes.parsed = map(genes.orig, ~parseGene(n=names(.x), .x[[1]])))

## FORMAT GENES - changed to map to entrez (Maria)
data.reactions <- data.reactions %>%
  mutate(genes = map(genes.parsed, function(gg){
    gn <- map(gg, function(kk)data.genes %>%
                filter(id %in% kk) %>%
                distinct(entrez_id) %>%
                filter(!is.na(entrez_id)) %>%
                arrange(entrez_id) %>%
                .$entrez_id)
    gn <- unique(gn)
    gn <- gn[ map_int(gn, length) > 0 ]
  })) 

## ADD GENE
data.reactions <- data.reactions %>%
  mutate(n.alt  = map_int(genes, length)) %>%
  #mutate(n.genes = map(genes  , ~map_int(.x, length))) %>%
  mutate(genes.all       = map(genes, function(g){
    if(length(g) == 0)
      return(integer())
    sort(unique(unlist(g)))
  })) %>%
  mutate(genes.essential = map(genes, function(g){
    if(length(g) == 0)
      return(integer())
    sort(Reduce(intersect, g))
  })) %>% 
  mutate(genes.parsed  = map(genes.parsed , ~ unlist(.x) %>% str_replace_all("__46__", ".") %>% str_remove("^G_")))


## ADD NOTES & ANNOTATION
data.reactions <- data.reactions %>%
  mutate(notes      = map(a$listOfReactions, ~.x$notes$body)) %>%
  mutate(notes      = map(notes            , ~paste(unlist(.x), collapse = "###"))) %>%
  mutate(annotation = map(a$listOfReactions, ~.x$annotation$RDF$Description$is$Bag)) %>%
  mutate(annotation = map(annotation     , function(l)map_chr(l, attr, "resource"))) %>%
  mutate(annotation = map_chr(annotation, paste, collapse = "###", sep = "###")) 

## ADD EXTERNAL IDs
data.reactions <- data.reactions %>%
  ## EC CODE
  mutate(annotation = map(annotation, str_replace_all, pattern = "ec-code/EC:", replacement = "ec-code/")) %>%
  mutate(ec.code    = str_extract_all(annotation, "ec-code/[0-9-]+\\.[0-9-]+\\.[0-9-]+\\.[0-9\\-]+"),
         ec.code    = map(ec.code, str_replace_all, pattern = "ec-code/", replacement = "")) %>%
  ## KEGG REACTION ID
  mutate(kegg       = str_extract_all(annotation, "kegg.reaction/R[0-9]+"),
         kegg       = map(kegg, str_replace_all, pattern = "kegg.reaction/", replacement = "")) %>%
  ## PUBMED ID
  mutate(pubmed     = str_extract_all(notes, "PMID\\s*[:-]*\\s*[0-9]+"),
         pubmed     = map(pubmed, str_replace_all, pattern = "PMID\\s*[:-]*\\s*", replacement = ""),
         pubmed     = map(pubmed, as.numeric)) %>%
  ## PUBMED CENTRAL ID
  mutate(pubmedC    = str_extract_all(notes, "PMCID\\s*[:-]\\s*PMC[0-9]+"),
         pubmedC    = map(pubmedC, str_replace_all, pattern = "PMCID\\s*[:-]\\s*", replacement = "")) %>%
  ## OMIM
  mutate(omim       = str_extract_all(notes, "OMIM\\s*[:-]\\s*[0-9]+")) %>%
  ## OMIM
  mutate(hmdb       = str_extract_all(notes, "HMDB[0-9]+")) %>%
  ## OTHER REF
  mutate(ref        = str_extract_all(notes, "References:.*"),
         ref        = map(ref, str_replace_all, pattern = "References:", replacement = ""),
         ref        = map(ref, str_replace_all, pattern = "###.*", replacement = ""),
         ref        = map(ref, str_replace_all, pattern = "PMID\\s*[:-]*\\s*[0-9]+\\s*[.,;]*", replacement = ""),
         ref        = map(ref, str_replace_all, pattern = "PMCID\\s*[:-]*\\s*PMC[0-9]+\\s*[.,;]*", replacement = ""), 
         ref        = map(ref, str_replace_all, pattern = "HMDB[0-9]+\\s*[.,;]*", replacement = ""), 
         ref        = map(ref, str_replace_all, pattern = "OMIM[:]*[0-9]+\\s*[.,;]*", replacement = ""), 
         ref        = map(ref, str_replace_all, pattern = "^\\s+", replacement = "")) %>%
  ## CONFIDENCE
  mutate(confidence = str_extract_all(notes     , "Confidence Level:\\s*[0-9]+"),
         confidence = ifelse(map_int(confidence, length) == 0, "", confidence),
         confidence = map_chr(confidence, str_replace_all, pattern = "Confidence Level:\\s*", replacement = ""),
         confidence = as.integer(confidence)) 


## ADD SUBSYSTEMS
data.reactions <- data.reactions %>%
  left_join(data.subsystems %>% transmute(id = members, subsystem = name))

## FORMAT
data.reactions <- data.reactions %>%
  mutate( id  = str_remove(id,"^R_"), genes.orig =genes.parsed) %>%
  select(id,
         reactants, reactants_stoi,
         products, products_stoi,
         name, subsystem,
         reversible, fast,
         lowerFluxBound, upperFluxBound,
         n.alt, genes.orig, genes, genes.all, genes.essential,
         ec.code, kegg, 
         pubmed, pubmedC, omim, hmdb, ref,
         confidence)

## REMOVE EXCHANGE/DEMAND REACTIONS (NO PRODUCTS)
data.reactions <- data.reactions %>%
  filter(map_int(products, length) > 0)

################################################################################
## WRITE
################################################################################
cat("SAVE\n")
save(data.metabolites,
     data.reactions,
     data.genes,
     file = file.out)


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