zap()
require(tidyverse)
require(xml2)
## GRAPHS
require(igraph)
## PARAMS
#args = commandArgs(trailingOnly = T)
#file.recon.sbml <- args[1]
#file.hgnc       <- args[2]
#file.out        <- args[3]
##
file.recon.sbml <- "~/Documents/ICB/neo4j/data/recon3d/Recon3D_301.xml"
file.hgnc       <- "~/Documents/ICB/neo4j/data/recon3d/genes_hgnc.rds"
file.out        <- "~/Desktop/recon3d.Rdata" 
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
                           function(x)ifelse(grepl("InChI"  , x, ignore.case = T), gsub('.*(InChI[^!"]+).*', "\\1", x), NA)),
         inchi2  = map_chr(annotation,
                           function(x)ifelse(grepl("InChI"  , x, ignore.case = T), gsub('.*(InChI[^!"]+).*', "\\1", x), NA))) %>%
  mutate(inchi = ifelse(is.na(inchi1), inchi2, inchi1)) %>%
  select(-inchi1, -inchi2) %>%
  select(-notes, -annotation)

## FORMAT
data.metabolites <- data.metabolites %>%
  select(-metaid) %>%
  mutate(charge = as.numeric(charge)) %>%
  rename(formula = chemicalFormula) 


################################################################################
## PARSE GENES
################################################################################
data.genes <- map_dfr(a$listOfGeneProducts, function(x)attributes(x)[names(attributes(x)) != "names"])
data.genes <- data.genes %>%
  mutate(entrez_id = as.integer(str_replace(label, "\\..*", ""))) %>%
  left_join(data.hgnc %>% transmute(entrez_id, symbol))


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
         products  = map(products , "species"))

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
  }))

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
  mutate(confidence = str_extract_all(notes, "Confidence Level:\\s*[0-9]+"),
         confidence = map(confidence, str_replace_all, pattern = "Confidence Level:\\s*", replacement = "")) 


## ADD SUBSYSTEMS
data.reactions <- data.reactions %>%
  left_join(data.subsystems %>% transmute(id = members, subsystem = name))

## FORMAT
data.reactions <- data.reactions %>%
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
## EDGES
################################################################################
data.edges.r <- data.reactions %>%
  transmute(v1 = reactants, v2 = id, reversible, fast) %>%
  unnest(v1) %>%
  as.data.frame() %>%
  select(v1, v2, everything()) 
data.edges.p <- data.reactions %>%
  transmute(v1 = id, v2 = products, reversible, fast) %>%
  unnest(v2) %>%
  as.data.frame() %>%
  select(v1, v2, everything()) 

data.edges <- bind_rows(data.edges.r, data.edges.p)


################################################################################
## VERTICES
################################################################################
data.vertices.m <- data.metabolites %>%
  transmute(id, type = "metab", label = name, compartment, charge) %>%
  as.data.frame()
data.vertices.r <- data.reactions %>%
  transmute(id, type = "reaction", label = name, confidence, subsystem,
            genes = map_chr(genes.all, paste, collapse = "|"), n.alt) %>%
  as.data.frame() 

data.vertices <- bind_rows(data.vertices.m,
                           data.vertices.r)

data.vertices <- data.vertices %>%
  mutate(shape = ifelse(type == "metab", "circle", "rectangle"),
         color = ifelse(type == "metab", "lightblue", "grey"),
         minor = FALSE)

################################################################################
## MINORS
################################################################################
# minors <- c("M_h", "M_h2o", "M_o2", "M_hco3", "M_h2o2", "M_co2",
#             "M_atp", "M_adp", "M_amp", "M_udp",
#             "M_dag_hs",
#             "M_na1", "M_nad", "M_nadh", "M_nadp", "M_nadph",
#             "M_fad", "M_fadh", "M_fadh2",
#             "M_pi", "M_ppi", "M_coa", "M_accoa")
# 
# for(m in minors){
#     for(c in data.compartments$id){
#         ## m = minors[1]
#         ## c = "c"
#         x = paste0(m, "__91__", c, "__93__")
#         ## UPDATE EDGES WITH RUNNIGN NUMBER
#         if(x %in% data.edges$v1){
#             data.edges.mod.v1 <- data.edges %>% filter(v1 == x) %>%
#                 mutate(v1 = paste(v1, 1:n(), sep ="_"))
#         }else{
#             data.edges.mod.v1 = data.frame()
#         }
#         if(x %in% data.edges$v2){
#             data.edges.mod.v2 <- data.edges %>% filter(v2 == x) %>%
#                 mutate(v2 = paste(v2, (nrow(data.edges.mod.v1) + 1):(nrow(data.edges.mod.v1) + n()), sep ="_"))
#         } else {
#             data.edges.mod.v2 = data.frame()
#         }
#         ##
#         ## UPDATE VERTICES WITH RUNNING NUMBER
#         if(nrow(data.edges.mod.v1) + nrow(data.edges.mod.v2) > 0){
#             data.vertices.new <- data.vertices %>% filter(id == x)
#             data.vertices.new <- data.vertices.new[ rep(1, nrow(data.edges.mod.v1) + nrow(data.edges.mod.v2)),] %>%
#                 mutate(id    = paste(id, 1:n(), sep = "_"),
#                        minor = T)
#             ##
#             ## UPDATE EDGES AND VERTICES LIST
#             data.edges.mod <- bind_rows(data.edges.mod.v1, data.edges.mod.v2)
#             data.edges.org <- data.edges %>% filter(v1 != x) %>% filter(v2 != x)
#             data.edges <- bind_rows(data.edges.mod, data.edges.org)
#             ##
#             data.vertices <- rbind(data.vertices %>% filter(id != x),
#                                    data.vertices.new)
#         }
#         ## CHECK
#         ## vs <- c(as.character(data.edges$v1), as.character(data.edges$v2))
#         ## unique(vs) %without% data.vertices$id
#     }
# }

################################################################################
## CREATE GRAPH
################################################################################
## reversible edges
data.edges.rev <- data.edges %>% filter(reversible) %>%
  mutate(tmp = v1) %>%
  mutate(v1 = v2) %>%
  mutate(v2 = tmp) %>%
  select(-tmp)
data.edges <- bind_rows(data.edges, data.edges.rev)

## CREATE GRAPH
graph <- graph_from_data_frame(data.edges, directed = TRUE, vertices = data.vertices)
##d = degree(graph)



################################################################################
## WRITE - for neo4j integration
################################################################################


neo4j.genes <- data.genes %>% group_by(entrez_id) %>% summarise(recon3d_id=paste(unique(id), collapse="|"),symbol=paste(unique(symbol), collapse="|"))

neo4j.metabolites <- data.metabolites %>% mutate(sid=id) %>% select(-id) %>% select(sid, everything())

##maybe add references??
neo4j.reactions <- data.reactions %>%
  transmute(sid=id, name, subsystem, reversible, rev=map_dbl(reversible, function(x)ifelse((x),1,0)), fast, genes = map_chr(genes.all, paste, collapse = "|"),genes.essential = map_chr(genes.essential, paste, collapse = "|")) %>%
  as.data.frame() 


neo4j.edges.r <- data.reactions %>%
  transmute(metBiGG = reactants, reBiGG = id, stoichiometry=reactants_stoi) %>%
  unnest(metBiGG,stoichiometry) %>% mutate(stoichiometry=stoichiometry*-1) %>%
  as.data.frame() %>%
  select(metBiGG, reBiGG, everything()) 

neo4j.edges.p <- data.reactions %>%
  transmute(metBiGG = products, reBiGG = id, stoichiometry=products_stoi) %>%
  unnest(metBiGG,stoichiometry) %>%
  as.data.frame() %>%
  select(metBiGG, reBiGG, everything()) 

neo4j.met.re <- bind_rows(neo4j.edges.r,
                          neo4j.edges.p)

neo4j.re.gene <- data.reactions %>%
  transmute(gene=genes.all,reBiGG = id) %>%
  unnest(gene) %>%
  as.data.frame() %>%
  select(gene, reBiGG, everything()) 

write.csv(file = "~/Documents/ICB/neo4j/data/recon3d/recon3d_metabolites.csv",neo4j.metabolites,row.names = F,quote = T,na="")

write.csv(file = "~/Documents/ICB/neo4j/data/recon3d/recon3d_reactions.csv",neo4j.reactions,row.names = F,quote = T,na="")

write.csv(file = "~/Documents/ICB/neo4j/data/recon3d/recon3d_genes.csv",neo4j.genes,row.names = F,quote = T,na="")

write.csv(file = "~/Documents/ICB/neo4j/data/recon3d/recon3d_met_re.csv",neo4j.met.re,row.names = F,quote = T,na="")

write.csv(file = "~/Documents/ICB/neo4j/data/recon3d/recon3d_re_gene.csv",neo4j.re.gene,row.names = F,quote = T,na="")
################################################################################
## WRITE
################################################################################
cat("SAVE\n")
save(data.metabolites,
     data.reactions,
     graph,
     file = file.out)
