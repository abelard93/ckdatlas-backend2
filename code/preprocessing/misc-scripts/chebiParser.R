require(ontologyIndex)
require(plyr)
require(data.table)


rm(list = ls())

wd="/home/icb/maria.woerheide/Documents/chebi/"
setwd(wd)

fix(get_ontology)
#change line 7 to
#m <- regexpr(text = raw_lines, pattern = "^([^!]+[^! \t])")
chebiT <- get_ontology("chebi.obo", extract_tags = "everything")
chebi <- chebiT
#save(chebi,file="chebiOntology.Rdata")
load("chebiOntology.Rdata")
#CHEBI:53308
#166

tmp=list()

for (i in names(chebi)){
  if(is.list(chebi[[i]])){
    tmp[[i]] <-t(sapply(chebi[[i]], paste0, collapse="|"))
  }else{
    tmp[[i]] <- t(chebi[[i]])
  }
}

result <- do.call(rbind,tmp)
result <- t(result)

colnames(result) <- names(chebi)

DTresult <- as.data.table(result)
DTresult[DTresult==""] <- NA
#drop cols with only na 
to_drop <- names(DTresult)[apply(DTresult,2,function(x){all(is.na(x))})]
DTresult[,eval(to_drop):=NULL]

save(DTresult,file = "chebiDT.Rdata")


pattern <- paste0("(?<=http://purl.obolibrary.org/obo/chebi/)*+(?=xsd:)")

properties = list()

for (i in 1:length(DTresult[,property_value])) {
  x <- DTresult[i,property_value]
  if (!is.na(x)) {
    entries <- unlist(strsplit(x, split = "\\|"))
    for (e in entries){
      tmp <- gsub("http://purl.obolibrary.org/obo/chebi/|\" xsd:.*", "", e)
      tmp <- strsplit(tmp, split = "\"")
      prop <- trimws(tmp[[1]][1])
      value <- trimws(tmp[[1]][2])
      names(value) <- DTresult[i,id]
      properties[[prop]] <- c(properties[[prop]], value)
      
    }
  }
}
save(properties,file = "properties.Rdata")

load("chebiOntology.Rdata")
load("properties.Rdata")
load("chebiDT.Rdata")


for (i in names(properties)){
  
  na_ids <- DTresult[(!DTresult[,id]%in%names(properties[[i]])),id]
  na_ids <- setNames(rep(NA,length(na_ids)),na_ids)
  properties[[i]] <- c(properties[[i]],na_ids)
  tmp <-data.table(prop=properties[[i]],id=names(properties[[i]]))
  properties[[i]] <- tmp[ , .(prop = paste(prop, collapse="|")), by = id]
  setnames(properties[[i]], "prop", i)
  
  
}
props <-Reduce(merge,properties)
props <-Reduce(function (...) merge(...,by="id"),list(DTresult,props))
DF <-  as.data.frame(props[match(chebi[["id"]],props[,id]),])
DF[DF=="NA"] <- NA
save(DF,file = "chebi_ontology180124.Rdata")
load("chebi_ontology180124.Rdata")



##prepare for neo4j import
colnames(DF)[1] <- "sid"
DF <- DF[-which(!startsWith(DF$sid,"CHEBI")),]
write.csv(DF,file="../../neo4j/data/metChEBI/chebi_ontology180124.csv",row.names=F,quote = T,na = "")

#DF[100675,"synonym"] <- "\"polyglycidol\" RELATED [ChEBI]|\"poly(glycidol)\" RELATED [SUBMITTER]|\"PGLD\" RELATED [SUBMITTER]|\"Polyglycidol\" RELATED [ChemIDplus]|\"Oxiranemethanol homopolymer\" RELATED [ChemIDplus]"

