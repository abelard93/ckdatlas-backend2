require(tidyverse)

##
##
## PARAMS
args = commandArgs(trailingOnly = T)
file.in  <- args[1]
file.out <- args[2]
##
## file.in <- "../raw/genes/hgnc_complete_set.txt"
##
################################################################################
## READ
################################################################################
data.entrez <- read_delim(file.in, delim = "\t", na = "-")# colClasses = list(character = c(32,48,41,33,39,35)))

## ONLY UPPER CASE
data.entrez <- data.entrez %>%
  mutate_at(vars(Symbol, Synonyms, Symbol_from_nomenclature_authority),toupper)

## FORMAT 
data.entrez <- data.entrez %>% 
  unite(col = symbols, Symbol,Synonyms,Symbol_from_nomenclature_authority, sep = "|", na.rm=T) %>% 
  transmute(sid = GeneID, 
            symbols, 
            chromosome, 
            location = map_location, 
            type_of_gene, 
            description, 
            name = Full_name_from_nomenclature_authority, 
            name_other = Other_designations,
            feature_type = Feature_type,
            tax_id = `#tax_id`,
            ensembl = map_chr(dbXrefs, function(x)ifelse(grepl("Ensembl", x, ignore.case = T), gsub('.*Ensembl:(ENSG[0-9]+).*', "\\1", x), NA))
            )

## REMOVE DUPLICATE SYMBOLS
data.entrez <- data.entrez %>%
  mutate(symbols = map_chr(
    symbols,
    ~ strsplit(., "\\|") %>% unlist %>% unique() %>% paste0(., collapse = "|")
  ))

################################################################################
## WRITE
################################################################################
write.csv(x = data.entrez, file = file.out, quote = T,na="",row.names = F)
