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
data.hgnc <- read_delim(file.in, delim = "\t")# colClasses = list(character = c(32,48,41,33,39,35)))

## REMOVE WITHDRAWN
data.hgnc <- data.hgnc %>%
    filter(name != "entry withdrawn")

## ONLY UPPER CASE
data.hgnc <- data.hgnc %>%
    mutate(symbol = toupper(symbol))

################################################################################
## HANDLE SYNONYMS
################################################################################
## SYNONYMS
data.hgnc.syn <- data.hgnc %>%
    filter(!is.na(prev_symbol)) %>%
    unnest(synonym = strsplit(prev_symbol, "\\|")) %>%
    select(hgnc_id,synonym)
## ALTERNATVE NAMES
data.hgnc.ali <- data.hgnc %>%
    filter(!is.na(alias_symbol)) %>%
    unnest(synonym = strsplit(alias_symbol, "\\|")) %>%
    select(hgnc_id, synonym)
## ACTUAL NAME
data.hgnc.self <- data.hgnc %>%
    transmute(hgnc_id, synonym = symbol)
## MERGE
data.hgnc.syn <- bind_rows(data.hgnc.self, data.hgnc.syn, data.hgnc.ali) %>%
    mutate(synonym = toupper(synonym)) %>%
    filter(!duplicated(.)) %>%
    filter(!duplicated(synonym)) %>%
    group_by(hgnc_id) %>%
    summarise(synonym = list(synonym))

## MERGE
data.hgnc <- data.hgnc %>%
    select(-alias_symbol, -prev_symbol) %>%
    left_join(data.hgnc.syn)

################################################################################
## WRITE
################################################################################
saveRDS(data.hgnc, file.out)
