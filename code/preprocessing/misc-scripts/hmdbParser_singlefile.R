# This script downloads the most recent version of HMDB, parses all XML files 
# in there and writes out a tabular Excel file. 
# 
# Warning: Will take several hours.
# Script does not clean up, i.e. you might have to delete the downloaded file.
#
# Jan K
# last update: 2017-02-16
#
# !!! adaptation of this script to the version of HDMB that contains only a single zipped XML


# Attributes that are NOT imported from the XML files at this point:
#
# taxonomy
# ontology
# experimental_properties
# predicted_properties
# spectra
# biofluid_locations
# tissue_locations>
#   normal_concentrations
# abnormal_concentrations
# diseases
# wikipidia
# nugowiki
# metagene
# synthesis_reference
# general_references
# protein_associations


#### script settings ----

hmdbfileurl = "http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip"
hmdbfile = "hmdb_metabolites.zip"
downloadfile = T # can be set to F if file already exists in current directory
parseTaxonomy = T # can be set to T if Chemical taxonomy should be parsed
prepNeo4j = T #renames colums + adds one

#### initialize ----
library(XML)
library(xml2)
library(tictoc)
#library(xlsx)
library(openxlsx)
library(curl)

printf <- function(...)cat(sprintf(...))
timestamp = function() format(Sys.time(), "%H:%M:%S")


#### download hmdb file ----
printf('%s - downloading XML...\n', timestamp())
if (downloadfile) {
  curl_download(hmdbfileurl,hmdbfile,quiet=F)
}

# read XML
printf('%s - parsing XML...\n', timestamp())
con = unz(hmdbfile, filename='hmdb_metabolites.xml')
X = as_list(read_xml(con))



############# FIX ABOVE!! X needs to exist as list



printf('%s - processing\n', timestamp())

# these are values that are just directly copied 1:1 from the XML file

directcopy =
  c('chemical_formula', 'average_molecular_weight', 'monisotopic_molecular_weight', 'iupac_name',
    'traditional_iupac', 'cas_registry_number','smiles', 'inchi', 'inchikey', 'drugbank_id',
    'phenol_explorer_compound_id', 'phenol_explorer_metabolite_id', 'foodb_id', 'knapsack_id',
    'chemspider_id', 'kegg_id', 'biocyc_id', 'bigg_id', 'metlin_id', 'pubchem_compound_id',
    'het_id', 'chebi_id')


directTaxonomy=
 c('description','direct_parent','kingdom','super_class','class','sub_class','molecular_framework')


collectlist=list()
n = length(X)
upto = n
#upto = 50

for (i in 1:upto) {
  
  # extract entry L
  L = X[[i]]
  
    # assemble data frame row, one by one
    if (!is.list(L$pathways[[1]])) {
      
      # no pathway annotations
      newrow = data.frame(
        accession = L$accession[[1]],
        name = L$name,
        # synonyms = paste0(unlist(L$synonyms),collapse="|"),
        pathways_smpdb = '',
        pathways_KEGG = ''
      )
    } else {
      # pathway annotations exist
      newrow = data.frame(
        accession = L$accession[[1]],
        name = L$name,
        #  synonyms = paste0(unlist(L$synonyms),collapse="|"),
        pathways_smpdb = paste0(unlist(sapply(L$pathways, function(x)x$smpdb_id)),collapse="|"),
        pathways_KEGG = paste0(unlist(sapply(L$pathways, function(x)x$kegg_map_id)),collapse="|")
      )
    }
    
    colnames(newrow)[2] = "name"
    # copy all others
    s=sapply(directcopy, function(x){u=unlist(L[[x]]);if(is.null(u)){""}else{u}})
    # add to data frame
    newrow = cbind(newrow,t(as.data.frame(s)))
    
    #add synonyms
    if(any(grepl("\n",L$synonyms,perl=T))){
        # no secondary_accessions annotations
        newrow = cbind(newrow,synonyms='')
    }else{
        # secondary_accessions exist
        newrow = cbind(newrow,synonyms=paste0(unlist(L$synonyms),collapse="|"))
    }
    
   
    # add secondary accessions
    #added by MW
    if(any(grepl("\n",L$secondary_accessions,perl=T))){
        # no secondary_accessions annotations
        newrow = cbind(newrow,secondary_accessions='')
    }else{
        # secondary_accessions exist
        newrow = cbind(newrow,secondary_accessions=paste0(unlist(L$secondary_accessions),collapse="|"))
    }

    # parse and add taxonomy to data frame row
    #added by MW
    if (parseTaxonomy){
        
        L=L[['taxonomy']]
        
        t=sapply(directTaxonomy, function(x){u=unlist(L[[x]]);if(is.null(u)){""}else{u}})
        
        newrow = cbind(newrow,t(as.data.frame(t)))
        
        #check if alternative_parents exist
        if(any(grepl("\n",unlist(L$alternative_parents),perl=T))){
            newrow = cbind(newrow,alternative_parents='')
        }else{
            newrow = cbind(newrow,alternative_parents = paste0(unlist(L$alternative_parents),collapse="|"))
        }
        
        #check if substituents exist
        if(any(grepl("\n",unlist(L$substituents),perl=T))){
            newrow = cbind(newrow,substituents ='')
        }else{
            newrow = cbind(newrow,substituents=paste0(unlist(L$substituents),collapse="|"))
        }
        
    }



    # add to full table
    collectlist[[i]] = newrow

    # status
    if (i%%1000==0) {
      printf('%s - %d of %d\n', timestamp(), i, n)
    }

}

s=toc()
printf('%.2f per second\n', upto/(s$toc-s$tic))
printf('%.2fs estimated for all %d files\n',  n/(upto/(s$toc-s$tic)), n )

# merge into one big data frame
fulltable = do.call(rbind,collectlist)

###prep for Neo4j
#fulltable$synonyms=gsub("\n  ","",fulltable$synonyms)

#### export to excel ----

printf('Exporting Excel\n', timestamp())

L = X[[upto]]
# generate filename, take version info from last metabolite still in memory
outname = sprintf("hmdb_v%s.csv", L$version[[1]])

if(prepNeo4j){
    colnames(fulltable)[which(colnames(fulltable)=="accession")] <- "sid"
    colnames(fulltable)[which(colnames(fulltable)=="secondary_accessions")] <- "sa"
    fulltable$acc <- with(fulltable, paste(sid, sa, sep="|"))
}

# export
tic()
#write.xlsx(fulltable, file=outname, row.names=F)
write.csv(fulltable,file=outname,row.names=F)
printf('Done\n', timestamp())
toc()

