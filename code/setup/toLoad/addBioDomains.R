
## BioDomains
addBioDomains <- function() {
  for (file in getFiles("bioDomains")) {
    #read info on dataset
    #source <- getInfo(file, c("sid", "file_name", "node_type"), err)
    #skip file if no/incomplete info
    #if (!all(is.na(source))) {
      
      #tic(paste0("Added BioDomain file ", source["file_name"]))
      #check/add source node for dataset
      #checkSource(source, log)
    tic(paste0("Added BioDomain file ", basename(file)))
    
      bioDOMs  <- readRDS(file) %>%
        tidyr::unnest(ensembl_id) %>%
        mutate(n_total=ensembl_id %>% unique%>% length) %>%
        group_by(Biodomain) %>%
        summarize(
          GOterm = list(sort(unique(GOterm_Name))),
          ensembl_id = list(sort(unique(ensembl_id))),
          n_total = unique(n_total))
    
      for (dom in bioDOMs$Biodomain){
        
        query <- paste("CALL apoc.periodic.iterate(
                    \"UNWIND $properties AS prop return prop as geneID\",
                     \"MATCH (g:gene {ensembl_id:geneID}) SET g:BIODOM_",str_replace_all(dom," ","_"),"\", {batchSize:500, parallel:true, iterateList:true, concurrency:8, params: {properties:$properties}})
                     ",sep="")
        
        runQuery(query, periodic = TRUE, param = TRUE, param_list = bioDOMs$ensembl_id[bioDOMs$Biodomain==dom] %>% unlist)
      }
     
      info(log, capture.output(toc()))
      
  #  }
  }
}


