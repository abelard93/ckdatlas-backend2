library(foreach)
library(doParallel)

getSignRels <- function(update=F){
  
  # if true, delete all significant relationships + created html edge information
  if(update){
    
    tic("Deleted precomputed relationships and redoing")
    
    # delete tabix file of edges
    if(file.exists("CKDatlas/html-edge-labels.csv.gz")) file.remove("CKDatlas/html-edge-labels.csv.gz")
    if(file.exists("CKDatlas/html-edge-labels.csv.gz.tbi")) file.remove("CKDatlas/html-edge-labels.csv.gz.tbi")
    
    # remove tmp files
    file.remove(getFiles("CKDatlas/tmp"))
    if(dir.exists("CKDatlas/tmp/snp_to_gene_tmp")) file.remove(getFiles("CKDatlas/tmp/snp_to_gene_tmp"))
    
    # remove COABUNDANCE files
    if(dir.exists("CKDatlas/coabundance")) file.remove(getFiles("CKDatlas/coabundance"))
    
    # remove COEXPRESSION files
    if(dir.exists("CKDatlas/coexpression")) file.remove(getFiles("CKDatlas/coexpression"))
    
    # remove COREGULATION files
    if(dir.exists("CKDatlas/coregulation")) file.remove(getFiles("CKDatlas/coregulation"))
    
    # remove GENETIC_ASSOCIATION files
    if(dir.exists("CKDatlas/genetic_association")) file.remove(getFiles("CKDatlas/genetic_association"))
    
    # remove METAVOLIC_ASSOCIATION files
    if(dir.exists("CKDatlas/metabolic_association")) file.remove(getFiles("CKDatlas/metabolic_association"))
    
    # remove PARTIAL_CORRELATION files
    if(dir.exists("CKDatlas/partial_correlation")) file.remove(getFiles("CKDatlas/partial_correlation"))
    
    info(log,capture.output(toc()))
  }
  
  tic("Extracted significant relationships")
  
  file.info <- data.makepath("data/QTLs/eQTL/cis") # path to data repository
  file.output <- data.makepath("data/CKDatlas/")
  

  # 1. GTEX tissue ---------------
  
  if(!file.exists(paste0(file.output,"coregulation/significant_eQTL.csv"))){ # is file there?
    print("Test 1")
    ## SNP mapping --------------------------
    ## get all snp to gene mappings from db (only for protein coding genes)
    ## includes mapping from snipa and eQTL mapping
    
    if(file.exists(paste0(file.output,"tmp/snp_to_gene_mapping.csv"))){ 
      file.mapping <- fread(paste0(file.output,"tmp/snp_to_gene_mapping.csv"))
    }else{
      file.mapping <- getSnp2gene()
    }
    
    # get all proteinCoding genes in database
    proteinCoding <- "MATCH (p:proteinCoding) return distinct p.ensembl_id as id" %>% 
      runQuery(read = TRUE) %>% 
      .$id
    
    # set key on SNP id for fast search 
    setkey(file.mapping, sid)
    # summarize genes so that there is only one entry for each rsID
    file.mapping <- file.mapping[, .(gene_snp = paste(gene, collapse="|")), by = sid]

    
    ## COREGULATION edges --------------------------
    # write as we go -> large file !!
    header <- matrix(nrow = 1,c("gene_snp","gene","tissue","pvalue","label", "label_tmp", "source", "relationship"))
    write.table(header,file = paste0(file.output,"coregulation/significant_eQTL.csv"),sep = ",",col.names = F,row.names = F,append = F,quote = T)
    
    ## get tissue types 
    files_info <- c()
    for(file in grep(list.files(file.info,full.names = T,recursive = T),pattern = "_info.csv",invert = T,value = T)){
      tissue <- getInfo(file, c("sample_type"), err) %>% 
                .$sample_type
      files_info <- rbind(files_info,c(file,tissue,"cis"))
    }
    
    for(file in grep(list.files("~/local/neo4j/ADatlas/import/data/QTLs/eQTL/trans",full.names = T,recursive = T),pattern = "_info.csv",invert = T,value = T)){
      tissue <- getInfo(file, c("sample_type"), err) %>% 
                .$sample_type
      files_info <- rbind(files_info,c(file,tissue,"trans"))
    }
    
    files_info <- files_info %>% data.frame(stringsAsFactors = F)
    colnames(files_info) <- c("file","tissue","eQTLtype")
    
    ## loop through all tissues
    for(tissue in unique(files_info$tissue)){
      
      # files for cis 
      files_cis <- files_info %>% 
        filter(tissue==!!tissue & eQTLtype=="cis") %>% 
        .$file
      # files for trans
      files_trans <- files_info %>% 
        filter(tissue==!!tissue & eQTLtype=="trans") %>% 
        .$file
      
      # df for cis
      gtex.cis <-  data.frame(
                      sid=character(), 
                      gene_snp=character(), 
                      gene=character(), 
                      pvalue = numeric(), 
                      type=character(), 
                      source=character(),
                      stringsAsFactors = F
                    )
      
      # df for trans
      gtex.trans <- data.frame(
                      sid=character(), 
                      gene_snp=character(), 
                      gene=character(), 
                      pvalue = numeric(), 
                      type=character(), 
                      source=character(),
                      stringsAsFactors = F
                    )
      
      ## cis eQTL ------------------------------------
      if(length(files_cis)>0){
  
        ## read all files for cis eQTLs in this tissue
        for(file in files_cis){
          
          # read source info
          source <- getInfo(file, c("sid", "file_name"), err) %>% .$publication
          ## remove info on pmid from source to make it shorter
          source <- ifelse(source=="GTExConsortium_2019","GTEx_2019",source) %>% 
            sub('^([^_]+(?:_[^_]+){1}).*', '\\1', .)

          # TISSUE
          tissue_id <- tissue
          ## READ file, rename and filter genes to proteincoding
          # disregard edges to non-protein coding genes (=>reduces complexity & not part of atlas)
          tmp <- fread(file,stringsAsFactors = F, na.strings = c("","NA")) 
          
          ## rename 
          if("pvalue_nominal"%in%colnames(tmp)){
            tmp <- tmp %>% 
              dplyr::rename(pvalue = pvalue_nominal)
          }
          
          ## rename 
          if("rs_id_dbSNP151_GRCh38p7"%in%colnames(tmp)){
            tmp <- tmp %>% 
              dplyr::rename(rsid = rs_id_dbSNP151_GRCh38p7)
          }
          
          ## rename 
          if("gene_id"%in%colnames(tmp)){
            tmp <- tmp %>% 
              dplyr::rename(gene = gene_id)
          }
          
          # construct final df and filter genes
          tmp <- tmp %>% 
            transmute(gene, pvalue, rsid, type="cis", source=source) %>% 
            filter(gene%in%proteinCoding) %>% 
            data.table
          
          # add genes that snp is annotated to
          setkey(tmp, rsid)
          tmp <- file.mapping[tmp] 
        
          # add to result df
          gtex.cis <- rbind(gtex.cis,tmp)
        }
        
        # order columns
        gtex.cis <- gtex.cis %>% 
          select(sid, pvalue, gene, gene_snp, type, source)
        
        # unwind df so each gene_snp - gene has one row
        if(nrow(gtex.cis)>0){
          
          gtex.cis <- gtex.cis %>% 
            separate_rows(gene_snp) %>% 
            distinct %>%
            data.table()
        }
      }
  
      ## trans eQTL ------------------------------------
      if(length(files_trans)>0){
        
        for(file in files_trans){
          #read source info
          source <- getInfo(file, c("sid", "file_name"), err) %>% .$publication
          ## remove info on pmid from source to make it shorter
          source <- ifelse(source=="GTExConsortium_2019","GTEx_2019",source) %>% 
            sub('^([^_]+(?:_[^_]+){1}).*', '\\1', .)
          
          # TISSUE
          tissue_id <- tissue
          ## READ file, rename and filter genes to proteincoding
          # disregard edges to non-protein coding genes (=>reduces complexity & not part of atlas)
          tmp <- fread(file,stringsAsFactors = F, na.strings = c("","NA"))
          
          ## rename 
          if("pvalue_nominal"%in%colnames(tmp)){
            tmp <- tmp %>% 
              dplyr::rename(pvalue = pvalue_nominal)
          }
          
          ## rename 
          if("rs_id_dbSNP151_GRCh38p7"%in%colnames(tmp)){
            tmp <- tmp %>% 
              dplyr::rename(rsid = rs_id_dbSNP151_GRCh38p7)
          }
          
          ## rename 
          if("gene_id"%in%colnames(tmp)){
            tmp <- tmp %>% 
              dplyr::rename(gene = gene_id)
          }
          
          ## READ file
          tmp <- tmp %>% 
            transmute(gene, pvalue ,rsid, type="trans", source=source) %>% 
            filter(gene%in%proteinCoding) %>% 
            data.table
          
          setkey(tmp, rsid)
          tmp <- file.mapping[tmp] 
          
          gtex.trans <- rbind(gtex.trans,tmp)
        }

        #### parallel MAPPING of gene to snp ----------------------
        # for each entry/hit/row in file
        # add gene that has been annotated by SNIPA to the SNP and add GENE as snp_gene
        # why? => SNPs that are associated with a significant eQTL edge are also assumed to "map" to that gene
        gtex.trans <- gtex.trans %>% 
          select(sid, pvalue, gene, gene_snp, type, source) 
        
        if(nrow(gtex.trans)>0){
          gtex.trans <- gtex.trans %>% 
            separate_rows(gene_snp) %>% 
            distinct %>%
            data.table()
        }
      }
    
      gtex.summary <- rbind(gtex.cis,gtex.trans) #%>% head(gtex.summary,n=1000)
      
      rm(gtex.cis)
      rm(gtex.trans)
   
      ## make HTML label ---------------------------------------
      gtex.summary <- gtex.summary[,.(pvalue_label= paste0(signif(pvalue,3),collapse="<br>"),source_label=paste0(source,collapse="<br>"), pvalue=min(pvalue), source=list(source)),by=list(gene,gene_snp,sid,type)]
     
      gtex.summary <- gtex.summary[,.(label=paste0("<table class='snipa-plot-tooltip'><thead><tr><th colspan='3' id='coreg'>COREGULATION</th></tr></thead><tbody><tr><td>Min. p-value&nbsp;</td><td width='33%'>",signif(min(pvalue),3),"<br>(",sid[which.min(pvalue)],")</td></tr><tr><td>#SNPs&nbsp;</td><td>",.N,"</td><td width='33%'>&nbsp;</td></tr><tr><td>Source(s)&nbsp;</td><td>",paste0(unique(source %>% unlist),collapse="<br>"),"</td></tr><tr><th class='tissue' colspan='3'>Tissue</th></tr><tr><td class='chead' width='33%'>SNP&nbsp;</td><td class='chead' width='33%'>P-value</td><td class='chead' width='33%'>source</td></tr></tbody></table><table class='snipa-plot-tooltip'><tr><th class='sub' colspan='3'>",unique(..tissue_id),"</th></tr>",paste0("<td width='33%'>",sid," (",type,")&nbsp;</td><td width='33%'>",pvalue_label,"</td><td width='33%'>",source_label,"</td></tr>",collapse = ""),"</table>"), 
                                      label_tmp = paste0("<table class='snipa-plot-tooltip'><tr><th colspan='3' class='sub'>",unique(..tissue_id),"</th></tr>",paste0("<td width='33%'>",sid," (",type,")&nbsp;</td><td width='33%'>",pvalue_label,"</td><td width='33%'>",source_label,"</td></tr>",collapse = ""),"</table>"), 
                                      pvalue=min(pvalue),source=paste0(unique(source %>% unlist),collapse="|")),by=list(gene,gene_snp)]
      
      gtex.summary$tissue <- tissue_id
      gtex.summary$relationship <- paste0("COREGULATION_", tissue_id %>% toupper %>% str_replace_all(" ","_") %>% str_replace("-","_"))
      
      gtex.summary <- gtex.summary %>%
        select(gene_snp, gene, tissue, pvalue, label, label_tmp, source, relationship)
      
      ## WRITE 
      write.table(gtex.summary,file = paste0(file.output,"coregulation/significant_eQTL.csv"),sep = ",",col.names = F,row.names = F,append = T,quote = T)
      rm(gtex.summary)
      message(tissue_id)
    }
    
    rm(file.mapping)
  }

  gc() 
  
  ## 2. GTEX all ----------------------------------------

  file.output <- data.makepath("data/CKDatlas/")
  
  if(!file.exists(paste0(file.output,"coregulation/significant_eQTL_all.csv"))){ # is file there?
    print("Test 2")
    ## read file with tissue specific eQTLs
    file.gtex <- fread(paste0(file.output,"coregulation/significant_eQTL.csv"))  
    
    ## SUMMARIZE -> add html edge code
    file.gtex.all <- file.gtex[,label:=NULL]
    rm(file.gtex)

    file.gtex.all <- file.gtex.all[,.(label=paste0("<table class='snipa-plot-tooltip'><thead><tr><th colspan='3' id='coreg'>COREGULATION</th></tr></thead><tbody><tr><td>Min p-value&nbsp;</td><td width='33%'>",signif(min(pvalue),3),"<br>(",tissue[which.min(pvalue)],")</td></tr><tr><td>#Tissues&nbsp;</td><td>",.N,"</td><td width='33%'>&nbsp;</td></tr><tr><td>Source(s)&nbsp;</td><td>",paste0(unique(source) %>% paste0(collapse = "|") %>% str_split("\\|") %>% unlist %>% unique ,collapse="<br>"),"</td></tr><tr><th class='tissue' colspan='3'>Tissue</th></tr><tr><td class='chead'>SNP&nbsp;</td><td class='chead' width='33%'>P-value</td><td class='chead' width='33%'>Source</td></tr></tbody></table>",paste0(label_tmp,collapse = "")), pvalue=min(pvalue), source = paste0(unique(source) %>% paste0(collapse = "|") %>% str_split("\\|") %>% unlist %>% unique ,collapse="|")),by=list(gene,gene_snp)]
    
    file.gtex.all$relationship <- "COREGULATION"
  
    file.gtex.all <- file.gtex.all %>%
      select(gene_snp, gene, pvalue, label, source, relationship) 
    
    ## WRITE
    write.table(file.gtex.all,file = paste0(file.output,"coregulation/significant_eQTL_all.csv"),sep = ",",col.names = T,row.names = F,quote = T)
  
    rm(file.gtex.all)
  }
  
  ## free memory
  gc()
  
  ## 3. GENE -> METABO --------------------------------------------------
  
  file.output <- data.makepath("data/CKDatlas/tmp/")

  if(!file.exists(paste0(file.output,"association_metabo_gene.csv"))){
  
    library(doParallel)
    library(foreach)
    library(reticulate)
    detectCores()
    
    myCluster <- makeCluster(10, # number of cores to use
                             type = "FORK") # type of cluster
    
    clusterEvalQ(myCluster, {
      #uri = "neo4j://metabo-interactive:7687"
      uri = "bolt://10.0.0.114:7687" # address of neo4j instance
      neo4j <- import(module = "neo4j",as = "neo4j")
      token <- neo4j$basic_auth("neo4j","sysdiab")
      driver <- neo4j$GraphDatabase$driver(uri, auth=token)
    })
    
    registerDoParallel(myCluster)
    
    genes <- "MATCH (g:proteinCoding) where exists((g)-[:MAP]->(:SNP)) RETURN g.ensembl_id as gene" %>% runQuery(read = TRUE) %>% .$gene
    
    total <-c()
    r<-c()
    counter=0
    for(i in split(1:length(genes), 1:length(genes)%/%70)){
      genes_test <- genes[i]
      counter=counter+1
      r <- foreach(gene = genes_test, .combine='rbind', .noexport = "driver") %dopar% { 
        
        "MATCH (g:proteinCoding) where g.ensembl_id=$properties MATCH (g)-[:MAP]->(s:SNP)-[:HIT]->(q:mQTL) where q.pvalue <= g.cutoff RETURN g.ensembl_id as gene, s.sid as snp, q.pvalue as pvalue, q.fluid as sample_type, [(q)-[:HIT]-(m:metMeta)|m.sid[0]][0] as metabolite, [(q)-[:FROM]-(so:source) | so.sid][0] as source, g.cutoff as cutoff" %>% runQuery(read = TRUE, param = TRUE, param_list = gene)
      }
      fwrite(r,file = paste0("CKDatlas/tmp/association_metabo_gene_",counter,".csv"))
      
    }
    
    stopCluster(myCluster)
    
    query <- paste0("CALL apoc.export.csv.query(\" MATCH (g:proteinCoding) where not exists(g.cutoff)
                     MATCH (g)-[:MAP]-(s:SNP)-[:HIT]-(q:mQTL) where q.pvalue <= 0.00000005 
                     RETURN g.ensembl_id as gene,s.sid as snp,q.pvalue as pvalue,q.fluid as sample_type, [(q)-[:HIT]-(m:metMeta)|m.sid[0]][0] as metabolite, [(q)-[:FROM]-(so:source) | so.sid][0] as source, '0.00000005' as cutoff\",\"CKDatlas/tmp/association_metabo_gene_genome_wide_tmp.csv\", {})")
    runQuery(query)
    
    ## READ tmp files and concatinate
    mQTL = lapply(getFiles(dir=paste0(file.output), pat="association_metabo_*"), fread, header=T)
    mQTL = do.call(rbind,mQTL)
    fwrite(mQTL,file = paste0("CKDatlas/tmp/association_metabo_gene.csv"))
    file.remove(getFiles(dir=paste0(file.output), pat="association_metabo_gene_genome_wide_tmp.csv"))
    file.remove(getFiles(dir=paste0(file.output), pat="association_metabo_gene_[0-9]*.csv"))
  }
  
  if(!file.exists(paste0(file.output,"../genetic_association/association_metabo_gene.csv"))){
    #all <- as.data.frame(mQTL)
    mQTL <- fread(file = "CKDatlas/tmp/association_metabo_gene.csv") 
    ## SUMMARIZE -> add html edge code
    all <- mQTL[,.(label_tmp=paste0("<tr><th class='sub' colspan='2'>",unique(snp),"</th></tr>",
                                    paste0("<td>",source[order(pvalue)],"</td><td width='50%'>",signif(pvalue[order(pvalue)],3),"</td></tr>",collapse = " "),collapse = " "),
                   cutoff = unique(cutoff),
                   pvalue = min(pvalue),
                   sample_type= list(sample_type)),by=list(gene,metabolite,snp)]
      
    all <- all[,.(label=paste0("<table class='snipa-plot-tooltip'><thead><tr><th colspan='2' id='gene_assoc'>GENETIC_ASSOCIATION</th></tr></thead>
                               <tbody><tr><td>Min. p-value&nbsp;</td><td width='50%'>",signif(min(pvalue),3),"<br>(",snp[which.min(pvalue)],")</td></tr>
                               <tr><td>#SNPs&nbsp;</td><td>",.N,"</td></tr>
                               <tr><td>Sign. cutoff&nbsp;</td><td>", signif(unique(cutoff),3),"</td></tr>
                               <tr><th class='tissue' colspan='2'>SNP</th></tr>
                               <tr><td class='chead'>Study&nbsp;</td><td class='chead' width='50%'>P-value</td></tr>
                               </tbody></table>
                               <table class='snipa-plot-tooltip'>",
                               paste0(unique(label_tmp),collapse = " "),"</table>",collapse = " "),
                  pvalue = min(pvalue),
                  sample_type=paste0(unique(sample_type %>% unlist), collapse ="|")), by=list(gene,metabolite)]
    
    all$relationship <- "GENETIC_ASSOCIATION"
    
    all <- all %>%
      select(gene, metabolite, pvalue, label, relationship, sample_type) %>% 
      mutate(label = str_replace_all(label,"\n","")) %>% 
      mutate(label = str_replace_all(label," {1,}"," ")) #%>% filter(!metabolite=="") # 
    
    ## WRITE
    write.table(all,file = paste0(file.output,"../genetic_association/association_metabo_gene.csv"),sep = ",",col.names = T,row.names = F,quote = T)
    
    rm(all)
    rm(mQTL)
  }
  

  ## 4. TRAIT -> GENE ----------------------------------------
  
  file.output <- data.makepath("data/CKDatlas/tmp/")
  
  if(!file.exists(paste0(file.output,"../genetic_association/association_gene_trait.csv"))){

    ## get all traits and loop through 
    traits <- "MATCH (t:trait) where (t)-[:HIT]-(:traitQTL) return distinct t.sid[0] as trait" %>% runQuery(read = TRUE) %>% .$trait
   
    for(trait in traits){
      
      "CALL apoc.export.csv.query(\"
          MATCH (g:proteinCoding) where exists(g.cutoff)
          MATCH (g)-[:MAP]->(s:SNP)-[:HIT]->(q:traitQTL)-[:HIT]->(:trait {sid:[$properties]}) 
          where q.pvalue <= g.cutoff 
          RETURN g.ensembl_id as gene, s.sid as snp,q.pvalue as pvalue,$properties as trait, [(q)-[:FROM]-(so:source) | so.cohort][0] as cohort, [(q)-[:FROM]-(sor:source) | sor.publication][0] as publication, g.cutoff as cutoff\",\"CKDatlas/tmp/association_gene_trait_\"+apoc.text.replace(toLower($properties), '[^A-Za-z0-9]', '')+\".csv\", {params: {properties:$properties, batchSize: 10}})" %>% runQuery(param = TRUE, param_list = trait)
      
      "CALL apoc.export.csv.query(\"
          MATCH (g:proteinCoding) where not exists(g.cutoff) 
          MATCH (g)-[:MAP]->(s:SNP)-[:HIT]->(q:traitQTL)-[:HIT]->(:trait {sid:[$properties]}) 
          where q.pvalue <= 0.00000005 
          RETURN g.ensembl_id as gene, s.sid as snp,q.pvalue as pvalue,$properties as trait, [(q)-[:FROM]-(so:source) | so.cohort][0] as cohort, [(q)-[:FROM]-(sor:source) | sor.publication][0] as publication, '0.00000005' as cutoff\",\"CKDatlas/tmp/association_gene_trait_\"+apoc.text.replace(toLower($properties), '[^A-Za-z0-9]', '')+\"_genome_wide.csv\", {params: {properties:$properties, batchSize: 10}})" %>% runQuery(param = TRUE, param_list = trait)
   
      message(trait)
    }
    
    df <- data.frame(matrix(ncol = 7,dimnames = list(c(),c("gene","snp","pvalue","trait","cohort","publication","cutoff"))))
    
    for (i in list.files(paste0(file.output),pattern="association_gene_trait_[a-z0-9]*_genome_wide.csv",full.names = T)){
      tmp_gene <- str_replace(i,"_genome_wide","")
      message(tmp_gene)
      tmp_gene <- fread(tmp_gene)
      tmp <- fread(i)
      names(tmp) <- c("gene","snp","pvalue","trait","cohort","publication","cutoff")
      names(tmp_gene) <- c("gene","snp","pvalue","trait","cohort","publication","cutoff")
      if(nrow(tmp)>0){
        df <- rbind(df,tmp)
      }
      if(nrow(tmp_gene)>0){
        df <- rbind(df,tmp_gene)
      }
    }
    df <- df[-1,] # first row is NA
    df <- mutate_all(df, list(~na_if(.,"")))
    
    df <- df %>% mutate(publication=map_chr(publication, function(x){ ifelse(is.na(x),NA,paste0("PMID: ",strsplit(x,"_")[[1]][3])) })) %>%
      mutate(source=map2_chr(publication,cohort,function(x,y) ifelse(is.na(x),y,x))) %>% data.table()
  
    ## make HTML label --------
    df <- df[,.(label_tmp=paste0("<tr><th class='sub' colspan='2'>",unique(snp),"</th></tr>", 
                                 paste0("<td>",source[order(pvalue)],"</td><td width='50%'>",signif(pvalue[order(pvalue)],3),
                                        "</td></tr>",collapse = " "),collapse = " "), 
                cutoff = unique(cutoff),
                pvalue = min(pvalue)), 
             by=list(gene,trait,snp)]
      
    df <- df[,.(label=paste0("<table class='snipa-plot-tooltip'><thead><tr><th colspan='2' id='gene_assoc'>GENETIC_ASSOCIATION</th></tr></thead>
                             <tbody><tr><td>Min p-value&nbsp;</td><td width='50%'>",signif(min(pvalue),3),"<br>(",snp[which.min(pvalue)],")</td></tr>
                             <tr><td>SNPs&nbsp;</td><td>",.N,"</td></tr>
                             <tr><td>Sign. cutoff&nbsp;</td><td>",signif(unique(cutoff),3),"</td></tr>
                             <tr><th class='tissue' colspan='2'>SNP</th></tr>
                             <tr><td class='chead'>Study&nbsp;</td><td class='chead' width='50%'>P-value</td></tr>
                             </tbody></table>
                             <table class='snipa-plot-tooltip'>",
                             paste0(unique(label_tmp),collapse = " "),"</table>",collapse = " "),pvalue = min(pvalue)), by=list(gene,trait)]
  
    df$relationship <- "GENETIC_ASSOCIATION"
    df <- df %>%
      select(gene, trait, pvalue, label, relationship) %>% 
      mutate(label = str_replace_all(label,"\n","")) %>% 
      mutate(label = str_replace_all(label," {1,}"," "))
    
    write.table(df,file = paste0(file.output,"../genetic_association/association_gene_trait.csv"),sep = ",",col.names = T,row.names = F,quote = T)
    
    rm(df)
    rm(tmp)
    rm(tmp_gene)
    file.remove(getFiles(dir=paste0(file.output), pat="association_gene_trait_[a-z0-9]*_genome_wide.csv"))
    file.remove(getFiles(dir=paste0(file.output), pat="association_gene_trait_[a-z0-9]*.csv"))
  }
  
  ## 5. TRAIT -> METABO #############################################################################################
  
  file.output <- data.makepath("data/CKDatlas/tmp/")

  if(!file.exists(paste0(file.output,"association_metabo_trait_metMeasured_00009090909"))){
    ## CYPHER TO CONVERT EDGE LIST TO CSV
    query <- paste0("CALL apoc.export.csv.query(\"MATCH (t:trait)<-[:HIT]-(w:mWAS)<-[:HIT]-(m:metMeta) where (w)-[:FROM]-(:source {metaboliteID:'metMeasured'})  and w.pvalue <= 0.0009090909 and w.sampleType='urine'
                   RETURN t.sid[0] as trait,w.pvalue as pvalue,m.sid[0] as metabolite, [(w)-[:FROM]-(so:source) | so.cohort][0] as source, w.sampleType as sample_type\",\"CKDatlas/tmp/association_metabo_trait_metMeasured_00009090909.csv\", {})")
    runQuery(query)
  }
  
  if(!file.exists(paste0(file.output,"association_metabo_trait_metMeasured_0000102459"))){
    ## CYPHER TO CONVERT EDGE LIST TO CSV
    query <- paste0("CALL apoc.export.csv.query(\"MATCH (t:trait)<-[:HIT]-(w:mWAS)<-[:HIT]-(m:metMeta) where (w)-[:FROM]-(:source {metaboliteID:'metMeasured'}) and w.pvalue <= 0.000102459 and w.sampleType='serum'
                   RETURN t.sid[0] as trait,w.pvalue as pvalue,m.sid[0] as metabolite, [(w)-[:FROM]-(so:source) | so.cohort][0] as source, w.sampleType as sample_type\",\"CKDatlas/tmp/association_metabo_trait_metMeasured_0000102459.csv\", {})")
    runQuery(query)
  }
  
  if(!file.exists(paste0(file.output,"association_metabo_trait_p180.csv"))){
    ## CYPHER TO CONVERT EDGE LIST TO CSV
    query <- paste0("CALL apoc.export.csv.query(\"MATCH (t:trait)<-[:HIT]-(w:mWAS)<-[:HIT]-(m:metMeta) where (w)-[:FROM]-(:source {metaboliteID:'p180'}) and w.pvalue <= 0.0003571429
                   RETURN t.sid[0] as trait,w.pvalue as pvalue,m.sid[0] as metabolite, [(w)-[:FROM]-(so:source) | so.cohort][0] as source, w.sampleType as sample_type\",\"CKDatlas/tmp/association_metabo_trait_p180.csv\", {})")
    runQuery(query)
  }
  
  if(!file.exists(paste0(file.output,"association_metabo_trait_p150.csv"))){
    ## CYPHER TO CONVERT EDGE LIST TO CSV
    query <- paste0("CALL apoc.export.csv.query(\"MATCH (t:trait)<-[:HIT]-(w:mWAS)<-[:HIT]-(m:metMeta) where (w)-[:FROM]-(:source {metaboliteID:'p150'}) and w.pvalue <= 0.0003311258
                   RETURN t.sid[0] as trait,w.pvalue as pvalue,m.sid[0] as metabolite, [(w)-[:FROM]-(so:source) | so.cohort][0] as source, w.sampleType as sample_type\",\"CKDatlas/tmp/association_metabo_trait_p150.csv\", {})")
    runQuery(query)
  }
  
  if(!file.exists(paste0(file.output,"../metabolic_association/association_metabo_trait.csv"))){
    
    ## READ files 
    #bile <- fread(paste0(file.output,"association_metabo_trait_bile.csv"))
    metMeasured1 <- fread(paste0(file.output,"association_metabo_trait_metMeasured_00009090909.csv"))
    metMeasured2 <- fread(paste0(file.output,"association_metabo_trait_metMeasured_0000102459.csv"))
    p180 <- fread(paste0(file.output,"association_metabo_trait_p180.csv"))
    p150 <- fread(paste0(file.output,"association_metabo_trait_p150.csv"))
    
    ## APPEND
    df <- rbind(metMeasured1,metMeasured2, p180, p150)
    #names(df) <- c("trait","pvalue","metabolite","source","cutoff")
    
    ## SUMMARIZE -> add html edge code
    df <- df %>% group_by(metabolite,trait) %>%
      dplyr::summarise(label=paste0("<table class='snipa-plot-tooltip'><thead><tr><th colspan='2' id='metabo_assoc'>METABOLIC_ASSOCIATION</th></tr></thead>
                             <tbody><tr><td>Min. p-value&nbsp;</td><td width='50%'>",signif(min(pvalue),3),"</td></tr>
                             <tr><td>#Studies&nbsp;</td><td>",n(),"</td></tr>
                             <tr><td>Sign. cutoff&nbsp;</td><td>study-dependen</td></tr>
                             <tr><td>Sample type(s)&nbsp;</td><td>",paste0(unique(sample_type),collapse="<br>"),"</td></tr>
                             <tr><th class='tissue' colspan='2'>-</th></tr>
                             <tr><td class='chead'>Study&nbsp;</td><td class='chead' width='50%'>P-value</td></tr>
                             </tbody></table>
                             <table class='snipa-plot-tooltip'><tr><th class='sub' colspan='2'>-</th></tr>",
                             paste0("<td>",source[order(pvalue)],"</td><td width='50%'>",signif(pvalue[order(pvalue)],3),"</td></tr>",collapse = " "),"</table>",collapse = " "),pvalue = min(pvalue), sample_type=paste0(unique(sample_type),collapse = "|"))
    
    df$relationship <- "METABOLIC_ASSOCIATION"
    df <- df %>%
      select(metabolite, trait, pvalue, label, relationship, sample_type) %>% 
      mutate(label = str_replace_all(label,"\n","")) %>% 
      mutate(label = str_replace_all(label," {1,}"," "))
    #df$id <- paste0("TQTLm_",c(1:nrow(df)))
    
    ## WRITE
    write.table(df,file = paste0(file.output,"../metabolic_association/association_metabo_trait.csv"),sep = ",",col.names = T,row.names = F,quote = T)
    
    #rm(bile)
    rm(metMeasured1)
    rm(metMeasured2)
    rm(df)
    rm(p180)
    rm(p150)
  }
  
  ## 6. GGM #############################################################################################
  
  file.output <- data.makepath("data/CKDatlas/tmp/")

  if(!file.exists(paste0(file.output,"ggm_tmp.csv"))){
    ## CYPHER TO CONVERT EDGE LIST TO CSV
    query <- paste0("CALL apoc.export.csv.query(\" MATCH (m1:metMeta)-[a:GGM]->(g:edgeGGM)<-[b:GGM]-(m2:metMeta) where g.significant='TRUE' and ID(m1)>ID(m2) and not a.metabo = b.metabo
                   RETURN m1.sid[0] as metabolite_a, m2.sid[0] as metabolite_b, g.pvalue as pvalue, [(g)-[:FROM]-(so:source) | so.cohort][0] as source, g.sampleType as sample_type, g.partialCorr as pCorr\",\"CKDatlas/tmp/ggm_tmp.csv\", {})")
    runQuery(query)
  }
  
  if(!file.exists(paste0(file.output,"../partial_correlation/ggm.csv"))){
    
    ## READ file
    df <- fread(paste0(file.output,"ggm_tmp.csv"))
   
    ## SUMMARIZE -> add html edge code
    df <- df[,.(label=paste0("<table class='snipa-plot-tooltip'><thead><tr><th colspan='2' id='ggm'>PARTIAL_CORRELATION</th></tr></thead>
                             <tbody><tr><td>Min. p-value&nbsp;</td><td width='50%'>", signif(min(pvalue),3),"</td></tr>
                             <tr><td>#Studies&nbsp;</td><td>",.N,"</td></tr>
                             <tr><td>Sample type(s)&nbsp;</td><td>",paste0(unique(sample_type),collapse="<br>"),"</td></tr>
                             </tbody></table>
                             <table class='snipa-plot-tooltip'><tr><th class='tissue' colspan='3'>-</th></tr>
                             <tr><td class='chead'>Study&nbsp;</td><td class='chead' width='33%'>P-value</td><td class='chead' width='33%'>pCorr.</td></tr>",
                             paste0("<td>",source[order(pvalue)],"</td><td width='33%'>",signif(pvalue[order(pvalue)],3),"</td><td width='33%'>",signif(pCorr[order(pvalue)],3),"</td></tr>",collapse = " "),
                             "</table>",collapse = " "),
                sample_type=paste0(unique(sample_type),collapse = "|"),
                pCorr=pCorr,
                pvalue = min(pvalue)),
             by=list(metabolite_a,metabolite_b)]
      
    df$relationship <- "PARTIAL_CORRELATION"
    
    df <- df %>%
      select(metabolite_a, metabolite_b, pvalue, label, relationship, sample_type, pCorr) %>% 
      mutate(label = str_replace_all(label,"\n","")) %>% 
      mutate(label = str_replace_all(label," {1,}"," "))
    
    ## WRITE
    write.table(df,file = paste0(file.output,"../partial_correlation/ggm.csv"),sep = ",",col.names = T,row.names = F,quote = T)
    
    rm(df)
  }
  
  ## 7. COEXPR #############################################################################################
  
  #file.output <- data.makepath("data/CKDatlas/tmp/")
  
  #if(!file.exists(paste0(file.output,"coexpr_tmp.csv"))){
    ## CYPHER TO CONVERT EDGE LIST TO CSV
  #  query <- paste0("CALL apoc.export.csv.query(\" MATCH (g1:proteinCoding)-[:COEXPR]->(c:edgeCOEXPR)<-[:COEXPR]-(g2:proteinCoding) where ID(g1)>ID(g2)
    #               RETURN g1.ensembl_id as gene_a, g2.ensembl_id as gene_b, c.pvalue as pvalue,c.sample_type as tissue, [(c)-[:FROM]-(so:source) | so.cohort][0] as source\",\"CKDatlas/tmp/coexpr_tmp.csv\", {})")
   # runQuery(query)
  #}
  
  #if(!file.exists(paste0(file.output,"../coexpression/coexpr.csv"))){
    ## READ file
   # df <- fread(paste0(file.output,"coexpr_tmp.csv"))
    
    ## SUMMARIZE -> add html edge code
   # df <- df[,.(label=paste0("<table class='snipa-plot-tooltip'><thead>
   #                          <tr><th colspan='2' id='coexpr'>COEXPRESSION</th></tr>
   #                          </thead><tbody>
    #                         <tr><td>Study&nbsp;</td><td>",paste0(unique(source),collapse="<br>"),"</td></tr>
     #                        <tr><td>Tissue&nbsp;</td><td>",paste0(unique(tissue),collapse="<br>"),"</td></tr>
     #                        </tbody></table>",collapse = " "), 
      #          tissue=paste0(unique(tissue),collapse = "|")), 
      #       by=list(gene_a,gene_b)]
    
    #df$relationship <- "COEXPRESSION"
    
    #df <- df %>%
     # select(gene_a, gene_b, label, relationship,tissue) %>% 
    #  mutate(label = str_replace_all(label,"\n","")) %>% 
    #  mutate(label = str_replace_all(label," {1,}"," "))
    
    ## WRITE
    #write.table(df,file = paste0(file.output,"../coexpression/coexpr.csv"),sep = ",",col.names = T,row.names = F,quote = T)
 # }
  
  ## 8. COABUNDANCE #############################################################################################
  
  file.output <- data.makepath("data/CKDatlas/tmp/")
 
  if(!file.exists(paste0(file.output,"ggm_prot_tmp.csv"))){
    ## CYPHER TO CONVERT EDGE LIST TO CSV
    query <- paste0("CALL apoc.export.csv.query(\"MATCH (g1:proteinCoding)-[:CODES]->(:transcript)-[:CODES]->(a:protein)-[c:COABUN]-(n:edgeCOABUN)-[d:COABUN]-(b:protein)<-[:CODES]-(:transcript)<-[:CODES]-(g2:proteinCoding) where ID(a)>ID(b) and not c.protein = d.protein 
                   RETURN distinct g1.ensembl_id as gene_a, g2.ensembl_id as gene_b, n.pvalue as pvalue,a.uniprot_id[0] as uniprot_a, b.uniprot_id[0] as uniprot_b, [(n)-[:FROM]-(so:source) | so.publication][0]as source, n.sample_type as sample_type, n.partialCorr as pCorr\",\"CKDatlas/tmp/ggm_prot_tmp.csv\", {})")
    runQuery(query)
  }

  if(!file.exists(paste0(file.output,"../coabundance/ggm_protein.csv"))){
    
    ## READ file
    df <- fread(paste0(file.output,"ggm_prot_tmp.csv")) %>% 
      mutate(protein=paste0(uniprot_a," - ",uniprot_b)) %>% 
      select(-uniprot_a,-uniprot_b) %>% data.table()
   
    ## SUMMARIZE -> add html edge code
    df <- df[,.(label_tmp=paste0("<tr><th class='sub' colspan='2'>",unique(protein),"</th></tr>",
                                 paste0("<td>",source[order(pvalue)],"</td><td width='33%'>",signif(pvalue[order(pvalue)],3),"</td><td width='33%'>",signif(pCorr[order(pvalue)],3),"</td></tr>",collapse = " ")
                                 ,collapse = " "),pvalue = min(pvalue), pCorr = pCorr[which.max(abs(pCorr))], sample_type = list(sample_type)),by=list(gene_a,gene_b,protein)]
      
    df <- df[,.(label=paste0("<table class='snipa-plot-tooltip'><thead><tr><th colspan='2' id='coabun'>COABUNDANCE</th></tr></thead>
                                    <tbody><tr><td>Min. p-value&nbsp;</td><td width='50%'>",signif(min(pvalue),3),"<br>(",protein[which.min(pvalue)],")</td></tr>
                                    <tr><td>#Protein pairs&nbsp;</td><td>",.N,"</td></tr>
                                    <tr><td>Sign. cutoff&nbsp;</td><td>Study-specific</td></tr>
                                    <tr><td>Sample type(s)&nbsp;</td><td>",paste0(unique(sample_type %>% unlist),collapse="<br>"),"</td></tr>
                                    <tr><th class='tissue' colspan='2'>Proteins</th></tr>
                                    <tr><td class='chead'>Study&nbsp;</td><td class='chead' width='33%'>P-value</td><td class='chead' width='33%'>pCorr.</td></tr>
                                    </tbody></table><table class='snipa-plot-tooltip'>",
                             paste0(unique(label_tmp),collapse = " "),"</table>",collapse = " "),pvalue = min(pvalue),pCorr = pCorr[which.max(abs(pCorr))],sample_type=paste0(unique(sample_type %>% unlist),collapse="|")),by=list(gene_a,gene_b)]
    
    df$relationship <- "COABUNDANCE"
    
    df <- df %>%
      select(gene_a,gene_b, pvalue, label, relationship, sample_type) %>% 
      mutate(label = str_replace_all(label,"\n","")) %>% 
      mutate(label = str_replace_all(label," {1,}"," "))
    
    ## WRITE
    write.table(df,file = paste0(file.output,"../coabundance/ggm_protein.csv"),sep = ",",col.names = T,row.names = F,quote = T)
  }
  
  info(log,capture.output(toc()))

}