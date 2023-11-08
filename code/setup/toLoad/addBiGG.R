#################
## BiGG model
#################

addBiGG <- function(){
  # create index on Symbol for faster search in functions
  # query <- "CREATE INDEX ON :metBiGG(metHMDB)"
  # cypher(db, query)
  
  queryConstraint <- "RETURN apoc.schema.node.constraintExists($properties, ['sid']) as exists"
  
  
  for (file in getFiles("BiGG", pat = "*_genes.csv")) {
    file <- gsub("_genes.csv",".csv",file)
    # read info on dataset
    source <- getInfo(file, c("sid", "organism", "file_name", "model"), err)
    # add model and put uniqueness constraint on nodel type
    if(!runQuery(queryConstraint,read = TRUE, param = TRUE, param_list = source$model)$exists=="TRUE"){
      add_constr(source$model, "sid")
    }
    
    # skip file if no/incomplete info
    if (!all(is.na(source))) {
      tic(paste0("Added BiGGModel ", source["sid"]))
      # check if node for source is present, if not create
      checkSource(source, log)
      # check if node for organism is present, if not create
      #checkOrganism(source,db,log)
      # general path to file
      tmp_path <- tools::file_path_sans_ext(file)
      model <- source$model
      ## ADD BIGG METABOLITES
      query <- paste0("
                        USING PERIODIC COMMIT 1000
                        LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file+'_metabolites.csv' AS props
                        MERGE (n:metBiGG {sid:props.sid})
                        SET n += props set n:",model," with n
                        MATCH (s:source {sid:$properties.source})
                        MERGE (n)-[:FROM]->(s) with n 
                        where exists(n.compartment)
                        MERGE (c:compBiGG {sid:n.compartment, id:n.`compartment.id`})
                        MERGE (n)-[:CELL_COMPARTMENT]->(c)
                        ")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=tmp_path,source=source["sid"]$sid))
      ## MAP METABOLITES TO HMDB
      query <- paste0("
                        USING PERIODIC COMMIT 1000
                        LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file+'_metabolites.csv' AS props
                        MATCH (n:metBiGG:",model," {sid:props.sid}) WHERE exists(n.hmdb) with n
                        MATCH (h1:accHMDB)-[]-(h:metHMDB) where h1.sid=n.hmdb 
                        with n,h
                        MERGE (n)-[:MAP]->(h)
                        ")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=tmp_path))
      ## ADD BIGG REACTIONS
      query <- paste0("
                        USING PERIODIC COMMIT 1000
                        LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file+'_reactions.csv' AS props
                        MERGE (n:reBiGG {sid:props.sid})
                        SET n += props set n:",model," with n
                        MATCH (s:source {sid:$properties.source})
                        MERGE (n)-[:FROM]->(s)
                        ")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=tmp_path,source=source["sid"]$sid))
      
      
      ## ADD BIGG GENES
      query <- paste0("
                        USING PERIODIC COMMIT 1000
                        LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file+'_genes.csv' AS props
                        MERGE (n:geBiGG {sid:props.sid})
                        SET n += props set n:",model," with n
                        MATCH (s:source {sid:$properties.source})
                        MERGE (n)-[:FROM]->(s)
                        ")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=tmp_path,source=source["sid"]$sid))
      
      ## ADD BIGG GENES - ENSEMBLE relationship
      query <- paste0("
                        USING PERIODIC COMMIT 1000
                        LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file+'_genes.csv' AS props
                        MATCH (n:geBiGG:",model," {sid:props.sid}) with n 
                        where exists(n.ensembl_gene)
                        MERGE (g:ensembl:gene {sid:n.ensembl_gene, ensembl_id:n.ensembl_gene})
                        ON CREATE set g.source=$properties.source
                        with g,n
                        MERGE (g)-[m:MAP]->(n) 
                        ON CREATE set m.source=$properties.source
                        ")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=tmp_path,source=paste0("'BiGG_",model,"'")))
      
      # query <- paste0("
      #                   USING PERIODIC COMMIT 1000
      #                   LOAD CSV WITH HEADERS FROM 'file:///",tmp_path,"_genes.csv' AS props
      #                   MATCH (n:geBiGG:",model," {sid:props.sid}) with n 
      #                   where exists(n.entrez_id)
      #                   MERGE (g:entrez:gene {sid:n.entrez_id, entrez_id:n.entrez_id})
      #                   ON CREATE set g.source=",paste0("'BiGG_",model,"'"),"
      #                   with g,n
      #                   MERGE (g)-[m:MAP]->(n) 
      #                   ON CREATE set m.source=",paste0("'BiGG_",model,"'"),"
      #                   ")
      # cypher(db, query)
      
      ## ADD "CATALYZES" REALTIONSHIP
      query <- paste0("
                        USING PERIODIC COMMIT 1000
                        LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file+'_re_gene.csv' AS props 
                        MATCH (g:geBiGG:",model,") where g.sid=props.gene
                        MATCH (r:reBiGG:",model,") where r.sid=props.reBiGG
                        MERGE (g)-[:CATALYZES {source:$properties.source}]->(r) 
                        ")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=tmp_path,source=source["sid"]$sid))
      ## ADD DIRECTED "PARTICIPATES" RELATIONSHIP
      # positiv stochiometry: reaction --> metabolite
      query <- paste0("
                        USING PERIODIC COMMIT 1000 
                        LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file+'_met_re.csv' AS props 
                        with props where toFloat(props.stoichiometry) >= 0
                        MATCH (m:metBiGG:",model," {sid:props.metBiGG}) 
                        MATCH (r:reBiGG:",model," {sid:props.reBiGG})
                        MERGE (r)-[:PARTICIPATES {source:$properties.source, stoichiometry:props.stoichiometry}]->(m)
                        ")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=tmp_path,source=source["sid"]$sid))
      # negative stochiometry: metabolite --> reaction
      query <- paste0("
                        USING PERIODIC COMMIT 1000 
                        LOAD CSV WITH HEADERS FROM 'file:///'+$properties.file+'_met_re.csv' AS props 
                        with props where toFloat(props.stoichiometry) < 0
                        MATCH (m:metBiGG:",model," {sid:props.metBiGG}) 
                        MATCH (r:reBiGG:",model," {sid:props.reBiGG})
                        MERGE (m)-[:PARTICIPATES {source:$properties.source, stoichiometry:props.stoichiometry}]->(r)
                        ")
      runQuery(query, param = TRUE, periodic = TRUE, param_list = list(file=tmp_path,source=source["sid"]$sid))
      
      info(log, capture.output(toc()))
    }
  }
}