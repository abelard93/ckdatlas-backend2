library(doParallel)
library(foreach)
library(reticulate)
detectCores()

myCluster <- makeCluster(10, # number of cores to use
                         type = "FORK") # type of cluster

clusterEvalQ(myCluster, {
  uri = "neo4j://metabo-interactive:7687"
  neo4j <- import(module = "neo4j",as = "neo4j")
  token <- neo4j$basic_auth("neo4j","sysdiab")
  driver <- neo4j$GraphDatabase$driver(uri, auth=token)
})

registerDoParallel(myCluster)

genes <- "MATCH (g:proteinCoding) where exists((g)-[:MAP]->(:SNP)) RETURN g.ensembl_id as gene" %>% readSimpleQuery() %>% .$gene

#genes <- head(genes, n=2)
 
total <-c()
r<-c()
counter=0
for(i in split(1:length(genes), 1:length(genes)%/%70)){
  genes_test <- genes[i]
  counter=counter+1
  r <- foreach(gene = genes_test, .combine='rbind', .noexport = "driver") %dopar% { 
    
    "MATCH (g:proteinCoding) where g.ensembl_id=$properties MATCH (g)-[:MAP]->(s:SNP)-[:HIT]->(q:mQTL) where q.pvalue <= g.cutoff RETURN g.ensembl_id as gene, s.sid as snp, q.pvalue as pvalue, [(q)-[:HIT]-(m:metMeta)|m.sid[0]][0] as metabolite, [(q)-[:FROM]-(so:source) | so.sid][0] as source, g.cutoff as cutoff" %>% readQueryWithParam(.,param = gene)
  }
  fwrite(r,file = paste0("ADatlas/misc/association_metabo_",counter,".csv"))
   
}

stopCluster(myCluster)





