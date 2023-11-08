

makeSchema <- function(graph,key,date,git){
  
  library(purrr)
  library(visNetwork)
  
  
  edge_query <- "
  MATCH (n)-[r]->(m)
  RETURN distinct
  labels(n)[0] as from
  , labels(m)[0] as to
  , type(r) as label
  "
  
  label_query <- "
  MATCH (n) RETURN distinct labels(n) as labels
  "
  

  label <- cypher(graph, label_query)
  
  
  #unique(unlist(label$labels))
  
  m <- data.frame(matrix(
    nrow=length(unique(unlist(label$labels))),
    ncol=2,
    dimnames= list(row=NULL,col=c("label","all")),data = "")
    )
  
  m[,1] <- unique(unlist(label$labels))
  
  for(i in m$label){
    all <- unique(unlist(label$labels[sapply(label$labels, function(x) any(str_detect(x,paste0("^",i,"$"))))]))
    m$all[m$label==i] <- paste(all,collapse ="\n")
  }
  
  
  
  
  node_query <- "CALL apoc.meta.data() YIELD label, property, type, elementType
  WHERE not elementType='relationship'
  AND not type='RELATIONSHIP'
  WITH label, collect(property + ' (' + type + ')') AS properties
  WITH apoc.map.fromLists(collect(label), collect(properties)) AS properties_map, label
  RETURN apoc.create.vNode([''],{name:label,properties:apoc.text.join(properties_map[label], '<br />')}) as node"
  
  # get neo4j data
  nodes <- cypherToList(graph, node_query)
  
  nodes <- data.table::rbindlist(lapply(nodes,data.table::rbindlist)) %>% 
    data.table::setnames(.,c("id","title")) #%>% mutate(group=id)
  nodes$title <- sapply(nodes$title, function(x) paste0("<font size=\"2\" color=\"black\">",x,"</font>"))
  nodes$font.size <-  11
  
  edges <- cypher(graph, edge_query) %>% 
    mutate(
      from=map_chr(from,function(x) m$all[m$label==x]),
      to=map_chr(to,function(x) m$all[m$label==x]))
  
  edges$font.size <- 9
  edges$dashes <- ifelse(edges$label=="FROM",TRUE,FALSE)
  
  
  nodes <- nodes %>% mutate(id=map_chr(id,function(x) m$all[m$label==x]))
  nodes <- nodes[nodes$id %in% edges$from | nodes$id %in% edges$to,]
  # calculate the graph
  network <-
    visNetwork(nodes, edges, main = "Data model",submain = paste0("Hash: ", key,"<br />Date: ",date), height = "90vh", width="100%") %>%
    visEdges(arrows = "to",color = list(color = "grey", highlight = "red"),width=0.5) %>%
    visOptions(highlightNearest = TRUE) %>%
    visNodes(shape = "box") %>%
    #visHierarchicalLayout(sortMethod="directed") %>%
    #visIgraphLayout(type = "full") %>% 
    visInteraction(tooltipStyle = 'position: fixed;visibility:hidden;padding: 5px;background-color: white;overflow:scroll;outline-style:solid;width: 250px; height: 150px;') %>%
    visLayout(randomSeed = 1233) %>%visPhysics(solver = "forceAtlas2Based",
                                               forceAtlas2Based = list(gravitationalConstant = -20))
  
  
  # plot the graph
  network
  
  # save as HTML
  visSave(network, file = paste0(git,"schema/",format(Sys.time(), "%Y%m%d"),"_dbSchema.html"))
  
}