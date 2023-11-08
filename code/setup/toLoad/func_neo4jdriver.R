## HELPER FUNCTIONS FOR NEO4J COMMUNICATION
## reticulate package enables use of python driver
library(reticulate)
library(data.table)

#lapply(function(x) lapply(x,function(y) paste0(y,collapse = "|")))

runQuery <- function(query, read=FALSE, param = FALSE, periodic = FALSE, network=FALSE, param_list=NULL){
  
  ## needed functions
  
  submitQuery <- function(tx,query){tx$run(query)}
  submitReadQuery <- function(tx,query){
    result <- tx$run(query) 
    result$data()}
  submitReadQueryWithParam <- function(tx,query,param){
    result <- tx$run(query, properties=param) 
    result$data()}
  submitQueryWithParam <- function(tx,query,param){tx$run(query, properties=param)}
 
  case <- c(ifelse(read,1,0),ifelse(param,1,0),ifelse(periodic,1,0)) %>% paste0(collapse = "")
  
  # start session
  session <- driver$session()
  
  if(case == "000"){ # read=F, param=F, periodic=F
    session$write_transaction(submitQuery,query) 
    session$close()
  }
   
  if(case == "100"){ # read=T, param=F, periodic=F
    result <- session$read_transaction(submitReadQuery,query) 
    if(!network){
      result <- result %>% rapply(., function(x) ifelse(all(x=="NA"),NA,paste0(x,collapse = "|")), how="replace") %>% rbindlist(fill = TRUE, idcol = F)
    }
    session$close()
    return(result) 
  }
  
  if(case == "010"){ # read=F, param=T, periodic=F
    session$write_transaction(submitQueryWithParam,query,param_list) 
    session$close()
  }
  
  if(case == "110"){ # read=F, param=T, periodic=T
    result <- session$read_transaction(submitReadQueryWithParam,query,param_list) 
    if(!network){
      result <- result %>% rapply(., function(x) ifelse(all(x=="NA"),NA,paste0(x,collapse = "|")), how="replace") %>% rbindlist(fill = TRUE, idcol = F)
    }
    session$close()
    return(result) 
  }
  
  if(case == "011"){ # read=F, param=T, periodic=T
    session$run(query,properties=param_list)
    session$close()
  }
  
  if(case == "001"){ # read=F, param=F, periodic=T
    session$run(query)
    session$close()
  }
  
  
}




