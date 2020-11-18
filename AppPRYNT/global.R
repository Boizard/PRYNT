#######Packages#######
####Packages####
library(igraph)
library(zoo)
library(RandomWalkRestartMH)
library(STRINGdb)
library(shinycssloaders)


#######Functions#######


shortest_path_ranking_function<-function(graph,seed,score=NULL,directed=T,distance_table=NULL){
  #calculate the closeness of each protein from seeds using the shortest path
  #graph : igraph object, containing the network to analyse
  #seed : character vector of DE proteins
  #score : NULL if we don't take weight in to account, or numeric vector containing weight of each relations (more reliable relations have a minimal score )
  #directed : boolean
  #distance_table : to have the possibility to calculate the distance table outside of the function
  
  if(is.null(distance_table)){
    if(directed){sens="out"}
    else{sens="all"}
    distance_table<-distances(graph,v = V(graph),to = V(graph), mode = sens,weights = score)
  }
  
  distance_table_to_seed<-distance_table[,colnames(distance_table)%in%seed]
  #0 means that the protein is DE. I put the distance of a DE to its self to 1.
  distance_table_to_seed[which(distance_table_to_seed==0,arr.ind = T)]<-Inf
  
  results_SP<-apply(X = distance_table_to_seed,MARGIN = 1,function(x){sum(1/x)})
  #resultat_closeness vector of adapted closeness centrality
  results_SP<-data.frame("node"=names(results_SP),"SP"=results_SP)
  results_SP<-results_SP[order(results_SP$node),]
  return(results_SP)
}

random_walk_ranking_function<-function(graph,seed,adjacency_matrix=NULL,multiplex_object=NULL){
  #calculate the closeness of each protein from seeds using the random walk 
  # using the RandomWalkRestartMH package
  #graph : igraph object, containing the network to analyse
  #seed : character vector of DE proteins
  
  #adjacency_matrix : to have the possibility to calculate the adjency matrix outside of the function
  #multiplex_object : to have the possibility to calculate the multiplex object outside of the function
  
  
  if(is.null(adjacency_matrix)){
    multiplex_object <- create.multiplex(graph,Layers_Name=c("PPI"))
  }
  if(is.null(adjacency_matrix)){
    adjacency_matrix <- compute.adjacency.matrix(multiplex_object)
  }
  adjacency_matrix_norm <- normalize.multiplex.adjacency(adjacency_matrix)
  seed<-seed[seed%in%names(V(graph))]
  ## We launch the algorithm with the default parameters (See details on manual)
  RWR_PPI_Results <- Random.Walk.Restart.Multiplex(adjacency_matrix_norm,
                                                   multiplex_object,seed)
  # We display the results
  results_random_walk<-RWR_PPI_Results$RWRM_Results
  colnames(results_random_walk)<-c("node","RW")
  node<-data.frame("node"=names(V(graph)))
  results_random_walk<-merge(x = node,y =results_random_walk,by="node",all.x = T )
  results_random_walk$RW[is.na(results_random_walk$RW)]<-0
  results_random_walk<-results_random_walk[order(results_random_walk$node),]
  
  return(results_random_walk)
}





########download the PPI network#####
###from string package
print("loading of string Database")
string_db <<- STRINGdb$new(version="11", species=9606,
                          score_threshold=0, input_directory="" )

annotation<<-string_db$get_proteins()


string_network<<-string_db$get_interactions(string_ids =annotation[,1] )
string_network_900<<-string_network[string_network$combined_score>=900,]
proteins_900<<-unique(c(string_network_900$from,string_network_900$to))