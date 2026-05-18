#' List all edges descendant from a node
#'
#' @param tree an object of class 'phylo'.
#' @param node number of interest.
#' @return vector of edge numbers.

#' @export

getEdges<-function(tree,node){
      Desc<-getDescendants(tree,node)
      Edges<-c()
      for(i in 1:length(Desc)){Edges<-c(Edges,which(tree$edge[,2]==Desc[i]))}
      Edges<-c(Edges,which(tree$edge[,2]==node))
      return(Edges)}
      
