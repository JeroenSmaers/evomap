#' List all nodes descendant from a node
#'
#' @param tree an object of class 'phylo'.
#' @param node number of interest.
#' @return vector of node numbers.

#' @export

getNodes<-function(tree,node){
c(getDescendants(tree,node)[which(getDescendants(tree,node)>(length(tree$tip.label)+1))],node)
                   }
