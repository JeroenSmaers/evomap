#' List all tips descendant from a node
#'
#' @param tree an object of class 'phylo'.
#' @param node number of interest.
#' @return vector of tip numbers.

#' @export

getTips<-function(tree,node){
getDescendants(tree,node)[which(getDescendants(tree,node)<(length(tree$tip.label)+1))]
                   }
