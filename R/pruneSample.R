#' Pruning data and a tree to a particular subset
#'
#' Prunes a comparative dataset down to a particular clade
#' @param data a vector of tip values for species; should be in the same order as tiplabels in the tree
#' @param tree an object of class "phylo".
#' @param targetGroup vector listing the tip labels to which the data should be pruned
#' @return pruned comparative data set

#' @export
pruneSample<-function(data,tree,targetGroup){
    SpeciesIn<-tree$tip.label[targetGroup]; SpeciesOut<-setdiff(tree$tip.label,SpeciesIn)
    treeTemp<-drop.tip(tree,SpeciesOut); dataTemp<-as.data.frame(treedata(treeTemp,data,sort=T,warnings=F)$data)
    output<-list(treeTemp,dataTemp); names(output)<-c("treePruned","dataPruned")
    return(output)}
