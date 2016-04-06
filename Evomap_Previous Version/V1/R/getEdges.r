getEdges<-function(tree,node){
      Desc<-getDescendants(tree,node)
      Edges<-c()
      for(i in 1:length(Desc)){Edges<-c(Edges,which(tree$edge[,2]==Desc[i]))}
      Edges<-c(Edges,which(tree$edge[,2]==node))
      return(Edges)}
      
