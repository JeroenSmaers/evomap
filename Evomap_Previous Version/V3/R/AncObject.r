#' Function to sort data according to ie output (following the ape tree structure)
#'
#' This function sorts extant and nodal values according the ape tree structure
#' @param tree an object of class "phylo".
#' @param data_extants named vector of tip values
#' @param data_nodes named vector of ancestral states
#' @return object with extant and nodal values sorted according to the ape tree structure
#' @export
AncObject<-function(tree,data_extants,data_nodes){

#Create general structure based on tree in ape format
            matrix=tree$edge
            phy.matrix=data.frame(tree$edge,tree$edge.length,data_extants[matrix[,2]],data_extants[matrix[,2]])
            names(phy.matrix)=c("node_anc","node_desc","BL","value_anc","value_desc")
#Add ancestors to value_desc
                  for(i in 1:length(names(data_nodes))){
                        phy.matrix$value_desc[which(phy.matrix$node_desc==names(data_nodes[i]))]<-data_nodes[[i]]
                  }
#Fill in value_anc
                  for(i in 1:length(names(data_nodes))){
                        phy.matrix$value_anc[which(phy.matrix$node_anc==names(data_nodes[i]))]<-data_nodes[[i]]
                  }
#Add rate and change
                  phy.matrix$rate<-(phy.matrix$value_desc-phy.matrix$value_anc)/sqrt(phy.matrix$BL)
                  phy.matrix$change<-(phy.matrix$value_desc-phy.matrix$value_anc)
#branching times of ancestral nodes
      times<-branching.times(tree)
      time_node_anc<-c()
      for (i in phy.matrix$node_anc){
      time_node_anc<-append(time_node_anc,times[which(names(times)==i)])}
      phy.matrix$time_node_anc<-time_node_anc

return(phy.matrix)
}
