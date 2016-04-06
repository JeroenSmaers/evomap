#' Function for the method of Independent Evolution
#'
#' This function allows computing rescaled branch lengths according to inferred phenotypic evolution and ancestral states using the method of Independent Evolution as described in Smaers & Vinicius (2009)
#' @param data numeric vector with elements in the same order as tiplabels in the tree
#' @return dataframe with rescaled branch lengths (rBL) for all branches in the tree
#' @details This function poses the following restrictions on the input data: no polytomies in the tree, no duplicated values, no zero values, and no negative values in the data. The algorithm was designed to deal with linear data (unlogged). This function is now deprecated; 'mvBM' should be used instead. 
#' @export

ie<-function(data,tree){
data_original<-data
            N=length(data)
      #Make phylo matrix and list nodes
            matrix=tree$edge
            phy.matrix=data.frame(tree$edge,tree$edge.length,data[matrix[,2]])
            names(phy.matrix)=c("Anc","Desc","Length","Value")
            phy.matrix.length=nrow(phy.matrix)
            nodes_extant=1:N
            nodes_extinct=(N+1):(N+(N-1))
            nodes_all=1:(N+(N-1))


#CODE

#1. CALCULATE AP BRANCH LENGTHS
            matrix_lengths=dist.nodes(tree)
#2. CALCULATE AP-values

      values=phy.matrix$Value[which(phy.matrix$Desc<=N)]
      extant_values=data.frame(values,nodes_extant)

      AP_values=c()

      for(j in nodes_extinct){
            #calculate nominator
                  nominator=c()
                        for(i in nodes_extant){
                              nominator=rbind(nominator,(extant_values[i,1]/matrix_lengths[i,j]))
                                              }
            #calculate denominator
                  denominator=c()
                        for(i in nodes_extant){
                              denominator=rbind(denominator,(1/matrix_lengths[i,j]))
                                              }
            #save AP
                  AP=sum(nominator)/sum(denominator)
                  AP_values=rbind(AP_values,AP)
                             }

      AP=c()
      AP=cbind(AP,AP_values)
      AP=cbind(AP,nodes_extinct)
      AP=as.data.frame(AP)
      colnames(AP)=c("value","nodes")

#3. CALCULATE ANCESTRAL STATES

nodes_extinct_reverse=sort(nodes_extinct,decreasing=TRUE)
            ancestral_states=c()
            ancestors=c()
            results=c()
            results_rBLs=c()
            results_branchlength=c()
            results_node_anc=c()
            results_node_desc=c()

for(i in nodes_extinct_reverse) {
        #calculating ancestral states
            sister_species=which(phy.matrix$Anc==i)

            X1=phy.matrix$Value[sister_species[1]]
            X2=phy.matrix$Value[sister_species[2]]
            AP_desc=AP$value[which(AP$nodes==i)]

            BL1=phy.matrix$Length[sister_species[1]]
            BL2=phy.matrix$Length[sister_species[2]]


            S1=abs(abs(X1-X2)/mean(c(X1,X2)))
            S2=abs(abs(X2-AP_desc)/mean(c(X2,AP_desc)))
            S3=abs(abs(X1-AP_desc)/mean(c(X1,AP_desc)))
            
            T1=((S1+S3)-S2)/2
            T2=((S1+S2)-S3)/2
            T3=((S2+S3)-S1)/2

            rBL1=T1*((BL1/(BL1+BL2))*2)
            rBL2=T2*((BL2/(BL1+BL2))*2)



            A=((X1/rBL1)+(X2/rBL2))/((1/rBL1)+(1/rBL2))

        #Saving results

            desc1=phy.matrix$Desc[sister_species[1]]
            desc2=phy.matrix$Desc[sister_species[2]]

            value_desc1=phy.matrix$Value[which(phy.matrix$Desc==desc1)]
            value_desc2=phy.matrix$Value[which(phy.matrix$Desc==desc2)]

            phy.matrix$Value[which(phy.matrix$Desc==i)]=A

        #building results dataframe
            results_node_anc=rbind(results_node_anc,i)
            results_node_anc=rbind(results_node_anc,i)
            results_node_desc=rbind(results_node_desc,desc1)
            results_node_desc=rbind(results_node_desc,desc2)
            results_rBLs=rbind(results_rBLs,rBL1)
            results_rBLs=rbind(results_rBLs,rBL2)
            results_branchlength=rbind(results_branchlength,BL1)
            results_branchlength=rbind(results_branchlength,BL2)
                                }

      #'results' dataframe
      results<-c()
      results=cbind(results,results_node_anc)
      results=cbind(results,results_node_desc)
      results=cbind(results,results_branchlength)
      results=cbind(results,results_rBLs)
      colnames(results)=c("node_anc","node_desc","BL","rBLs")

      #ancestors
      results<-results[order(match(results[,2],phy.matrix[,2])),]
      rownames(results)<-c(1:length(results[,1]))
      results<-as.data.frame(results)

      return(results)

}
