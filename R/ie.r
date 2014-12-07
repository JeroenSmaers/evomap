#' Function for the method of Independent Evolution
#'
#' This function allows computing rates and ancestral states using the method of Independent Evolution
#' @param data numeric vector with elements in the same order as tiplabels in the tree
#' @return dataframe with rates for all branches and ancestral states for all nodes in the tree
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

      AP<-c()
      AP<-cbind(AP,AP_values)
      AP<-cbind(AP,nodes_extinct)
      colnames(AP)<-c("value","nodes")
      rownames(AP)<-1:length(nodes_extinct)
      AP<-as.data.frame(AP)

      AP$value<-exp(AP$value)
      phy.matrix$Value[which(phy.matrix$Desc<=N)]<-exp(phy.matrix$Value[which(phy.matrix$Desc<=N)])
#3. CALCULATE ANCESTRAL STATES

nodes_extinct_reverse=sort(nodes_extinct,decreasing=TRUE)
            ancestral_states=c()
            ancestors=c()
            results=c()
            results_rates_abs=c()
            results_rates_rel=c()
            results_branchlength=c()
            results_value_anc=c()
            results_value_desc=c()
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

            R1=T1*((BL1/(BL1+BL2))*2)
            R2=T2*((BL2/(BL1+BL2))*2)

            A=exp(((log(X1)/R1)+(log(X2)/R2))/((1/R1)+(1/R2)))

        #Saving results
            R1_rel=(T1*((BL1/(BL1+BL2))*2))/BL1
            R2_rel=(T2*((BL2/(BL1+BL2))*2))/BL2

            desc1=phy.matrix$Desc[sister_species[1]]
            desc2=phy.matrix$Desc[sister_species[2]]

            value_desc1=phy.matrix$Value[which(phy.matrix$Desc==desc1)]
            value_desc2=phy.matrix$Value[which(phy.matrix$Desc==desc2)]

            rate1_abs=R1
            rate2_abs=R2

            rate1_rel=R1_rel
            rate2_rel=R2_rel

            phy.matrix$Value[which(phy.matrix$Desc==i)]=A

        #building results dataframe
            ancestral_states=rbind(ancestral_states,A)
            results_rates_abs=rbind(results_rates_abs,rate1_abs)
            results_rates_abs=rbind(results_rates_abs,rate2_abs)
            results_rates_rel=rbind(results_rates_rel,rate1_rel)
            results_rates_rel=rbind(results_rates_rel,rate2_rel)
            results_branchlength=rbind(results_branchlength,BL1)
            results_branchlength=rbind(results_branchlength,BL2)
            results_value_anc=rbind(results_value_anc,A)
            results_value_anc=rbind(results_value_anc,A)
            results_value_desc=rbind(results_value_desc,value_desc1)
            results_value_desc=rbind(results_value_desc,value_desc2)
            results_node_anc=rbind(results_node_anc,i)
            results_node_anc=rbind(results_node_anc,i)
            results_node_desc=rbind(results_node_desc,desc1)
            results_node_desc=rbind(results_node_desc,desc2)
                                }

      #'ancestors' dataframe
      ancestors=cbind(ancestors,ancestral_states)
      ancestors=cbind(ancestors,nodes_extinct_reverse)
      colnames(ancestors)=c("ancestor","node")
      rownames(ancestors)=paste(1:(N-1))

      #branching times of ancestral nodes
      times<-branching.times(tree)
      time_node_anc<-c()
      for (i in results_node_anc){
      time_node_anc<-append(time_node_anc,times[which(names(times)==i)])}

      #'results' dataframe
      results<-c()
      results=cbind(results,results_node_anc)
      results=cbind(results,results_node_desc)
      results=cbind(results,results_branchlength)
      results=cbind(results,time_node_anc)
      results=cbind(results,results_value_anc)
      results=cbind(results,results_value_desc)
      results=cbind(results,results_rates_abs)
      colnames(results)=c("node_anc","node_desc","BL","time_node_anc","value_anc","value_desc","phi")

      #ancestors
      results<-results[order(match(results[,2],phy.matrix[,2])),]
      rownames(results)<-c(1:length(results[,1]))
      results<-as.data.frame(results)

      results$value_desc<-log(results$value_desc)
      results$value_anc<-log(results$value_anc)

      return(results)

}
