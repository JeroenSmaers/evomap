#' Function for mapping ie results onto a phylogenetic tree
#'
#' This function maps the results of an ie analysis onto a tree
#' @param ieObject object comprising the output of an ie analysis
#' @return phylogenetic tree with mapped change and ancestral states
#' @export

iePlot<-function(ieObject,tree,branchWidth,branchCutoff,labelSize,scaling){

                                    par(bg="black")
            #Branches
                  cutoff<-(sd(ieObject$change)*branchCutoff)
                  branch_widths_ieObject<-sqrt(ieObject$change*ieObject$change)*branchWidth
                  branch_col_ieObject<-rep("white",length(ieObject$change))
                  for(i in 1:length(branch_col_ieObject)){
                  if(ieObject$change[i]>cutoff)branch_col_ieObject[i]<-"green" else (if(ieObject$change[i]<=-cutoff)branch_col_ieObject[i]<-"red")}
                        plot(tree,edge.col=branch_col_ieObject,edge.width=branch_widths_ieObject,label.offset=3,cex=0.9,tip.col="white")
                                    axisPhylo(cex=2,col="white",col.axis="white")
                                          mtext("Mya",side=1,adj=-0.05,col="white")
            #Add labels
                        #Internal nodes
                        rootvalue<-ieObject$value_anc[1]
                        results_ordered<-ieObject[order(ieObject$node_desc),]
                        internalnodesize<-results_ordered$value_desc[(length(tree$tip.label)+1):((length(tree$tip.label)*2)-2)]
                        nodelabelsize<-(c(rootvalue,internalnodesize)*labelSize)^scaling
                                    nodelabels(pch=21,node=which(nodelabelsize>=0)+length(tree$tip.label),cex=abs(nodelabelsize[which(nodelabelsize>=0)]),bg="white")
                                    nodelabels(pch=21,node=which(nodelabelsize<0)+length(tree$tip.label),cex=abs(nodelabelsize[which(nodelabelsize<0)]),bg="white")
                        #Terminal nodes
                        tiplabelsize<-(results_ordered$value_desc[1:length(tree$tip.label)]*labelSize)^scaling
                                    tiplabels(pch=21,tip=which(tiplabelsize>=0),cex=abs(tiplabelsize[which(tiplabelsize>=0)]),bg="white")
                                    tiplabels(pch=21,tip=which(tiplabelsize<=0),cex=abs(tiplabelsize[which(tiplabelsize<=0)]),bg="white")
                                    par(bg="white")
                                                                }
