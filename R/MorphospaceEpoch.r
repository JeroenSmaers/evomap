#' Function to calculate ancestral values at particular snapshots in evolutionary time
#'
#' This function calculates ancestral values at a particular snapshot in evolutionare time
#' @param X ieObject
#' @return dataframe with rates for all branches and ancestral states for all nodes in the tree up to the indicated epoch
#' @export
MorphoSpaceEpoch<-function(X,Mya){

#Mya<-90

#Which branches fall onto the Mya line?
value_epoch_all<-c()
MorphoSpace_epoch<-X[which((X$time_node_anc>Mya)&((X$time_node_anc-X$BL)<=(Mya+.00000000001))),]
data_epoch<-split(MorphoSpace_epoch,1:nrow(MorphoSpace_epoch))

#Calculate ancestral state at the Mya line
for (i in 1:length(data_epoch)){

changePerUnit<-(data_epoch[[i]]$value_desc-data_epoch[[i]]$value_anc)/data_epoch[[i]]$BL
timeleft<-data_epoch[[i]]$BL-(Mya-(data_epoch[[i]]$time_node_anc-data_epoch[[i]]$BL))
value_epoch<-data_epoch[[i]]$value_anc+(changePerUnit*timeleft)
value_epoch_all<-rbind(value_epoch_all,value_epoch)

                                }
MorphoSpace_epoch<-cbind(MorphoSpace_epoch,value_epoch_all)
            return(MorphoSpace_epoch)

                                }
