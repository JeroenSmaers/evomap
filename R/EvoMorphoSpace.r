EvoMorphoSpace<-function(tree,X,Y,sequence,colors){
            #Plot the EvoSpace
                  for(i in sequence){
                  jpeg(sprintf("./EvoSpace_%05.2f.jpg",i),res=300,height=8,width=8,units="in")
                  layout(matrix(c(1), 1, 1, byrow = TRUE))
                        Epoch<-i
                              EvoSpace_X<-MorphoSpaceEpoch(X,Epoch)
                              EvoSpace_Y<-MorphoSpaceEpoch(Y,Epoch)
                              par(bg="black")
                                    #UnEqual limits
                                          lim_min_x<-min(c(X$value_desc,X$value_anc))
                                          lim_min_y<-min(c(Y$value_desc,Y$value_anc))
                                          lim_max_x<-max(c(X$value_desc,X$value_anc))
                                          lim_max_y<-max(c(Y$value_desc,Y$value_anc))
                                          plot(EvoSpace_X$value_desc,EvoSpace_Y$value_desc,pch=19,col="red",bg="red",xlim=c(lim_min_x,lim_max_x),ylim=c(lim_min_y,lim_max_y),ylab="",xlab="",col.axis="white",fg="white")
                                          title(main=paste0((i),"Mya"),col.main="white")
                  dev.off()
                                     }
                                     }
#