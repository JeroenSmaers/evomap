#' Function to allocate the rate distribution between 2 variables to rate quadrants
#'
#' This function allocates the rate distribution between 2 variables to rate quadrants
#' @param ieObject object comprising the output of an ie analysis
#' @return list of ratedistribution elements and branches for each quadrant
#' @export

ieRateDistr<-function(Y,X,neutral){
            #relative rate
                    YvX<-Y[1:3]
                    D<-Y$rate-X$rate
                    YvX<-cbind(YvX,Y$rate)
                    YvX<-cbind(YvX,X$rate)
                    YvX<-cbind(YvX,D)
                    MA<-sqrt((D^2)-((sqrt(2*(D^2))/2)^2))
                    MAresidual<-MA
                    YvX<-cbind(YvX,MAresidual)
                          negatives<-which(YvX$D<0)
                          for(i in negatives){
                          YvX$MAresidual[i]<-YvX$MAresidual[i]*(-1)}
                    colnames(YvX)<-c("node_anc","node_desc","BL","Y_rate","X_rate","D","MAresidual")

            #neutral
                  neutral_max<-neutral
                  neutral_min<-(-neutral)
            #accellerated increase
                  AI<-which(YvX$MAresidual>neutral_max & Y$rate>0 & X$rate>0 & Y$rate>neutral)
            #decellerated increase
                  DI<-which(YvX$MAresidual<neutral_min & Y$rate>0 & X$rate>0 & X$rate>neutral)
            #decellerated decrease
                  DD<-which(YvX$MAresidual>neutral_max & Y$rate<0 & X$rate<0 & X$rate<(-neutral))
            #accellerated decrease
                  AD<-which(YvX$MAresidual<neutral_min & Y$rate<0 & X$rate<0 & Y$rate<(-neutral))
            #separation
                  S1<-c(which(YvX$MAresidual>neutral_max & Y$rate>0 & X$rate<0 & Y$rate>neutral),
                        which(YvX$MAresidual>neutral_max & Y$rate>0 & X$rate<0 & Y$rate<=neutral & X$rate< (-(neutral))))
                  S2<-c(which(YvX$MAresidual<neutral_min & Y$rate<0 & X$rate>0 & Y$rate<(-neutral)),
                        which(YvX$MAresidual<neutral_min & Y$rate<0 & X$rate>0 & Y$rate>=(-neutral) & X$rate> neutral))
            #all
                  all<-c(S1,AI,DI,S2,AD,DD)

      #Return
            results<-list(YvX,S1,AI,DI,S2,AD,DD,all)
            names(results)<-c("YvX","S1","AI","DI","S2","AD","DD","All")
            return(results)

}
