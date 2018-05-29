#' pGLS plot of a particular grade
#'
#' Plots the pGLS model of a particular subset of a comparative dataset
#' @param Yvar column name of the dependent variable
#' @param Xvar column name of the independent variable
#' @param data a vector of tip values for species; should be in the same order as tiplabels in the tree
#' @param tree an object of class "phylo".
#' @param model of evolution to be used in pGLS ("BM" or "lambda").
#' @param group vector listing the tip labels for which the pGLS model should be plotted
#' @return pGLS abline and points of the group of interest
#' @examples see smaerslab.com/software

#' @export
pGLS.plotGrade<-function(Yvar,Xvar,data,tree,model,group,...){
    dataTemp<-pruneSample(na.omit(data),tree,group)$data
    treeTemp<-pruneSample(dataTemp,tree,group)$tree
    Y<-dataTemp[,which(colnames(dataTemp)==paste(Yvar))]
    X<-dataTemp[,which(colnames(dataTemp)==paste(Xvar))]
    dataGLS<-as.data.frame(cbind(Y,X)); rownames(dataGLS)<-rownames(dataTemp)
switch(model,
  BM={pGLSTemp<-gls(Y~X,dataGLS,correlation=corBrownian(phy=treeTemp))},
  lambda={pGLSTemp<-gls(Y~X,dataGLS,correlation=corPagel(1,phy=treeTemp,fixed=FALSE))})
      a<-summary(pGLSTemp)$tTable[1,1]
      b<-summary(pGLSTemp)$tTable[2,1]
      lines(c(min(X),max(X)),c((a+b*min(X)),(a+b*max(X))),...)
      points(X,Y,...)
  }
