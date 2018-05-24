#' GLS confidence intervals
#'
#' Computes confidence intervals for a GLS regression. When the variance covariance matrix is based on a phylogeny, this computes the confidence intervals for a phylogenetic regression (PGLS).
#' @param Y Dependent variable.
#' @param X Independent variable.
#' @param Sigma Matrix describing the covariance among data points. For the phylogenetic regression, this is based on a phylogeny (vcv(tree)).
#' @return GLS results.
#' @return Confidence intervals of the estimates ('CI').
#' @return Confidence intervals of unobserved values of Y ('CI.plot').
#' @references Smaers & Rohlf (2016) Testing species' deviations from allometric preductions using the phylogenetic regression. Evolution. In Press.
#' @examples See supplementary information Smaers & Rohlf (2016).
#' @export
gls.ci<-function(Y,X,Sigma){
      n<-length(X)
#Standardize sigma to reduce rounding errors
      tr<-sum(diag(Sigma))
      Sigma<-n*Sigma/tr
#Input
      invSigma<-solve(Sigma)
#pgls mean and variance
      X1<-rep(1,n)
      q<-2          # adjust if multivariate!!
      C1<-solve(t(X1)%*%invSigma%*%X1)
      Y_PGLSmean<-C1%*%t(X1)%*%invSigma%*%Y
      Y_PGLSdeviations = Y - Y_PGLSmean
      Y_PGLSvariance = (t(Y_PGLSdeviations)%*%invSigma%*%Y_PGLSdeviations)/(n-1)
      SE_Y_mean = sqrt(Y_PGLSvariance/n)
#Regression
      XX<-cbind(rep(1,n),X)
      C<-solve(t(XX)%*%invSigma%*%XX)
      w<-C%*%t(XX)%*%invSigma
      B<-w%*%Y
      Yhat<-XX%*%B
      Yresid=Y-Yhat
      Y_MSEresid<-(t(Yresid)%*%invSigma%*%Yresid)/(n-q)
      a<-B[1]
      b<-B[2]
#SE of statistics
      SEa<-sqrt(diag(C)*Y_MSEresid)[1]
      SEb<-sqrt(diag(C)*Y_MSEresid)[2]
#Save model
      intercept<-cbind(a,SEa)
      slope<-cbind(b,SEb)
      model<-rbind(intercept,slope)
      colnames(model)<-c("Estimate","Std.Error")
      rownames(model)<-c("intercept","slope")
#confint
      SEYhat<-sqrt(diag(XX%*%C%*%t(XX))%*%Y_MSEresid)
      CI<-cbind(X,Yhat,SEYhat)
            Lower2.5<-Yhat-qt(0.975,n)*SEYhat
            Lower5<-Yhat-qt(0.95,n)*SEYhat
            Upper5<-Yhat+qt(0.95,n)*SEYhat
            Upper2.5<-Yhat+qt(0.975,n)*SEYhat
      CI<-cbind(CI,Lower2.5)
      CI<-cbind(CI,Lower5)
      CI<-cbind(CI,Upper5)
      CI<-cbind(CI,Upper2.5)
      CI<-CI[order(CI[,1]),]
      colnames(CI)<-c("X","Yhat","SEYhat","Lower2.5","Lower5","Upper5","Upper2.5")
      CI<-as.data.frame(CI)
#confint extrapolated
      Xi<-seq(c(min(X)-abs(max(X))),to=c(abs(max(X))*5),length.out=100)
      Z<-c(Xi)
      ZZ<-cbind(rep(1,length(Z)),Z)
      SEYhat.Xi<-sqrt(diag(ZZ%*%C%*%t(ZZ))%*%Y_MSEresid)
      Yhat.Xi<-a+b*Xi
            Lower2.5.Yhat.Xi<-Yhat.Xi-qt(0.975,n)*SEYhat.Xi
            Lower5.Yhat.Xi<-Yhat.Xi-qt(0.95,n)*SEYhat.Xi
            Upper5.Yhat.Xi<-Yhat.Xi+qt(0.95,n)*SEYhat.Xi
            Upper2.5.Yhat.Xi<-Yhat.Xi+qt(0.975,n)*SEYhat.Xi
      CI.plot<-cbind(Xi,Yhat.Xi,SEYhat.Xi,Lower2.5.Yhat.Xi,Lower5.Yhat.Xi,Upper5.Yhat.Xi,Upper2.5.Yhat.Xi)
      colnames(CI.plot)<-c("X","Yhat","SEYhat","Lower2.5","Lower5","Upper5","Upper2.5")
      CI.plot<-as.data.frame(CI.plot)
Output<-list(model,CI,CI.plot)
names(Output)<-c("model","CI","CI.plot")
return(Output)
}

