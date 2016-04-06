#' GLS model comparison
#'
#' Computes F-ratio test of GLS models that include different sets of indicator variables. See supplementary information Smaers & Rohlf (2016) for more information.
#' @param Y Dependent variable. 
#' @param X Independent variable. 
#' @param Sigma Matrix describing the covariance among data poitns. For the phylogenetic regression, this is based on a phylogeny (vcv(tree)). 
#' @param ReducedModel The GLS model that includes the fewest indicator variables. 
#' @param FullModel The GLS model that includes the fewest indicator variables. 
#' @return Anova table.
#' @references Smaers & Rohlf (2016) Testing species' deviations from allometric preductions using the phylogenetic regression. Evolution. In Press.
#' @examples See supplementary information Smaers & Rohlf (2016).

#' @export

gls.ancova<-function(Y,X,Sigma,ReducedModel,FullModel){
            n<-length(X); tr<-sum(diag(Sigma)); Sigma<-n*Sigma/tr; invSigma<-solve(Sigma)
            df_ReducedModel<-(n-length(ReducedModel[1,])); df_FullModel<-(n-length(FullModel[1,]))
        #Model 1
              C<-solve(t(ReducedModel)%*%invSigma%*%ReducedModel)
              B<-C%*%t(ReducedModel)%*%invSigma%*%Y
            SS_unexpl_ReducedModel<-(t(Y-(ReducedModel%*%B))%*%invSigma%*%(Y-(ReducedModel%*%B)))
            MS_unexpl_ReducedModel<-SS_unexpl_ReducedModel/df_ReducedModel
        #Model 2
              C<-solve(t(FullModel)%*%invSigma%*%FullModel)
              B<-C%*%t(FullModel)%*%invSigma%*%Y
            SS_unexpl_FullModel<-(t(Y-(FullModel%*%B))%*%invSigma%*%(Y-(FullModel%*%B)))
            MS_unexpl_FullModel<-SS_unexpl_FullModel/df_FullModel
        #F
            Fs<-((SS_unexpl_ReducedModel-SS_unexpl_FullModel)/(df_ReducedModel-df_FullModel))/(SS_unexpl_FullModel/df_FullModel)
            P<-1-pf(Fs,(df_ReducedModel-df_FullModel),df_FullModel)

results<-as.data.frame(cbind(rbind(df_FullModel,df_ReducedModel),rbind(round(SS_unexpl_FullModel,4),round(SS_unexpl_ReducedModel,4)),rbind(round(MS_unexpl_FullModel,4),round(MS_unexpl_ReducedModel,4)),rbind(round(Fs,4),""),rbind(round(P,4),"")))
colnames(results)<-c("df","Sum Sq","Mean Sum Sq","F value","Pr(>F)")
rownames(results)<-c("FullModel","ReducedModel")
return(results)
}

