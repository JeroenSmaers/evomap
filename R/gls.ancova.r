#' pGLS model comparison using phylogenetic analysis of covariance (gls.ancova)
#'
#' Computes F-ratio test of GLS models that include different sets of indicator variables.
#' @param form Model formula. 
#' @param Sigma Matrix describing the covariance among data points. For the phylogenetic regression, this is based on a phylogeny (vcv(tree)). 
#' @param ReducedModel The GLS model that includes the fewest indicator variables. 
#' @param FullModel The GLS model that includes the fewest indicator variables. 
#' @return Anova table.
#' @references Smaers & Rohlf (2016) Testing species' deviations from allometric preductions using the phylogenetic regression. Evolution. 70 (5): 1145-1149.
#' @examples see smaerslab.com/gls.ancova/ 
#' @export

gls.ancova<-function (form, Sigma, ReducedModel, FullModel) 
{
    Y<-c(data[,which(colnames(data)==all.vars(form)[1])])
    n <- length(Y); tr <- sum(diag(Sigma)); Sigma <- n * Sigma/tr; invSigma <- solve(Sigma)
    df_ReducedModel <- c(length(ReducedModel[1, ])); df_FullModel <- c(length(FullModel[1, ]))

    C <- solve(t(ReducedModel) %*% invSigma %*% ReducedModel)
    B <- C %*% t(ReducedModel) %*% invSigma %*% Y
    SS_unexpl_ReducedModel <- (t(Y - (ReducedModel %*% B)) %*% invSigma %*% (Y - (ReducedModel %*% B)))
    MS_unexpl_ReducedModel <- SS_unexpl_ReducedModel/c(n-df_ReducedModel)

    C <- solve(t(FullModel) %*% invSigma %*% FullModel)
    B <- C %*% t(FullModel) %*% invSigma %*% Y
    SS_unexpl_FullModel <- (t(Y - (FullModel %*% B)) %*% invSigma %*% (Y - (FullModel %*% B)))
    MS_unexpl_FullModel <- SS_unexpl_FullModel/c(n-df_FullModel)

    Fs <- ((SS_unexpl_ReducedModel - SS_unexpl_FullModel)/(c(n-df_ReducedModel) - c(n-df_FullModel)))/(SS_unexpl_FullModel/c(n-df_FullModel))
    P <- 1 - pf(Fs, (c(n-df_ReducedModel) - c(n-df_FullModel)), c(n-df_FullModel))

    Output <- as.data.frame(cbind(rbind(df_FullModel, df_ReducedModel), 
        rbind(round(SS_unexpl_FullModel, 4), round(SS_unexpl_ReducedModel, 4)), rbind(round(MS_unexpl_FullModel, 4), round(MS_unexpl_ReducedModel, 4)), rbind(round(Fs, 4), ""), rbind(round(P, 4), "")))
    colnames(Output) <- c("df", "Sum Sq", "Mean Sum Sq", "F value", "Pr(>F)")
    rownames(Output) <- c("FullModel", "ReducedModel")
    return(Output)
}


