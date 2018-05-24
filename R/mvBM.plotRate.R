#' Lineage-specific rate estimation using multiple variance Brownian motion
#'
#' Plots lineage-specific rates using mvBM 
#' @param tree an object of class "phylo".
#' @param tree_mvBM an object of class "phylo". The rescaled tree from an mvBM procedure.
#' @param branches vector listing the branch numbers ('edge' numbers) for which the mvBM rate should be computed
#' @param sigma2Distr the MCMC distribution of sigma2 from an MCMC mvBM procedure
#' @return probability density distribution of mvBM rate
#' @references Smaers, Mongle & Kandler (2016) A multiple variance Brownian motion framework for estimating variable rates and inferring ancestral states. Biological Journal of the Linnean Society. 118 (1): 78-94.
#' @examples see smaerslab.com/mvBM/ 

#' @export
    mvBM.plotRate<-function(tree,tree_mvBM,branches,sigma2Distr,adjust,...){
    rateDistr<-sigma2Distr*(sum(tree_mvBM$edge.length[branches])/sum(tree$edge.length[branches]))
    Density_rateDistr<-density(rateDistr,adjust=adjust)
    plot(Density_rateDistr,,xlab=expression("Rate of evolution ("~sigma~"^2)"),...)
    polygon(Density_rateDistr,...)
    }
