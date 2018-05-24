#' Multiple variance Brownian motion estimation
#'
#' Computes rescaled branch lengths using multiple variance Brownian motion as described in Smaers et al. (2016) 
#' @param data a vector of tip values for species; should be in the same order as tiplabels in the tree
#' @param tree an object of class "phylo".
#' @param sigma sig2 value from a Bayesian MCMC run using standard Brownian motion (e.g. ?anc.Bayes)
#' @return dataframe with rescaled branch lengths (rBL) for all branches in the tree
#' @references Smaers, Mongle & Kandler (2016) A multiple variance Brownian motion framework for estimating variable rates and inferring ancestral states. Biological Journal of the Linnean Society.
<<<<<<< HEAD
#' @examples tree<-pbtree(n=50) # simulate tree ('pbtree' requires the 'phytools' package)
#' x<-fastBM(tree,sig2=2) # simulate values ('pbtree' requires the 'phytools' package)
#' BMsigma2<-ace(x,tree,method="REML")$sigma2[1] # get BM sigma2 ('ace' requires the 'ape' package)
#' mvBMresults<-mvBM(x,tree,BMsigma2) # calculate rescaled branch lengths using mvBM
#' tree_mvBM<-tree
#' tree_mvBM$edge.length<-mvBMresults$rBL # create new tree with rescaled branch lengths
#' plot(tree_mvBM) # plot mvBM tree
#' ace(x,tree_mvBM,method="REML") # get ancestral estimates using mvBM tree 
#' # To obtain the mvBM branch specific rate estimate: calculate sig2 for the trait in question using the mvBM tree; multiply the mvBM sig2 with the mvBM rescaled branch length of the lineage of interest; divide this value by the phylogenetic branch length of the lineage of interest. 
#' # MCMC posterior distributions can be obtained by using anc.Bayes for the above calculations. ('anc.Bayes' requires the 'phytools' package)

#' @export
mvBM <- function (data, tree, sigma2) 
{
    data_original <- data
    N = length(data)
    phy.matrix = data.frame(tree$edge, tree$edge.length, data[tree$edge[, 
        2]])
    names(phy.matrix) = c("Anc", "Desc", "Length", "Value")
    nodes_extant = 1:N
    nodes_extinct = (N + 1):(N + (N - 1))
    extant_values <- data.frame(data, nodes_extant)
    Pk_values <- c()
    
    dist.tree <- dist.nodes(tree)
    
    for (j in nodes_extinct) {
        nominator = c()
        denominator = c()

        for (i in nodes_extant) {
            nominator = rbind(nominator, (extant_values[i, 1]/dist.tree[i, 
                j]^2))
            denominator = rbind(denominator, (1/dist.tree[i, 
                j]^2))
        }
 
        Pk = sum(nominator)/sum(denominator)
        Pk_values = rbind(Pk_values, Pk)
    }
    
    
    Pk <- c()
    Pk <- cbind(Pk, Pk_values)
    Pk <- cbind(Pk, nodes_extinct)
    colnames(Pk) <- c("value", "nodes")
    rownames(Pk) <- 1:length(nodes_extinct)
    Pk <- as.data.frame(Pk)
    nodes_extinct_reverse = sort(nodes_extinct, decreasing = TRUE)
    rBL <- c()
    rBL_edge <- c()
    for (i in nodes_extinct_reverse) {
        sister_branches <- which(phy.matrix$Anc == i)
        X1 <- phy.matrix$Value[sister_branches[1]]
        X2 <- phy.matrix$Value[sister_branches[2]]
        Y1 <- sqrt(phy.matrix$Length[sister_branches[1]])
        Y2 <- sqrt(phy.matrix$Length[sister_branches[2]])
        Pk_anc <- Pk$value[which(Pk$nodes == i)]
        Ax <- (Pk_anc + X1 + X2)/3
        T1 <- ((X1 - Ax)^2)/sigma2 + Y1^2
        T2 <- ((X2 - Ax)^2)/sigma2 + Y2^2
        phy.matrix$Value[which(phy.matrix$Desc == i)] <- Ax
        rBL <- c(rBL, T1)
        rBL <- c(rBL, T2)
        rBL_edge <- c(rBL_edge, sprintf("%05.f", sister_branches[1]))
        rBL_edge <- c(rBL_edge, sprintf("%05.f", sister_branches[2]))
        names(rBL) <- rBL_edge
    }
    Output <- data.frame(tree$edge, tree$edge.length, rBL[sort(names(rBL))])
    names(Output) <- c("node_anc", "node_desc", "BL", "rBL")
    return(Output)
=======
#' @examples tree<-pbtree(n=50)
#' x<-fastBM(tree,sig2=2) # simulate using fastBM
#' X<-anc.Bayes(tree,x,ngen=10000) # sample ancestral states
#' estimates<-colMeans(X[21:nrow(X),]) # get estimates, excluding burnin
#' sigma2<-mean(X[21:nrow(X),2]) # get sigma2, excluding burnin
#' mvBMresults<-mvBM(x,tree,sigma2) # calculate rescaled branch lengths using mvBM
#' tree_mvBM<-tree
#' tree_mvBM$edge.length<-mvBMresults$rBL # create new tree with rescaled branch lengths
#' X_mvBM<-anc.Bayes(tree_mvBM,x,ngen=10000)  # use new tree in Bayesian MCMC run
#' estimates_mvBM<-colMeans(X_mvBM[21:nrow(X_mvBM),]) # get estimates, excluding burnin

#' @export

mvBM <- function(data,tree,sigma2){

#Tree structure
  data_original<-data
            N=length(data)
      #Make phylo matrix and list nodes
            phy.matrix=data.frame(tree$edge,tree$edge.length,data[tree$edge[,2]])
            names(phy.matrix)=c("Anc","Desc","Length","Value")
            nodes_extant=1:N
            nodes_extinct=(N+1):(N+(N-1))

#Global estimate

      extant_values<-data.frame(data,nodes_extant)
      Pk_values<-c()

      for(j in nodes_extinct){
      #calculate nominator
            nominator=c()
                  for(i in nodes_extant){
                        nominator=rbind(nominator,(extant_values[i,1]/(dist.nodes(tree)[i,j]^2)))
                                        }
      #calculate denominator
            denominator=c()
                  for(i in nodes_extant){
                        denominator=rbind(denominator,(1/(dist.nodes(tree)[i,j]^2)))
                                        }
      #save Pk
            Pk=sum(nominator)/sum(denominator)
            Pk_values=rbind(Pk_values,Pk)
                             }

      Pk<-c()
      Pk<-cbind(Pk,Pk_values)
      Pk<-cbind(Pk,nodes_extinct)
      colnames(Pk)<-c("value","nodes")
      rownames(Pk)<-1:length(nodes_extinct)
      Pk<-as.data.frame(Pk)

#Rescaled branch lengths
      nodes_extinct_reverse=sort(nodes_extinct,decreasing=TRUE)
            rBL<-c()
            rBL_edge<-c()

      for(i in nodes_extinct_reverse) {
            sister_branches<-which(phy.matrix$Anc==i)

            X1<-phy.matrix$Value[sister_branches[1]]
            X2<-phy.matrix$Value[sister_branches[2]]
            Y1<-sqrt(phy.matrix$Length[sister_branches[1]])
            Y2<-sqrt(phy.matrix$Length[sister_branches[2]])

            Pk_anc<-Pk$value[which(Pk$nodes==i)]

            Ax<-(Pk_anc+X1+X2)/3

            T1<-((X1-Ax)^2)/sigma2+Y1^2
            T2<-((X2-Ax)^2)/sigma2+Y2^2

            phy.matrix$Value[which(phy.matrix$Desc==i)]<-Ax

            rBL<-c(rBL,T1)
            rBL<-c(rBL,T2)
            rBL_edge<-c(rBL_edge,sprintf("%05.f",sister_branches[1]))
        rBL_edge<-c(rBL_edge,sprintf("%05.f",sister_branches[2]))
            names(rBL)<-rBL_edge

                                }
#Results
        results<-data.frame(tree$edge,tree$edge.length,rBL[sort(names(rBL))])
            names(results)<-c("node_anc","node_desc","BL","rBL")
        return(results)

>>>>>>> parent of 1504e9a... Package Update
}
