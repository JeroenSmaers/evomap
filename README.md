# evomap

An R package for evolutionary mapping of continuous traits on phylogenies.

## Installation

```r
# install.packages("devtools")
devtools::install_github("JeroenSmaers/evomap")
```

## Functions

### Multiple variance Brownian motion (mvBM)

| Function | Description |
|---|---|
| `mvBM()` | Estimates rescaled branch lengths under a multiple-variance Brownian motion model |
| `mvBM.getRate()` | Extracts branch-specific evolutionary rates from mvBM output |
| `mvBM.plotRate()` | Plots branch-specific rates on a phylogeny |

```r
tree <- pbtree(n = 50)
x <- fastBM(tree, sig2 = 2)
BMsigma2 <- ace(x, tree, method = "REML")$sigma2[1]
results <- mvBM(x, tree, BMsigma2)

tree_mvBM <- tree
tree_mvBM$edge.length <- results$rBL
ace(x, tree_mvBM, method = "REML")
```

### Phylogenetic GLS

| Function | Description |
|---|---|
| `gls.ci()` | Confidence intervals for a pGLS regression |
| `gls.pi()` | Prediction intervals for a pGLS regression |
| `gls.ancova()` | F-ratio comparison of pGLS models (phylogenetic ANCOVA) |
| `pGLS.plotGrade()` | Overlay a pGLS regression line for a subset (grade) of taxa |

### Tree utilities

| Function | Description |
|---|---|
| `getEdges()` | All edges descendant from a node |
| `getNodes()` | All nodes descendant from a node |
| `getTips()` | All tips descendant from a node |
| `pruneSample()` | Prune a tree and dataset to a target group of taxa |

## Data

The package includes an example dataset of primate brain and body mass from Isler et al. (2008):

```r
data(primateData)
data(primateTree)
```

## Citation

If you use the mvBM framework, please cite:

> Smaers, Mongle & Kandler (2016) A multiple variance Brownian motion framework for estimating variable rates and inferring ancestral states. *Biological Journal of the Linnean Society* 118(1): 78–94.

## License

MIT
