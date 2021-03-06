% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Linarius.R
\name{col.edge.phylo}
\alias{col.edge.phylo}
\title{Coloring terminal branch of trees}
\usage{
col.edge.phylo(phy, groups, color = c(2, 3, 4, 5), root.color = 1)
}
\arguments{
\item{phy}{A tree of class phylo}

\item{groups}{A vector indicating what should be of the same color}

\item{color}{A vector indicating the color to choose}

\item{root.color}{A single color, to use for non terminal edge}
}
\description{
Generate colors for terminal edge branch of a phylogenic tree
}
\note{
This function should be used with the package ade4
}
\examples{
data(Birch)
require(vegan)
require(ade4)
require(ape)
dist=vegdist(Betula, method="jaccard", binary=TRUE, diag=FALSE, upper=FALSE, na.rm = FALSE)
dendro=hclust(dist,"complete")
dendro2<-as.phylo(dendro,cex=0.5)
plot(dendro2,,type="u",lab4ut="axial",font=1, adj=1,label.offset = 0.01,edge.color = col.edge.phylo(dendro2,code$Ploidy),tip.color = code$Ploidy)
}
\seealso{
ade4
}

