\name{random.gene.sets}
\alias{random.gene.sets}
\title{
Generate random gene sets from a list of genes that share similar structure with an input gene set
}
\description{
This function is used to generate randomized sets of genes that can be used to evaluate test results. This was mainly introduced to calibrate our test results for the giving input data. Need to investigate/explore more about its usage. Might want to remove it from our package if we are not going to use it.
}
\usage{
random.gene.sets(gene.sets, genes)
}
\arguments{
  \item{gene.sets}{The actual gene sets that would be used in a MegaCollapsing analysis.
}
  \item{genes}{A list of genes that would be used to generate randomzied gene sets.
}
}
\details{
Explain/explore why we had this function.
}
\value{
A randomized sets of genes that share the same set length with the input gene sets but with each gene set populated with random genes.
}
