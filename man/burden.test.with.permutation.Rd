\name{burden.test.with.permutation}
\alias{burden.test.with.permutation}
\title{
Permutation based Burden Test
}
\description{
Don't do this unless we really want to do it. This will run permutation based on Burden Test for given gene list and input mega/collapsing matrix.
}
\usage{
burden.test.with.permutation(sample.list, gene.set, syn.mat, non.syn.mat,
                            X = NULL, Y = NULL, n.permutations = 100)
}
\arguments{
    \item{sample.list}{ An object that hold information/covaraites about samples.[Copied from burden test]
}
  \item{gene.set}{The list of gene sets to be anallize.
}
  \item{syn.mat}{The syn collsping matrix to be used. We should have an option not to use this one.
}
  \item{non.syn.mat}{The non syn collsping matrix to be used.
}
  \item{X}{a charactor list of covariates to be used in regression models.
}
  \item{Y}{The response varaible. 1 for cases and 0 for controls.
}
  \item{n.permutations}{Number of permutations to use in this test. Default to 100.
}
}

\value{
A list that holds the actual burden test results and permutated results. Add more if we are to keep this function.
}
