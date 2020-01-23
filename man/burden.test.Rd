\name{burden.test}
\alias{burden.test}
\title{
Run regression/test using collapsing matrix with covaraiates.
}
\description{
Add a lot of explainations here to justify why we set up the modeld this way.
}
\usage{
burden.test(sample.list, gene.sets, syn.mat, non.syn.mat, X = NULL, Y = NULL)
}
\arguments{
  \item{sample.list}{ An object that hold information/covaraites about samples.
}
  \item{gene.sets}{The list of gene sets to be anallize.
}
  \item{syn.mat}{The syn collsping matrix to be used. We should have an option not to use this one.
}
  \item{non.syn.mat}{The non syn collsping matrix to be used.
}
  \item{X}{a charactor list of covariates to be used in regression models.
}
  \item{Y}{The response varaible. 1 for cases and 0 for controls.
}
}
\details{
Explain.
}
\value{
Return a list of test results and summary statistics. Will document this further once we decide what to keep.
}
