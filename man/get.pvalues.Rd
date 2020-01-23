\name{get.pvalues}
\alias{get.pvalues}
\title{
Run permutation based Fisher Test on (Mega)collapsing matrix.
}
\description{
This function allows one to run permutation based FET tests. Not sure if we need this one in this package though.
}
\usage{
get.pvalues(matrix, is.case, n.permutations = 1000)
}
\arguments{
  \item{matrix}{A (mega) collapsing matrix.}
  \item{is.case}{A boolean vector indicating case/control status.}
  \item{n.permutations}{Number of permutations to run the tests.}
}
\details{
Running permutation based FET test from (mega)Collapsing matrix. Will explain more if we will keep this function.
}
\value{
A list that holds both FET p values and permutatied p-values. These p-values can be used to do permutation based QQ plot
}
