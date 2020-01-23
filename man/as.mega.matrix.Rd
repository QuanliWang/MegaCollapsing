\name{as.mega.matrix}
\alias{as.mega.matrix}
\title{Generate mega collapsing matrix/matrices for a list of gene sets.
}
\description{
Taking collapsing matrix/matrices and a list of gene sets, this function "collapses" gene/sample level variant information to gene set level varaint information, which can further be used for (mega) collapsing analysis.
}
\usage{
as.mega.matrix(gene.sets, input.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene.sets}{ The gene sets object that hold list of gene sets. It is typically obtained by calling function read.gene.sets.
%%     ~~Describe \code{gene.sets} here~~
}
  \item{input.data}{The objects that holds the collapsing matrix/matrices. It is typically obtained by calling function read.collapsing.data.
}
}
\details{
A lot of explaination here if we wanted. Or provide a link to the method.
}
\value{
Return a list of objects that represent the mega collapsing matrix/matices.
}
