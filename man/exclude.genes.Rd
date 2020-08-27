\name{exclude.genes}
\alias{exclude.genes}
\title{
Exclude a set of genes from input collapsing matrices before running MegaCollapsing Analysis.
}
\description{
It might be desire to remove genes that are causing trouble before running MegaCollapsing and this function can be used for that purpose.
}
\usage{
exclude.genes(input.data, exclude.file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.data}{The output from function call to read.collapsing.data.
}
  \item{exclude.file}{The file contains the list of genes to be excluded before analysis.
}
}
\details{
Add more about how it was implemented and the reason why we might want to do this.
}
\value{
A list that has the same structure as input.data.
}
