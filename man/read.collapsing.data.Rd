\name{read.collapsing.data}
\alias{read.collapsing.data}
\title{
Read all data needed for a "paired burden test.
}
\description{
This function is currently used by CGR to run MegaCollapsing analsyis that that uses both synonamous and "non-synonamous" collapsing matrices. We need to come up with better names and descriptions for this kind of analysis.
}
\usage{
read.collapsing.data(samples, syn, non.syn, sample.column = "IID")
}
\arguments{
  \item{samples}{The input PED file to define the list of samples to be used in analysis.
}
  \item{syn}{The collapsing matrix for syn model.
}
  \item{non.syn}{The collapsing matrix for matrix to be acutally tested.
}
  \item{sample.column}{The column name that is used to identify the sample ID from input sample file.
}
}
\details{
Need a lot of description here to describe we are trying to do. Obvious bug: the passed in "sample.column"" was not used.
}
\value{
Return a lost of objects that are parsed and matched for "paired burden test"
  \item{sample.list }{The final sample list used for analysis}
  \item{syn}{The final matching collapsing matrix for syn model}
  \item{non.syn}{The final matching collapsing matrix for non-syn model}
}
