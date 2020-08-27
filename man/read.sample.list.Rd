\name{read.sample.list}
\alias{read.sample.list}
\title{
Read and format a sample list file.
}
\description{
This is a helper function to parse and standardize the sample names from sample list file. It assumes that the sample list file is tab demilited, has a column header and the sample ID column is labeled by column header "IID". To do: make this more generic and allow user to change the format.
}
\usage{
read.sample.list(sample_file, sample_ID)
}
\arguments{
  \item{sample_file}{ The input sample list file. Tab delimited, with column header and the sample ID column is indicated by column name IID.}
  \item{sample_ID}{ The sample column ID }
}
\details{
This is a helper function used only when the sample file has column headers and the some sample names are just integers. In this case, we will add an "X" in front the integer to make it a proper sample name. We might want to remove this function later.
}
\value{
Retrun a data frame represent the sample list.
}
