\name{read.data}
\alias{read.data}
\title{
Read gene/sample variant count matrix.
}
\description{
This is a helper function to extract gene/sample variant count matrix based on an input sample list and also optionally a list of genes to be used as filter.
}
\usage{
read.data(sample.file, matrix.file, filter.list = NULL, sample.column = 2,
          case.control.column = 6)
}
\arguments{
  \item{sample.file}{An input sample list file. The default format is a standard PED file that is tab delimited, without column header and has second column indicating the sample ID/name and the 6 column indicating the control/case status, with 1 for control sample and 2 for case sample. User can also provide a tab delimited text file and indicating the column index for sample ID and case/control status.
}
  \item{matrix.file}{The gene/sample count (collapsing) matrix to be read. A gene/sample count matrix is a tab delimited spreadsheet like text file with each column represents a sample, and each row represents a gene. Gene names will be stored in the first column an the sample names will be stored in the first row.}
  \item{filter.list}{Optional. A flat text file with each row gives a gene that is to be extracted from the input matrix file. By default, all genes found in the matrix file will be used.  }
  \item{sample.column}{Optional. If the sample ID column from the sample file is not the second column, user will need to provide the index(one based) for actutal sample ID column. }
  \item{case.control.column}{Optional. If the case/control status column from the sample file is not the second column, user will need to provide the index(one based) for actutal sample ID column.}
}

\value{Return a list containing the extracted case/control status and the selected gene/sample count matrix.
    \item{data }{A matrix for the extracted gene/sample count}
    \item{is.case }{A boolean vector indicating the case/control status}
}
