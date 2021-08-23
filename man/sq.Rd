\name{sq}

\alias{sq}

\title{Build matrix of binary vectors of a fixed lenght}

\description{
It creates all possible binary (with elements equal to 0 or 1) vectors of a specified length; these vectors are casted in a matrix having a number of columns equal to the vector length. If required, the binary vectors are constrained to have the same total (sum of their elements).
}

\usage{sq(J, t = NULL)}

\arguments{
  \item{J}{legnth of binary vectors}
  \item{t}{possible fixed total}
}
\value{
  \item{M}{matrix containing all binary vectors}
}

\author{Francesco Bartolucci, Antonio Forcina}

\examples{
# build matrix of all possible binary vectors of legnth 5
  M = sq(5)
  
# build matrix of all possible binary vectors of legnth 5 with total equal 2
  M = sq(5,2)
}
