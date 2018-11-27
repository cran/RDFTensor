\name{spt_mttkrp}
\alias{spt_mttkrp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matricized tensor times Khatri-Rao product for ktensor
}
\description{
V = MTTKRP(X,U,n) efficiently calculates the matrix product of the
   n-mode matricization of X with the Khatri-Rao product of all
   entries in U, a cell array of matrices, except the nth.  
   The tensor is considered to be sparse.
}
\usage{
    spt_mttkrp(X, U, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
 sparse tensor
}
  \item{U}{
list of matrices to be multiplied by X
}
  \item{n}{
  the mode used
}
}
\value{
a dense matrix of size nrow=X$size[n],ncol=R
}
\references{
    -Brett W. Bader, Tamara G. Kolda and others. 
     MATLAB Tensor Toolbox, Version [v3.0]. Available online at https://www.tensortoolbox.org, 2015.
}
\author{
Abdelmoneim Amer Desouki
}
