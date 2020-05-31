\name{tensor2mat}
\alias{tensor2mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Reduce tensor to matrix
}
\description{
   Reduce the tensor to one matrix by summing all slices efficiently.
}
\usage{
tensor2mat(X,binary=FALSE,symmetrize=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
    List of tensor slices (Latent folding).
}
  \item{binary}{
    when true returns binary matrix (i.e. ignores multiple connections 
    between entities), default FALSE.
}
  \item{symmetrize}{
    when true all relations are considered symetric before adding the slices, 
    default FALSE.
}
}
\value{
  Adjacency matrix according to flags binary and symmetrize.
}

\author{
Abdelmoneim Amer Desouki
}
