\name{rescal_01}
\alias{rescal_01}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
transformation of the real-number-RESCAL-decompositions to Binary
}
\description{
applies different thresholds to get a number of triples from A and R (result of RESCAL)
decomposition as scale_fact*the_number_of_triples_in_original_tensor.
}
\usage{
rescal_01(X, A, R, scale_fact = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
is a sparse tensor as set of sparse matrices, one for every relation (predicate). (a LIST of SparseMatrix )
}
  \item{A}{
the A matrix returned by RESCAL factorization
}
  \item{R}{
the R LIST returned by RESCAL factorization
}
  \item{scale_fact}{
    scale of the number of triples to be considered in the result. 
    When it is 1 then a threshold will be taken to get the same number of triples in each slice as the original tensor.
}
}

\value{
 a LIST
 \item{X_ }{The reconstructed tensor as a set of frontal slices.}
 \item{tp }{the number of true positives}
 \item{fp }{the number of false positives}
 \item{fn }{the number of false negatives}
 \item{sr }{the details of each slice i.e the number of tp, fn, fp, etc}

}
\references{
     -Maximilian Nickel, Volker Tresp, Hans-Peter-Kriegel,
      "A Three-Way Model for Collective Learning on Multi-Relational Data",
      ICML 2011, Bellevue, WA, USA
}
\author{
Abdelmoneim Amer Desouki
}
\examples{
    X1=matrix(c(1,0,0,0,0, 0,1,0,0,0, 0,0,1,1,0, 0,0,0,0,1, 1,0,0,0,0),byrow=TRUE,nrow=5,ncol=5)
    X2=matrix(c(0,1,0,1,1, 1,0,0,1,0, 0,1,0,1,1, 0,0,0,0,1, 0,0,1,0,0),byrow=TRUE,nrow=5,ncol=5)
    X2_=matrix(c(0,1,0,1,1, 1,0,0,1,0, 0,0,0,0,0, 0,0,0,0,1, 0,0,1,0,0),byrow=TRUE,nrow=5,ncol=5)
    X=list(t(X1),t(X2),t(X2_))

    N=nrow(X1)
    Xs=list()
    for(s in 1:length(X)){
      aa=which(X[[s]]==1,arr.ind=TRUE)
      Xs[[s]]=Matrix::sparseMatrix(x=rep(1,nrow(aa)),i=aa[,1],j=aa[,2],dims=c(N,N))
    }

    print(Xs)

    rf=rescal(Xs,2)
    A=rf$A
    R=rf$R
    tmp=rescal_01(Xs,A,R,scale_fact=1.5)#generate 1.5*original number of triples
    print(sprintf('Precision:\%.4f, Recall:\%.4f',tmp$tp/(tmp$tp+tmp$fp),tmp$tp/(tmp$tp+tmp$fn)))
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{rescal}} 
    \code{\link{CP_01}} 
}
