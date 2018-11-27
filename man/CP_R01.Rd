\name{CP_R01}
\alias{CP_R01}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
inverse of real-number-CP-decompositions to Binary
}
\description{
gets the best reconstruction of a CP-factorization by trying different thresholds 
on the result. If the size of the tensor is too big sampling is done to get estimates
of TP, FP.
}
\usage{
CP_R01(X, P, 
    pthr = c(1e-06,1e-04,0.001,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.55,0.6, 0.65, 0.7, 0.8),
    cntNnz = 200, startSize = 1e+07)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
the original tensor as sptensor list(subs,vals,size,nnz)
}  
\item{P}{
a LIST containing the CP transformation
}
  \item{pthr}{
list of threshold values to be tried
}

  \item{cntNnz}{
  If |X| is the number of non-zeros in data X, the sampled locations will be cntNnz*|X| of 0s.
}
  \item{startSize}{
if size of tensor < startSize all values are calculated, default 1e7.
}
}
\value{
A LIST containing the TP, FP, FN and threshold value also the result of best threshold
 }

\author{
Abdelmoneim Amer Desouki
}

\seealso{
 \code{\link{cp_apr}} 
    \code{\link{serial_parCube}} 
    \code{\link{rescal}} 
    \code{\link{rescal_01}} 
    \code{\link{cp_nmu}} 
}
\examples{
trp=rbind(
    cbind('Alex',  'loves', 'Don'),
    cbind('Alex',  'loves', 'Elly'),
    cbind('Alex',  'hates', 'Bob'),
    cbind('Don',   'loves', 'Alex'),
    cbind('Don',   'hates', 'Chris'),
    cbind('Chris', 'hates', 'Bob'),
    cbind('Bob',   'hates', 'Chris'),
    cbind('Elly',  'hates', 'Chris'),
    cbind('Elly',  'hates', 'Bob'),
    cbind('Elly',  'loves', 'Alex')
    )
######
# form tensor as a set of frontal slices (Predicate mode)
    tnsr=getTensor(trp)
    subs=getTnsrijk(tnsr$X)
    X=list(subs=subs,vals=rep(1,nrow(subs)),size=c(5,2,5))
    normX=sqrt(sum(X$vals))
    set.seed(123)
    # NMU decomposition with rank 2
    P1=cp_nmu(X,2)
    res=CP_R01(X,P1[[1]])
}