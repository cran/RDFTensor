\name{cp_nmu}
\alias{cp_nmu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute nonnegative CP with multiplicative updates(NMU)
}
\description{
computes an estimate of the best rank-R PARAFAC model of a tensor X with 
nonnegative constraints on the factors.
This version uses the Lee & Seung multiplicative updates from their NMF algorithm.
Translated from cp_nmu.m : MATLAB Tensor Toolbox
}
\usage{
cp_nmu(X, R, opts = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
is a sparse tensor (a LIST containing subs, vals and size)
}
  \item{R}{
    The rank of the factorization
}
  \item{opts}{
    a list containing the options for the algorithm like maxiters:maximum iterations, tol:tolerance .. etc.
}
}
\value{
%%  ~Describe the value returned
  \item{P }{the factorization of X as a LIST representing Kruskal Tensor (lambda and u)}
   \item{Uinit}{the initial solution}
   \item{stats}{statistics about the solution like the running time of each step and the error.}
   \item{fit}{fraction explained by the model.}
}
\references{
    -Brett W. Bader, Tamara G. Kolda and others. 
     MATLAB Tensor Toolbox, Version [v3.0]. Available online at https://www.tensortoolbox.org, 2015.

    -Lee, Daniel D., and H. Sebastian Seung. "Algorithms for non-negative matrix factorization." 
     In Advances in neural information processing systems, pp. 556-562. 2001.
}
\author{Abdelmoneim Amer Desouki}

\seealso{
    \code{\link{cp_apr}} 
    \code{\link{serial_parCube}} 
    \code{\link{rescal}} 
    \code{\link{cp_als}} 
}
\examples{
subs=matrix(c(5,1,1,
              3,1,2,
              1,1,3,
              2,1,3,
              4,1,3,
              6,1,3,
              1,1,4,
              2,1,4,
              4,1,4,
              6,1,4,
              1,2,1,
              3,2,1,
              5,2,1),byrow=TRUE,ncol=3)

X=list(subs=subs,vals=rep(1,nrow(subs)),size=c(6,2,4))
set.seed(123)
P1=cp_nmu(X,2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ NMF }% use one of  RShowDoc("KEYWORDS")
\keyword{ NMU }% __ONLY ONE__ keyword per line
