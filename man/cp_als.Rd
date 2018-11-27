\name{cp_als}
\alias{cp_als}
\encoding{utf8}
\title{
    Compute a CP decomposition using an alternating least-squares algorithm(als)
}
\description{
computes an estimate of the best rank-R PARAFAC model of a tensor X using an 
alternating least-squares algorithm
Translated from cp_als.m : MATLAB Tensor Toolbox
}
\usage{
cp_als(X, R, opts = list())
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
   \item{iters}{number of iterations.}
   \item{fit}{fraction explained by the model.}
}
\references{
    -Brett W. Bader, Tamara G. Kolda and others. 
     MATLAB Tensor Toolbox, Version [v3.0]. Available online at https://www.tensortoolbox.org, 2015.
}
\author{Abdelmoneim Amer Desouki}

\seealso{
    \code{\link{cp_apr}} 
    \code{\link{serial_parCube}} 
    \code{\link{rescal}} 
    \code{\link{cp_nmu}} 
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
P1=cp_als(X,2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ als }% __ONLY ONE__ keyword per line
