\name{cp_apr}
\alias{cp_apr}
\title{
Compute nonnegative CP with alternating Poisson regression(CP_APR)}
\description{
computes an estimate of the best rank-R CP model of a
   nonnegative tensor X using an alternating Poisson regression. This is
   most appropriate for sparse count data (i.e., nonnegative integer 
values) because it uses Kullback-Liebler divergence.}
\usage{
cp_apr(X, R, opts = list())
}
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
\details{
Different algorithm variants are available (selected by the 'alg'
   parameter):
    'pqnr' - row subproblems by projected quasi-Newton (default)
    'pdnr' - row subproblems by projected damped Hessian
    'mu'   - multiplicative update
   Additional input parameters for algorithm 'mu':
      'kappa'         - Offset to fix complementary slackness {100}
      'kappatol'      - Tolerance on complementary slackness {1.0e-10}

   Additional input parameters for algorithm 'pdnr':
      'epsActive'     - Bertsekas tolerance for active set {1.0e-8}
      'mu0'           - Initial damping parameter {1.0e-5}
      'precompinds'   - Precompute sparse tensor indices {TRUE}
      'inexact'       - Compute inexact Newton steps {TRUE}

   Additional input parameters for algorithm 'pqnr':
      'epsActive'     - Bertsekas tolerance for active set {1.0e-8}
      'lbfgsMem'      - Number vector pairs to store for L-BFGS {3}
     'precompinds'   - Precompute sparse tensor indices {TRUE}

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{M }{the factorization of X as a LIST representing Kruskal Tensor (lambda and u)}
   \item{Minit}{the initial solution}
   \item{output}{statistics about the solution like the running time of each step and the error.}
}
\references{
    -Brett W. Bader, Tamara G. Kolda and others. 
     MATLAB Tensor Toolbox, Version [v3.0]. Available online at https://www.tensortoolbox.org, 2015.

    -E. C. Chi and T. G. Kolda. On Tensors, Sparsity, and Nonnegative
        Factorizations, SIAM J. Matrix Analysis,  33(4):1272-1299, Dec. 2012,
        http://dx.doi.org/10.1137/110859063
   -S. Hansen, T. Plantenga and T. G. Kolda, Newton-Based Optimization
    for Kullback-Leibler Nonnegative Tensor Factorizations, 
    Optimization Methods and Software, 2015, 
    http://dx.doi.org/10.1080/10556788.2015.1009977        
}
\author{Abdelmoneim Amer Desouki}

\seealso{
    \code{\link{cp_nmu}} 
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
    set.seed(12345)#for reproducability
    P1=cp_apr(X,2,opts=list(alg='mu'))
    print(P1$M)
    set.seed(12345)#for reproducability
    P2=cp_apr(X,2,opts=list(alg='pdnr'))
    print(P2$M)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ APR }% use one of  RShowDoc("KEYWORDS")
\keyword{ Tensor Factorization }% __ONLY ONE__ keyword per line
