\name{serial_parCube}
\alias{serial_parCube}
\title{
Serial 3-mode ParCube algorithm for memory resident tensors
}
\description{
ParCube uses sampling to reduce the problem then apply one of the tensor factorization methods: 
    \code{\link{cp_apr}} , \code{\link{cp_als}}, \code{\link{cp_nmu}} on the small tensor. 
    It is suitable for large sparse tensors. 
}
\usage{
serial_parCube(X, R, sample_factor, repetitions, opts = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
is a sparse tensor (a LIST containing subs, vals and size)
}
  \item{R}{
    The rank of the factorization
}
  \item{sample_factor}{
[s1 s2 s3] such that each sampled tensor is of size   [I/s1 J/s2 K/s3]
}
  \item{repetitions}{
    number of sampling repetitions
}
   \item{opts}{
    a list containing the options for the algorithm like maxiters:maximum iterations, tol:tolerance .. etc.
}

}
\details{
   opts: structure that stores options of the algorithm. For default,
       leave blank or use 'default_parcube_options()'.
           opts.p: percentage of common indices
           opts.nonneg: nonnegativity constraint enforced (binary)
           opts.loss: loss function (opts: 'fro' for Frobenius norm
               'kl' for KL-divergence)
           opts.weights: function of calculating sampling weights
               (opts: 'sum_abs' for sum of absolute values or 'sum_squares'
               for sum of squares)
           opts.normalize: normalize the factor matrices to unit norm per column (binary);
           opts.tolerance: the numerical tolerance of the algorithm
               (everything smaller than that is considered zero)
           opts.internal_tolerance: the tolerance of the solvers
               used interally
           opts.num_restarts: number of repetitions of each decomposition of each sample
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{lambda }{lambdas of Kruskal tensor}
\item{u }{LIST of A, B and C factor matrices for mode 1, 2 and 3 respectively.}
%% ...
}
\references{
   -Brett W. Bader, Tamara G. Kolda and others. 
     MATLAB Tensor Toolbox, Version [v3.0]. Available online at https://www.tensortoolbox.org, 2015.

   -Papalexakis, Evangelos E., Christos Faloutsos, and Nicholas D. Sidiropoulos. 
   "Parcube: Sparse parallelizable tensor decompositions." Machine Learning and Knowledge Discovery in Databases. 
   Springer Berlin Heidelberg, 2012. 521-536.
}
\author{Abdelmoneim Amer Desouki}

\seealso{
    \code{\link{cp_apr}} 
    \code{\link{cp_als}} 
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
 opts = default_parcube_options();
 opts[['loss']]='fro'
 opts[['nonneg']]='1'#nmu
 P1=serial_parCube(X,2,1,2,opts=opts)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ parCube }% use one of  RShowDoc("KEYWORDS")

