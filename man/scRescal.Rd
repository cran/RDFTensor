\name{scRescal}
\alias{scRescal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scalable RESCAL Factorization
}
\description{
implements improvements of RESCAL to be able to factorize graaphs with millions of entities.
Use parallelization and compact representation of sparse matrices.
Modified error calculation criteria that is more efficient and is based on the change of A & R instead
of calulating tensor difference.
}
\usage{
scRescal(X, rnk, ainit = "nvecs", verbose = 2, Ainit = NULL, Rinit = NULL, lambdaA = 0,
lambdaR = 0, lambdaV = 0, epsilon = 0.001, maxIter = 100, minIter = 1, P = list(), 
orthogonalize = FALSE, func_compute_fval = "compute_fit", retAllFact = FALSE, 
useQR = FALSE, ncores = 0, OS_WIN = FALSE, savePath = "", oneCluster = TRUE, 
useXc = FALSE, saveARinit = FALSE, saveIter = 0, dsname = "", maxNrows = 50000, 
generateLog = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{X}{
    is a sparse tensor as set of sparse matrices, one for every relation (predicate). (a LIST of SparseMatrix )
}
  \item{rnk}{
    The rank of the factorization
}
  \item{ainit}{
  the method used to initialize matrix A
}
  \item{verbose}{
 the level of messages to be displayed, 0 is minimal.
}

  \item{Ainit}{
 the initial value of matrix A, can be used to continue factorization from 
 previous results.
}
  \item{Rinit}{
    the initial value of R (the core tensor, as LIST of frontal slices)
}
  \item{lambdaA}{
    Regularization parameter for A factor matrix. 0 by default
}
  \item{lambdaR}{
    Regularization parameter for R_k factor matrices. 0 by default
}
  \item{lambdaV}{
    Regularization parameter for R_k factor matrices. 0 by default
}
  \item{epsilon}{
     error threshold
}
  \item{maxIter}{
  Maximum number of iterations
}
  \item{minIter}{
Minimum number of iterations
}
  \item{P}{
Not implemented
}
  \item{orthogonalize}{
Not implemented
}
  \item{func_compute_fval}{
function used to compute fit.
}

\item{retAllFact}{
       flag to return intermediate values of A & R
       }
\item{useQR}{ was suggested by Nickel in Factorizing Yago and 
implemented by Michail Huffman;    found to be converging more slowly.
}
\item{ncores}{ number of cores used to run in parallel, 0 means no paralellism
}
\item{useXc}{ compact the sparse tensor: require more space but much faster.
}
\item{saveIter}{ iterations in which A&R to be saved, default 0 :none
}
\item{saveARinit}{ option to save init A and R}
\item{maxNrows}{ used as limit to decide if the compact form of the 
sparse matrix is to be dense matrix (#rows <maxNrows) or to be sparse
used in updateA to consider the predicate having rows (in compact form) 
more than the given number as Big and hence return dense matrix.
Default 50000.
}
\item{OS_WIN }{ True when the operating system is windows, 
used to allow using Fork when running in parallel
}

  \item{savePath}{
 optional path to save intermediate results  into it.
}
  \item{oneCluster}{
  Boolean flag indicating that one cluster will be used in different steps when running in parallel 
}

%  \item{saveARinit}{
 % flag used to when the initial matrix A and R are required to be saved
%}

  \item{dsname}{
  the dataset name
}
  \item{generateLog}{
  save output when running in parallel to log file in current directory.
}
}

\value{

list(A=A, R=R, all_err, nitr=itr + 1, times=as.vector(exectimes)
Returns a LIST of the following:
\item{A}{The matrix A of the factorization ( n by r)}
\item{R}{The core tensor R the factorization as r (rank) matrices of ( r by r)}
\item{nitr}{number of iterations}
\item{times}{list of running times of each step.}
}
\references{
    -Maximilian Nickel, Volker Tresp, Hans-Peter-Kriegel,
      "A Three-Way Model for Collective Learning on Multi-Relational Data",
      ICML 2011, Bellevue, WA, USA

    -Maximilian Nickel, Volker Tresp, Hans-Peter-Kriegel,
      "Factorizing YAGO: Scalable Machine Learning for Linked Data"
      WWW 2012, Lyon, France
      
    -SynthG: mimicking RDF Graphs Using Tensor Factorization
}

\author{Abdelmoneim Amer Desouki}

\seealso{
    \code{\link{cp_apr}} 
    \code{\link{serial_parCube}} 
    \code{\link{cp_nmu}} 
    \code{\link{cp_als}} 
    \code{\link{rescal}} 
    \code{\link{rescal_Trp_Val}} 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
\dontrun{
   data('umls_tnsr')
   ntnsr=umls_tnsr
   tt=scRescal(ntnsr$X,rnk=10,ainit='nvecs',verbose=1,lambdaA=0,epsilon=1e-4,
                   lambdaR=0,ncores = 2,OS_WIN = TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RESCAL }% use one of  RShowDoc("KEYWORDS")
\keyword{ scalable RESCAL }% __ONLY ONE__ keyword per line
