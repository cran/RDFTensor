\name{rescal}
\alias{rescal}

\title{
    RESCAL: Tensor Factorization.
}
\description{
RESCAL-ALS algorithm to compute the RESCAL tensor factorization.
The solution is a matrix and a core tensor.
RESCAL factors a (usually sparse) three-way tensor X such that each frontal slice X_k is factored into
    X_k = A * R_k * A^T
}
\usage{
rescal(X, rnk,ainit = 'nvecs',verbose=2,Ainit=NULL,Rinit=NULL,lambdaA=0,lambdaR=0, 
      lambdaV=0,epsilon=1e-3,maxIter=100, minIter=1, P = list(),orthogonalize=FALSE,
      func_compute_fval='compute_fit')
}

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
 the initial value of matrix A.
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
}
\author{Abdelmoneim Amer Desouki}

\seealso{
    \code{\link{cp_apr}} 
    \code{\link{serial_parCube}} 
    \code{\link{cp_nmu}} 
    \code{\link{cp_als}} 
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
    Tensor_error(Xs,A,R)
    tmp=rescal_01(Xs,A,R,scale_fact=1.5)#generate 1.5*original number of triples
    print(sprintf('Precision:\%.4f, Recall:\%.4f',tmp$tp/(tmp$tp+tmp$fp),tmp$tp/(tmp$tp+tmp$fn)))

#using RESCAL for prediction missing relations.

aa=read.table(file = paste0(path.package("RDFTensor"), '/toy_vicePresident.nt'),
                                sep=' ',header=FALSE,stringsAsFactors=FALSE)
trp=aa[,1:3]

tnsr=getTensor(trp)
    r=4
    
    sr=NULL
    t0=proc.time()
    X_=list()
    library(Matrix)
    tt=rescal(tnsr$X,rnk=r,ainit='nvecs',verbose=1,lambdaA=0,epsilon=1e-4,lambdaR=0)
    R=tt$R
    A=tt$A
s1=A\%*\%R[[1]]\%*\%Matrix::t(A)
s2=A\%*\%R[[1]]\%*\%Matrix::t(A)
#predict the party of AlGore (no explicit info is given in the nt file)
print(s1[tnsr$SO=='<http://example.com/AlGore>',tnsr$SO=='<http://example.com/RepublicanParty>'])
print(s1[tnsr$SO=='<http://example.com/AlGore>',tnsr$SO=='<http://example.com/DemocraticParty>'])
partyOf=data.frame(tnsr$SO,Repub=s1[,tnsr$SO=='<http://example.com/RepublicanParty>'],
           Democ=s1[,tnsr$SO=='<http://example.com/DemocraticParty>'],
    GivenRepub=tnsr$X[[1]][,tnsr$SO=='<http://example.com/RepublicanParty>'],
          GivenDemoc=tnsr$X[[1]][,tnsr$SO=='<http://example.com/DemocraticParty>'],
    stringsAsFactors=FALSE)  
    
 print(partyOf)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RESCAL }% use one of  RShowDoc("KEYWORDS")

