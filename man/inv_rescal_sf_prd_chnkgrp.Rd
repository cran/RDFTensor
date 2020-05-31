\name{inv_rescal_sf_prd_chnkgrp}
\alias{inv_rescal_sf_prd_chnkgrp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reconstruct predicate(s) via chunks and grouping
}
\description{reconstruct a predicate or set of predicates from RESCAL factorization result
by returning top scores in each predicate.
}
\usage{
inv_rescal_sf_prd_chnkgrp(R, A, cnt_prd, scale_fact = 1, verbose = 1, ncores = 3,
ChnkLen = 1000, saveRes = FALSE, OS_WIN = FALSE, pve = 1e-10, grpLen = NULL,
dsname = "", rmx = NULL, readjustChnkLen = TRUE, TotalChnkSize = 5e+08, 
chTrpCntTol = 1000, generateLog = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{R}{
  core tensor resulting from RESCAL factorization (r by r by m ).
}
\item{A}{
    Embedding matrix part resulting from RESCAL factorization.
}
\item{cnt_prd}{
vector with length of number of predicate giving number of triples in each predicate. 
When 0 the predicate will not be processed.
}

  \item{scale_fact}{
scale factor, the generated number of triples in a predicate s is sf*prd_cnt[s]
}
  \item{verbose}{
 the level of messages to be displayed, 0 is minimal.
}
\item{ncores}{ number of cores used to run in parallel, 0 means no paralellism
}
\item{OS_WIN }{ True when the operating system is windows, 
used to allow using Fork when running in parallel
}
  \item{pve}{
 positive value: representing the smallest value of allowed score of reconstructed triple.
}
  \item{grpLen}{length of one group of iterations, 
  when running iterations in parallel results are collected for all iterations to be summarized after last iteration. 
  Thus more memory is required. To avoid that iterations are divided to groups with summaries calculated for each group. 
  Default 15.
}
  \item{ChnkLen}{
number of rows in one chunk. Default 1000.
}
  \item{generateLog}{
  save output when running in parallel to a log file in current directory.
}

  \item{saveRes}{
  opionally save each predicate
}

  \item{dsname}{
  optional:name of dataset
}
  \item{rmx}{
  optionally give the max absolute value in A to save its calculation in multiple calls
}
  \item{readjustChnkLen}{
  automatically increase chunk length when possible
}
  \item{TotalChnkSize}{
  instead of defining ChnkLen define the number of pairs of entities to processed
  every chunk. Equal to ChnkLen*N, where N is the number of entities (i.e. nrows(A))
}
  \item{chTrpCntTol}{
  tolerance in number of triples returned in one chunk, will be eliminated at end of group
}
}
\details{
Multiplication of A*R[[p]]*A^T can be impossible in typical 16GB RAM machine when the number of entities is more than 50K.
(needs ~25GB RAM) To deal with that we use chunks to obtain the top scores from such multiplication.
The main idea is to get the required number of triples from each chunk then choose from them again the top scores.

}
\value{
LIST of three components:
    \item{ijk}{A data frame containing the reconstructed triples using indexes of entities}
    \item{val}{A vector containing the scores of the reconstructed triples}
    \item{act_thr}{the minimum score in each predicate (minimum score of a triple, threshold)}
  }
\references{
    -Maximilian Nickel, Volker Tresp, Hans-Peter-Kriegel,
      "Factorizing YAGO: Scalable Machine Learning for Linked Data"
      WWW 2012, Lyon, France
      
    -SynthG: mimicking RDF Graphs Using Tensor Factorization, Desouki et al
}
\author{
Abdelmoneim Amer Desouki
}

\seealso{
    \code{\link{rescal}} 
    \code{\link{scRescal}} 
}

\examples{
\dontrun{
  library(RDFTensor)
  data('umls_tnsr')
   ntnsr=umls_tnsr
    
   #Calculate Factorization
    tt0=proc.time()
    tt=rescal(ntnsr$X,rnk=10,ainit='nvecs',verbose=2,lambdaA=0,epsilon=1e-4,lambdaR=0)
    ttq1=proc.time()
    A=tt$A
    R=tt$R
    
    # reconstruct second predicate (slice) in tensor

  p=2
  prd_cnt=rep(0,length(ntnsr$X))#Zero counts will not be reconstructed
  prd_cnt[p]=sum(ntnsr$X[[p]])
  tmpRes=inv_rescal_sf_prd_chnkgrp(R,A,cnt_prd=prd_cnt,ChnkLen=50,grpLen=6,OS_WIN=TRUE,ncores=1,
            chTrpCntTol=1000,  TotalChnkSize=1e4)

    ijk=tmpRes[[1]]
    ix=which(ntnsr$X[[p]]==1,arr.ind=TRUE)
    oijk=cbind(ix[,1],p,ix[,2])#Original
    flag= paste(oijk[,1],oijk[,2],oijk[,3]) \%in\% paste(ijk[,1],ijk[,2],ijk[,3])
    print(table( flag))#True positives
  
    pTrp=cbind(ntnsr$SO[ijk[,1]],ntnsr$P[ijk[,2]],ntnsr$SO[ijk[,3]])
    
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RESCAL }% use one of  RShowDoc("KEYWORDS")
