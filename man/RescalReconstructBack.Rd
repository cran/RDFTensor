\name{RescalReconstructBack}
\alias{RescalReconstructBack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
reconstruct a tensor from RESCAL factorization result
}
\description{
 reconstruct a tensor from RESCAL factorization result
}
\usage{
RescalReconstructBack(A, R, trpo = NULL, otnsr = NULL, chkZR = FALSE, prd_cnt = 
NULL, sf = 1, verbose = 1, ncores = 3, OS_WIN = FALSE, pve = 1e-10, 
grpLen = NULL, ChnkLen = 1000, generateLog = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{A}{
    Embedding matrix part resulting from RESCAL factorization.
}
\item{R}{
  core tensor resulting from RESCAL factorization (r by r by m ).
}
\item{trpo}{
optional original tensor in for of triples(data frame of three columns subject, predicate, object)
}
  \item{otnsr}{
sparse tensor containing the triples to be reconstructed (i.e. to calculate scores for).
}
  \item{chkZR}{
flag to check if there is any slices in core tensor which are entirely zeros
}
\item{prd_cnt}{
optional. Giving number of triples in each predicate.
}
  \item{sf}{
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
}
\details{
The function finds triples with top scores from RESCAL factorization.
When the scale factor is 1 each predicate the resulting number of triples will be equal to that of the original tensor.
Each predicate will get the same number of triples as that of the original graph.
Since direct multiplication of A by A^T may not fit into RAM even for mid-size graphs, chunks are used to get the results.
}
\value{
The result is a LIST of three items:
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
    
    # reconstruct tensor and calulate
    RecRes=RescalReconstructBack(R=R,A=A,otnsr=ntnsr,ncore=3,OS_WIN=TRUE)
    print(sprintf('True positive rate:\%.2f \%\%',100*sum(RecRes$TP)/length(RecRes$TP)))
    
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RESCAL }% use one of  RShowDoc("KEYWORDS")
