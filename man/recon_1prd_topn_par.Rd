\name{recon_1prd_topn_par}
\alias{recon_1prd_topn_par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Reconstruct predicate top scores on very large graphs
}
\description{
Reconstruct predicate top scores on very large graphs from RESCAL Factorization A & R.
NB: A and R should be already loaded in the environment calling the function.
Calculates top scores of A\%*\%R[[p]]\%*\%A^T.
Uses chunks for rows and columns and constraints on maximum possible value of scores 
to avoid calculations of too small values
}
\usage{
recon_1prd_topn_par(A,R,p, pcnt, rchLen = 1000, cchLen = 200, pve = 1e-10, 
mxrIter = 5, mxcIter = 25, grpLen = 40, OS_WIN = FALSE, dsname = "", ncores = 8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{A}{
    Embedding matrix part resulting from RESCAL factorization.
}
\item{R}{
  core tensor resulting from RESCAL factorization (r by r by m ).
}
  \item{p}{
predicate number.
}
  \item{pcnt}{
Number of triples in the predicate (to be reconstructed)
}
  \item{rchLen}{
row chunk length.
}
  \item{cchLen}{
  column chunk length
}
  \item{pve}{
 positive value: representing the smallest value of allowed score of reconstructed triple.
}
  \item{mxrIter}{
maximum number of iterations (chunks) in rows}
  \item{mxcIter}{
maximum number of iterations (chunks) in columns
}
 \item{grpLen}{length of one group of iterations, 
  when running iterations in parallel results are collected for all iterations to be summarized after last iteration. 
  Thus more memory is required. To avoid that iterations are divided to groups with summaries calculated for each group. 
  Default 40.
}
\item{ncores}{ number of cores used to run in parallel, 0 means no paralellism
}
\item{OS_WIN }{ True when the operating system is windows, 
used to allow using Fork when running in parallel
}
  \item{dsname}{
  optional:name of dataset
}
}

\value{
The result is a LIST of three items:
    \item{ikv}{A data frame containing the reconstructed triples (subject, Object, value, rchnk,chnk) 
    using indexes of entities, rchnk and colun chunk are the chunk in which the triple is generated
    note predicate is an argument}
    \item{minThr}{the minimum score in each predicate (minimum score of a triple, threshold)}
    \item{Iter}{Number of iterations done}
  }

\references{
    -Maximilian Nickel, Volker Tresp, Hans-Peter-Kriegel,
      "Factorizing YAGO: Scalable Machine Learning for Linked Data"
      WWW 2012, Lyon, France
      
    -SynthG: mimicking RDF Graphs Using Tensor Factorization,
	 Desouki et al.  IEEE ICSC 2021
}
\author{
Abdelmoneim Amer Desouki
}

\seealso{
    \code{\link{rescal}} 
    \code{\link{scRescal}} 
    \code{\link{RescalReconstructBack}} 
    \code{\link{inv_rescal_sf_prd_chnkgrp}} 
}
\examples{
\dontrun{
#Reconstructing one predicate from DBpedia factorization
#print(load('dbpi42_r100.RData'))
lst=recon_1prd_topn_par(A=A,R=R,p,pcnt=prd_cnt1[p],rchLen=250000,cchLen=750,
                     mxrIter=30,mxcIter=70,ncores=12,grpLen=40,pve=pve,dsname=name)
    Res=cbind(S=lst$ikv[,1],P=p,lst$ikv[,2:5])#triples , value
 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RESCAL }% use one of  RShowDoc("KEYWORDS")
\keyword{ SynthG }% __ONLY ONE__ keyword per line
