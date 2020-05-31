\name{rescal_Trp_Val}
\alias{rescal_Trp_Val}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    RESCAL scores of triples 
}
\description{
   Calculate score values corresponding to triples in tensor from RESCAL factorization.
}
\usage{
rescal_Trp_Val(R, A, ntnsr, Sl_ix = NULL, verbose = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{
    core tensor resulting from RESCAL factorization (r by r by m ).
}
  \item{A}{
    Embedding matrix part resulting from RESCAL factorization.
}
  \item{ntnsr}{
    sparse tensor containing the triples to be reconstructed (i.e. to calculate scores for).
}
  \item{Sl_ix}{
    Optional: indexes of slices to be processed.
}
  \item{verbose}{
 the level of messages to be displayed, 0 is minimal.
}
}
\value{
  Returns a data frame containing four columns:
\item{Subj}{index(es) of subject(s)}
\item{P}{index(es) of predicate}
\item{Obj}{index(es) of object(s)}
\item{val}{the score resuling from A[Subj,]\%*\% R_p * A[Obj,].}

}

\author{
Abdelmoneim Amer Desouki
}

\seealso{
    \code{\link{rescal}} 
    \code{\link{scRescal}} 
     \code{\link{RescalReconstructBack}}
}

\examples{
\dontrun{
  library(RDFTensor)
  data('umls_tnsr')
  ntnsr=umls_tnsr
   
   #Calculate RECAL Factorization
    tt=rescal(ntnsr$X,rnk=10,ainit='nvecs',verbose=2,lambdaA=0,epsilon=1e-4,lambdaR=0)
    A=tt$A
    R=tt$R
    
    #Calculate scores
    res=rescal_Trp_Val(tt$R,tt$A,ntnsr)
    
    % plot(density(res[,'val']),main='RESCAL Factorization rank=10, density of triples of UMLS')
    print(summary(res[,'val']))
}
}