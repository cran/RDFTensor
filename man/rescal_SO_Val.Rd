\name{rescal_SO_Val}
\alias{rescal_SO_Val}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calaculate scores Subj,Predicate,Obj
}
\description{
calaculate scores from RESCAL factorization for a pairse of subjects and objects indices
on one predicate
}
\usage{
rescal_SO_Val(R, A, Subj, P, Obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{
    core tensor resulting from RESCAL factorization (r by r by m ).
}
  \item{A}{
    Embedding matrix part resulting from RESCAL factorization.
}
  \item{Subj}{
integer vector containing indices of Subjects
}
  \item{P}{
  index of predicate  (slice in tensor)
}
  \item{Obj}{
  integer vector containing inedices of Objects
}
}
\details{calaculate scores from RESCAL factorization for a pairse of subjects and objects indices
on one predicate. Subj and Obj parameters must have the same length.
rowSums(A[Subj,] \%*\% R[[p]] * A[Obj,])
}
\value{
Data frame of four columns: Subj, Predicate, Obj, and val.
val column is the score
}
\author{
Abdelmoneim Amer Desouki
}

\seealso{
    \code{\link{rescal}} 
    \code{\link{scRescal}} 
    \code{\link{rescal_Trp_Val}}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (R, A, Subj, P, Obj) 
{
    if (length(Subj) != length(Obj)) {
        stop("rescal_SO_Val: List of subjects and objects must be of the same length")
    }
    val = rowSums((A[Subj, , drop = FALSE] \%*\% R[[P]]) * (A[Obj, 
        , drop = FALSE]))
    return(data.frame(Subj, P, Obj, val))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RESCAL }% use one of  RShowDoc("KEYWORDS")

