\name{tnsr2trp}
\alias{tnsr2trp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
sparse tensor to triples
}
\description{
convert a sparse tensor (e.g. created using \code{\link{getTensor}}) to set of triples.
}
\usage{
tnsr2trp(tnsr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tnsr}{
a sparse tensor (e.g. created using \code{\link{getTensor}})
}
}
\value{
a matrix of three columns containing the triples
}
\author{
Abdelmoneim Amer Desouki}

\seealso{
 \code{\link{getTnsrijk}}
 \code{\link{getTensor}}
}
\examples{
trp=rbind(
    cbind('Alex',  'loves', 'Don'),
    cbind('Alex',  'loves', 'Elly'),
    cbind('Alex',  'hates', 'Bob'),
    cbind('Don',   'loves', 'Alex'),
    cbind('Don',   'hates', 'Chris'),
    cbind('Chris', 'hates', 'Bob'),
    cbind('Bob',   'hates', 'Chris'),
    cbind('Elly',  'hates', 'Chris'),
    cbind('Elly',  'hates', 'Bob'),
    cbind('Elly',  'loves', 'Alex')
    )
######
# form tensor as a set of frontal slices (Predicate mode)
    tnsr=getTensor(trp)
    trp_=tnsr2trp(tnsr)
    #print(any(! paste(trp_[,1],trp_[,2],trp_[,3]) %in% paste(trp[,1],trp[,2],trp[,3])))
    #print(any(! paste(trp[,1],trp[,2],trp[,3]) %in% paste(trp_[,1],trp_[,2],trp_[,3])))
}