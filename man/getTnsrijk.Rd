\name{getTnsrijk}
\alias{getTnsrijk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 tensor frontal slices to indices
}
\description{
convert tensor frontal slices to indices with three columns (i.e i ,j , k)
}
\usage{
getTnsrijk(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
LIST of sparce matrices, with each entry representing one predicate.
}
}
\value{
a matrix of three columns representing indices of 1 values in the tensor
}
\author{
Abdelmoneim Amer Desouki}

\seealso{
 \code{\link{tnsr2trp}}
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
    subs=getTnsrijk(tnsr$X)
    print(subs)
}
