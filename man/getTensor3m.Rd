\name{getTensor3m}
\alias{getTensor3m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    getting tensor from triples
}
\description{
   gets a tensor data structure as a set of frontal slices represented as sparse matrices.
   The entities as subject has a different id than object and the size of 
   subject and object modes can be different.
}
\usage{
getTensor3m(trp, S = NULL, P = NULL, O = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trp}{
    three columns representing the set of triples (subject, predicate, object).
    Each row representes one triple.
}
  \item{S}{
    the predifined list of subjects
}
  \item{P}{
    The set of predicates
}
 \item{O}{
    the predifined list of objects
}
}
\value{
  A LIST containing:
  \item{X }{set of sparse matrices with ones in the position of triples,
  each sparse matix represents one predicate (relationship)}
  \item{S }{the set of entities in the subject domain}
  \item{O }{the set of entities in the object domain}
  \item{P }{the set of predicates}
}

\author{
Abdelmoneim Amer Desouki
}
