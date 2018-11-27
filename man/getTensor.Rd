\name{getTensor}
\alias{getTensor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    getting tensor from triples
}
\description{
   gets a tensor data structure as a set of frontal slices represented as sparse matrices.
   The entities have the same id as subject and/or object.
}
\usage{
getTensor(trp, SO = NULL, P = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trp}{
    three columns representing the set of triples (subject, predicate, object).
    Each row representes one triple.
}
  \item{SO}{
    the predifined list of entities
}
  \item{P}{
    The set of predicates
}
}
\value{
  A LIST containing:
  \item{X }{set of sparse matrices with ones in the position of triples,
  each sparse matix represents one predicate (relationship)}
  \item{SO }{the set of entities in the subject and object domains}
  \item{P }{the set of predicates}
}

\author{
Abdelmoneim Amer Desouki
}
