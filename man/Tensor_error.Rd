\name{Tensor_error}
\alias{Tensor_error}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RESCAL Tensor error
}
\description{
Calculates error (Eculidean distance) between X:sparse tensor (frontal slices) 
and the approximations as R :core tensor and matrix :A.
}
\usage{
Tensor_error(X, A, R)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{

}
  \item{A}{
matrix A from RESCAL factorization (n by r).
}
  \item{R}{
core tensor from RESCAL factorization as set of frontal slices ( r by r by r)
}
}
\value{
return numeric value as the Eculidean distance between the two tensors.
}

\author{
Abdelmoneim Amer Desouki}

\seealso{
    \code{\link{rescal}} 
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

}
