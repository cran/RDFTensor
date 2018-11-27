\name{RDFTensor-package}
\alias{RDFTensor-package}
\alias{RDFTensor}
\encoding{utf8}
\docType{package}
\title{
\packageTitle{RDFTensor}
}
\description{
\packageDescription{RDFTensor}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{RDFTensor}
\packageIndices{RDFTensor}
Read the tensor using getTensor or use NTriples file and parseNT.
Use one of the methods to factorize it cp_nmu, cp_apr, parCube or rescal.
}
\author{Abdelmoneim Amer Desouki}

\references{
   -Brett W. Bader, Tamara G. Kolda and others. 
     MATLAB Tensor Toolbox, Version [v3.0]. Available online at https://www.tensortoolbox.org, 2015.

   -Papalexakis, Evangelos E., Christos Faloutsos, and Nicholas D. Sidiropoulos. 
   "Parcube: Sparse parallelizable tensor decompositions." Machine Learning and Knowledge Discovery in Databases. 
   Springer Berlin Heidelberg, 2012. 521-536.

   NMU
   -Lee, Daniel D., and H. Sebastian Seung. "Algorithms for non-negative matrix factorization." 
   In Advances in neural information processing systems, pp. 556-562. 2001.
   
   -Anh Huy Phan, Petr Tichavský, Andrzej Cichocki, On Fast Computation of Gradients for 
   CANDECOMP/PARAFAC Algorithms, arXiv:1204.1586, 2012
   
   CP_APR
   -E. C. Chi and T. G. Kolda. On Tensors, Sparsity, and Nonnegative
    Factorizations, SIAM J. Matrix Analysis,  33(4):1272-1299, Dec. 2012,
    http://dx.doi.org/10.1137/110859063  
   
   -S. Hansen, T. Plantenga and T. G. Kolda, Newton-Based Optimization
    for Kullback-Leibler Nonnegative Tensor Factorizations, 
    Optimization Methods and Software, 2015, 
    http://dx.doi.org/10.1080/10556788.2015.1009977
    
    RESCAL
   -Nickel, Maximilian, and Volker Tresp. "Tensor factorization for multi-relational learning."
    In Joint European Conference on Machine Learning and Knowledge Discovery in Databases, 
    pp. 617-621. Springer, Berlin, Heidelberg, 2013.
    
}
\keyword{ Tensor Factorization }
\keyword{ Tensor Decomposition }
\keyword{ RDF Tensor }
\seealso{
%~~ Optional links to other man pages, e.g. ~~
%\code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
    \code{\link{cp_nmu}} 
    \code{\link{cp_apr}} 
    \code{\link{serial_parCube}} 
    \code{\link{rescal}}
    \code{\link{cp_als}}     
}
\examples{
%~~ simple examples of the most important functions ~~
\dontrun{
 data(umls_tnsr)
 ntnsr=umls_tnsr
    r=10
    
    sr=NULL
    t0=proc.time()
    X_=list()
    library(Matrix)
    tt=rescal(umls_tnsr$X,rnk=10,ainit='nvecs',verbose=1,lambdaA=0,epsilon=1e-2,lambdaR=0)
    R=tt$R
    A=tt$A
     for(s in 1:length(R)){
      print(sprintf('-------------=======s=\%d==========---------',s))
      t1=proc.time()
      t1n=A\%*\%R[[s]]\%*\%Matrix::t(A)
      mx=max(t1n@x)
      Xb=Matrix::spMatrix(i=ntnsr$X[[s]]@i+1,j=ntnsr$X[[s]]@j+1,x=ntnsr$X[[s]]@x==1,
      nrow=nrow(ntnsr$X[[s]]),ncol=ncol(ntnsr$X[[s]]))
      
      all_scale_fact=c(0.7,0.8,0.9,0.95,1,1.05,1.1,1.2,1.3,1.4,1.5,1.8)
      for(scale_fact in all_scale_fact){
          qntl=scale_fact*sum(Xb)/(nrow(ntnsr$X[[s]])*ncol(ntnsr$X[[s]]))
            thr=quantile(t1n@x,1-qntl)
            print(sprintf('-------------\%f---------',thr))
                 aa=which(t1n>thr,arr.ind=TRUE)
                 if(length(aa)>0){
                    X_[[s]]=Matrix::spMatrix(i=aa[,1],j=aa[,2],x=rep(1,nrow(aa)),
                    nrow=nrow(A),ncol=nrow(A))#tt > threshold[i],'sparseMatrix')
                }else{
                   X_[[s]]=Matrix::spMatrix(i=1,j=1,x=0,nrow=nrow(A),ncol=nrow(A))
                }
        #---
            li=Xb@i[Xb@x]+1
            lj=Xb@j[Xb@x]+1
            tp=sum(X_[[s]][cbind(li,lj)])
            fn=sum(Xb@x)-tp#sum(!X_[cbind(li,lj)])
            # incase of scale_fact=1 fp=fn as number of 1's in X_ and X is the same
            fp=sum(X_[[s]]@x)-tp
            sr=rbind(sr,cbind(s=s,r=r,scale_fact=scale_fact,mx=mx,thr=thr,nnz=sum(Xb),
            tp=tp,fn=fn,fp=fp,R=tp/(tp+fn),P=tp/(tp+fp)))
        #    if(tp==0) break;
        }
        t2=proc.time()
        print(t2-t1)  
    }
 tf=proc.time()
 print(tf-t0)

 Res=NULL
    for(sf in all_scale_fact){
        sr.sf=sr[sr[,'scale_fact']==sf ,]
         R=sum(sr.sf[,'tp'])/(sum(sr.sf[,'tp'])+sum(sr.sf[,'fn']))
         P=sum(sr.sf[,'tp'])/(sum(sr.sf[,'tp'])+sum(sr.sf[,'fp']))
        cnt=nrow(sr.sf)
        Res=rbind(Res,cbind(sf=sf,P,R,cnt))
    }
   print(Res)

   stats=Res

     plot(stats[,'sf'],stats[,'R']*100,type='b',col='red',lwd=2,
    main=sprintf('RESCAL, choosing scale factor (sf):(ntrp*sf), dataset: \%s, 
        #Slices:\%d\n #Known facts:\%d','UMLS',length(ntnsr$X),
        sum(sr.sf[,'tp']+sr.sf[,'fn'])),ylab="",xlab='Scale Factor',
    xlim=c(0,max(sf)),ylim=c(0,100))
    HM=apply(stats,1,function(x){2/(1/x['P']+1/x['R'])})
    points(stats[,'sf'],stats[,'P']*100,col='blue',lwd=2,type='b')
    points(stats[,'sf'],100*HM,col='green',lwd=2,type='b')
    grid(nx=10, lty = "dotted", lwd = 2)
    legend(legend=c('Recall','Precision','Harmonic mean'),col=c('red','blue','green'),
    x=0.6,y=20,pch=1,cex=0.75,lwd=2)

    max(HM)


    hist(sr[sr[,'scale_fact']==1,'thr'],col='grey',
    main='UMLS Reconstring the same number of triples, Actual threshold',
    xlab='threshold',cex.main=0.75)
}


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
    X=list(subs=subs,vals=rep(1,nrow(subs)),size=c(5,2,5))
    normX=sqrt(sum(X$vals))
    set.seed(123)
    # NMU decomposition with rank 2
    P1=cp_nmu(X,2)
   
    ###
    # find best CP boolean Factorization based on NMU
    res=CP_01(X,P1[[1]])
    Fac=res$sol$u # The factorization
    # TP,FP,FN
    print(sprintf("TP=\%d, FP=\%d, FN=\%d, Threshold=\%f",res$sol$TP,res$sol$FP,res$sol$FN,res$sol$thr))

    #################   CP_APR   ###################
    res=cp_apr(X,R=2,opts=list(alg='pdnr',printinneritn=1))

    set.seed(12345)
    res1=CP_01(X,res[[1]])
    res2=CP_R01(X,res[[1]])
    #res3=CP_01ext(X,res[[1]])
    ###################   RESCAL  #############################

    tt=rescal(tnsr$X,2,ainit='random',verbose=3)
    R=tt[['R']]
    A=tt[['A']]

    Tensor_error(tnsr$X,A,R)
    t1n= A%*%R[[1]]%*%Matrix::t(A)
    t2n= A%*%R[[2]]%*%Matrix::t(A)
}
