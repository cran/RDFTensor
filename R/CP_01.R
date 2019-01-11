# 22/5/2018
# cp_nmu01: is thr applied to factors or reconstruction?

CP_01 <- function(X,P,pthr=c(0.001,0.01,0.05,0.1,0.2,0.3,0.4, 0.5,0.6,0.8),testAll=FALSE){
#P is the factorization as KTensor:list(lambda,u)
#X sptensor list(subs,vals,size,nnz)
# testAll:try all thr values or stop after having a zero in fp or tp
#1-try applying threshold to all Factors
#2-try applying threshold to multiplication
# return best thr,Factors or Reconstruction
 Res=NULL
 A=P$u[[1]]
 B=P$u[[2]]
 C=P$u[[3]]
 r=ncol(A)
 TP_plus_FN=nrow(X$subs)
 pthr=sort(pthr)#break when no more 1's
 # absorb lambda's
 # for(i in 1:r){	A[,i]=P$lambda[i]*A[,i] }
 # browser()
 A = A%*%diag(P$lambda^(1/3)); B = B%*%diag(P$lambda^(1/3)); C = C%*%diag(P$lambda^(1/3));
 # browser()
 bestSol=NULL
 finished=FALSE;
 minErr=Inf## minimum Error so far
 for(thr in pthr){
  print(thr)
	A01 = A>=thr
	B01 = B>=thr
	C01 = C>=thr
	
	sl=list()
	for(s in 1:nrow(B01))sl[[s]]=Matrix::spMatrix(x=0,j=1,i=1,nrow=nrow(A),ncol=nrow(C));
	nnz=0
	for(i in 1:r){# rank
	# ai[X]bi[X]ci
		print(i)
		# sl0=outer(xx$A[,i],xx$C[,i])
		if((sum(A01[,i])*sum(C01[,i]))>=(2*TP_plus_FN)){#Error more than 100%
		  print("Error>100%")
		  break;
		}
		tt=as.matrix(expand.grid(which(A01[,i]),which(C01[,i])))
		if(nrow(tt)==0) next;
		nnz=nnz+nrow(tt)#not precise
		for(s in which(B01[,i])) sl[[s]][tt] = 1#sl[[s]] + sl0
		gc()
	}	
	if(i<r) next;
	if(nnz>0 ){
		subs=getTnsrijk(sl)
		if(nrow(subs) > 2*TP_plus_FN){
			print('Error>nnz')
			next;
		}
		tp=sum(paste(X$subs[,1],X$subs[,2],X$subs[,3]) %in% paste(subs[,1],subs[,2],subs[,3]) )
		fp=nrow(subs)-tp
		fn=nrow(X$subs)-tp
		Res=rbind(Res,cbind(thr=thr,tp=tp,fn=fn,fp=fp,Error=fp+fn,nnz=nrow(subs),R=tp/(tp+fn),P=tp/(tp+fp)))
		print(sprintf("Thr:%.3f, TP:%d, FN:%d, FP:%d, Recall:%.3f, Precision:%.3f",thr,tp,fn,fp,tp/(tp+fn),tp/(tp+fp)))
		if((fp==0 || tp==0) && !testAll) finished=TRUE;
	}else{
			subs=NULL
			tp=0
			fn=nrow(X$subs)-tp
			fp=0
			Res=rbind(Res,cbind(thr=thr,tp=tp,fn=fn,fp=fp,Error=fp+fn,nnz=0,R=0,P=NA))
			break;
	}	
	Error=fn+fp
	if(Error<minErr) {
	print("-------min Error------")
		bestSol=list(u=list(A=A01,B=B01,C=C01),sl=sl,subs=subs,thr=thr,Error=Error,TP=tp,FP=fp,FN=fn)
		minErr=Error
	}
	if(finished) break;
  }
  return(list(Res,sol=bestSol))
}

# res=CP_01(X,P1[[1]])
# res$sol$u
# If |X| is the number of non-zeros in data X, we sampled 200 |X| locations of 0s, and computed the error
# the methods make for every non-zero and for the sampled zero elements. The sampled results were then 
# extrapolated to the size of the full tensor. We had to use sampling, as the factor matrices were too dense 
# to allow reconstructing the full tensor.
#######

CP_R01 <- function(X,P,pthr=c(1E-6,1e-4,0.001,0.01,0.05,0.1,0.2,0.3,0.4, 0.5,0.55,0.6,0.65,0.7,0.8),cntNnz=200,startSize=1e7){
#P is the factorization as KTensor:list(lambda,u)
#X sptensor list(subs,vals,size)
# if size of tensor < startSize all values are calculated
# Using sampling
# CP_01: converts A,B,C matrices to binary but here the result of the multiplication is converted to binary.
 Res=NULL
 A=P$u[[1]]
 B=P$u[[2]]
 C=P$u[[3]]
 r=ncol(A)
 TP_plus_FN=nrow(X$subs)
 pthr=sort(pthr)#break when no more 1's
 # absorb lambda's
 # for(i in 1:r){	A[,i]=P$lambda[i]*A[,i] }
 A = A%*%diag(P$lambda^(1/3)); B = B%*%diag(P$lambda^(1/3)); C = C%*%diag(P$lambda^(1/3));
 
 bestSol=NULL
 minErr=Inf## minimum Error so far
 if(startSize > prod(X$size)){
	print('checking all values...')
	grd=expand.grid(1:X$size[1],1:X$size[2],1:X$size[3])
	zi=grd[!(paste(grd[,1],grd[,2],grd[,3])%in% paste(X$subs[,1],X$subs[,2],X$subs[,3])),]
	# rzi=rowSums(A[zi[,1],]*B[zi[,2],]*C[zi[,3],])
	ssz=nrow(zi)
	sample_sz=0
	rzi=rowSums(A[zi[,1],]*B[zi[,2],]*C[zi[,3],])#what type of product is this!!
	}else{
	 sample_sz=min(cntNnz*TP_plus_FN,prod(X$size));#*(1+1/prod(X$size)):to coup with 1's in rand
	 print('sampling mode 1...')
	 V1=sample(1:X$size[1],sample_sz,replace=TRUE)
	 print('sampling mode 2...')
	 V2=sample(1:X$size[2],sample_sz,replace=TRUE)
	 print('sampling mode 3...')
	 V3=sample(1:X$size[3],sample_sz,replace=TRUE)
	 zi=cbind(V1,V2,V3)[!(paste(V1,V2,V3)%in% paste(X$subs[,1],X$subs[,2],X$subs[,3])),]
	 ssz=nrow(zi)
	 print(sprintf('sample_sz:%d, sszeros:%d',sample_sz,ssz))

	 print('rzi..')	
	 rzi=numeric(cntNnz*TP_plus_FN)
	 for(ch in 1:cntNnz){
		print(paste0("ch:",ch))
		ix=(1+(ch-1)*cntNnz):(ch*cntNnz)
		rzi[ix]=rowSums(A[zi[ix,1],]*B[zi[ix,2],]*C[zi[ix,3],])#what type of product is this!!
	
	 }
	}
 print('r1s..')
 r1s=rowSums(A[X$subs[,1],]*B[X$subs[,2],]*C[X$subs[,3],])
 for(thr in pthr){
	print(thr)
	tp=sum(r1s>=thr)
	fn=TP_plus_FN-tp
	sszErr=sum(rzi>=thr)
	#extrapolate
	fp=round(sszErr*(prod(X$size)-TP_plus_FN)/ssz);
	print(sprintf("tp:%d, sszErr:%d, fp(extrap.):%d",tp,sszErr,fp))
	Error=fp+fn
	Res=rbind(Res,cbind(thr=thr,tp=tp,fn=fn,fp=fp,sszErr=sszErr,ssz=ssz,sample_sz=sample_sz,Error=Error))
	if(fp==0 || tp==0) break;
 }
return(Res)
}

#######################################

CP_01ext <- function(X,P,pthr=c(1E-6,1e-4,0.001,0.01,0.05,0.1,0.2,0.3,0.4, 0.5,0.55,0.6,0.65,0.7,0.8),minR=0.01){
#P is the factorization as KTensor:list(lambda,u)
#X sptensor list(subs,vals,)
# if size of tensor < startSize all values are calculated
# calc. All
# return best thr,Factors or Reconstruction
 Res=NULL
 A=P$u[[1]]
 B=P$u[[2]]
 C=P$u[[3]]
 r=ncol(A)
 TP_plus_FN=nrow(X$subs)
 pthr=sort(pthr)#break when no more 1's
 # absorb lambda's
 A = A%*%diag(P$lambda^(1/3)); B = B%*%diag(P$lambda^(1/3)); C = C%*%diag(P$lambda^(1/3));
 
 bestSol=NULL
 minErr=Inf## minimum Error so far
 print('r1s..')
 r1s=rowSums(A[X$subs[,1],]*B[X$subs[,2],]*C[X$subs[,3],])
 ttp=NULL
 for(thr in pthr){
	print(thr)
	tp=sum(r1s>=thr)
	fn=TP_plus_FN-tp
	ttp=rbind(ttp,cbind(thr=thr,tp=tp,fn=fn))
	if(tp==0) break;
 }
 
 finished=FALSE
 pthr1=pthr[pthr %in% ttp[(ttp[,'tp']/TP_plus_FN)>=minR,'thr']];#remaining
 grd=expand.grid(1:X$size[1],1:X$size[2])
 AB=A[grd[,1],]*B[grd[,2],]
 print(sprintf('grd size:%d',nrow(grd)))
 for(i in 1:X$size[3]){
		print(i)
		t1=proc.time()
		# ix=cbind(grd[,1],i,grd[,2])
		# resi=rowSums(A[grd[,1],]*B[rep(i,nrow(grd)),]*C[grd[,2],])
		ix=which(C[i,]!=0)#in case of sparse factors
		resi=rowSums(AB[,ix,drop=FALSE]*C[rep(i,nrow(grd)),ix,drop=FALSE])
		# aa=cbind(grd[,1],i,grd[,2],resi)
		gc()
		for(ti in pthr1){
			resinnz=sum(resi>=ti)
			Res=rbind(Res,cbind(i=i,ti=ti,nnz=resinnz))
			ti_tp=ttp[ttp[,'thr']==ti,'tp']
			if(resinnz>2*ti_tp){
				print(sprintf("thr %f excluded, ti_tp:%d, resinnz:%d ",ti,ti_tp,resinnz))
				pthr1=pthr1[pthr1!=ti]#exclude ti
				if(length(pthr1)==0) finished=TRUE;
				break;
			}
		}
		if(i%%20==0){
				tnnz1=aggregate(Res[,'nnz'],list(Res[,'ti']),sum)
				th=which(tnnz1[,'x']>=2*TP_plus_FN)
				pthr1=pthr1[!(pthr1 %in% tnnz1[th,1])]
				print(pthr1)
				print(tnnz1[tnnz1[,1] %in% pthr1,])
				# print(tnnz1)
				t2=proc.time()
				print(t2-t1)
			}
 }

 tnnz=aggregate(Res[,'nnz'],list(Res[,'ti']),sum)
 rr=merge(ttp,tnnz,by.x='thr',by.y=1,all.x=TRUE)
 Error=rr[,'x']-rr[,'tp']+rr[,'fn']
 bestind=which.min(Error)
 
   return(list(sol=list(bestind,thr=rr[bestind,1],tp=rr[bestind,2],nnz=rr[bestind,'x']),rr,Res,ttp))
 }
 
 rescal_01<-function(X,A,R,scale_fact=1){
 # find best quantization of the resacl factorization A,R
    X_=list()
    sr=NULL
    for(s in 1:length(R)){
          t1n=A%*%R[[s]]%*%Matrix::t(A)
          qntl=scale_fact*sum(X[[s]])/(nrow(X[[s]])*ncol(X[[s]]))
          Xb=as(X[[s]],'TsparseMatrix')#Matrix::spMatrix(i=X[[s]]@i+1,j=X[[s]]@j+1,x=X[[s]]@x==1,nrow=nrow(X[[s]]),ncol=ncol(X[[s]]))
          thr=quantile(t1n@x,1-qntl)
         aa=which(t1n>thr,arr.ind=T)
         if(length(aa)>0){
            X_[[s]]=Matrix::spMatrix(i=aa[,1],j=aa[,2],x=rep(1,nrow(aa)),nrow=nrow(A),ncol=nrow(A))#tt > threshold[i],'sparseMatrix')
        }else{
           X_[[s]]=Matrix::spMatrix(i=1,j=1,x=0,nrow=nrow(A),ncol=nrow(A))
        }
        li=Xb@i[Xb@x==1]+1
        lj=Xb@j[Xb@x==1]+1
        tp=sum(X_[[s]][cbind(li,lj)])
        fn=sum(Xb@x)-tp#sum(!X_[cbind(li,lj)])
        fp=sum(X_[[s]]@x)-tp# incase of scale_fact=1 fp=fn as number of 1's in X_ and X is the same
        # print(sprintf('s,r,nnz,tp,fn,fp:%d %d %d %d %d %d %.2f %.2f',s,r,sum(Xb),tp,fn,fp,R=tp/(tp+fn),P=tp/(tp+fp)))
        sr=rbind(sr,cbind(s=s,scale_fact=scale_fact,thr=thr,nnz=sum(Xb),tp=tp,fn=fn,fp=fp,R=tp/(tp+fn),P=tp/(tp+fp)))
    }
    tp=sum(sr[,'tp'])
    fn=sum(sr[,'fn'])
    fp=sum(sr[,'fp'])
    return(list(X_,tp=tp,fn=fn,fp=fp,sr=sr)) 
 }

 ### 22/11/2018

getCP_val<-function(P,ijk){
#P is the factorization as KTensor:list(lambda,u)
#ijk indexes on the 3-mode tensor
# get the values corresponding to triples from the factorization

 A=P$u[[1]]
 B=P$u[[2]]
 C=P$u[[3]]
 
 # absorb lambda's
 A = A%*%diag(P$lambda^(1/3)); B = B%*%diag(P$lambda^(1/3)); C = C%*%diag(P$lambda^(1/3));
 r1s=rowSums(A[ijk[,1],]*B[ijk[,2],]*C[ijk[,3],])
 return(cbind(ijk,r1s))
}


absorb_lambdas<-function(P){
#P is the factorization as KTensor:list(lambda,u)
# make lambda 1's

 A=P$u[[1]]
 B=P$u[[2]]
 C=P$u[[3]]
 
 # absorb lambda's
    A = A%*%diag(P$lambda^(1/3)); B = B%*%diag(P$lambda^(1/3)); C = C%*%diag(P$lambda^(1/3));
    
    return(list(lambda=rep(1,3),u=list(A,B,C)))
 }
 
 CP_get_Frontal_slice<-function(P,s=1){
   return(P$u[[1]]%*%diag(P$u[[3]][s,,drop=TRUE])%*%t(P$u[[2]]))
 }