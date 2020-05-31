# Mar-2020
# Reconstruct predicate top scores on very large graphs from RESCAL Factorization A & R.

#  Copyright (C) 2020 Abdelmoneim Amer Desouki, 
#   Data Science Group, Paderborn University, Germany.
#  All right reserved.
#  Email: desouki@mail.upb.de
#
#  This file is part of RDFTensor.
#
#  RDFTensor is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  RDFTensor is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with RDFTensor.  If not, see <http://www.gnu.org/licenses/>.

###----------------------
recon_1prd_topn_par<-function(A,R,p,pcnt,rchLen=1000,cchLen=200,pve=1e-10,mxrIter=5,mxcIter=25,grpLen=40,OS_WIN=FALSE,dsname='',ncores=8){
  
  # Reconstruct predicate top scores on very large graphs from RESCAL Factorization A & R.
  # Calculates top scores of A\%*\%R[[p]]\%*\%A^T.
  # Uses chunks for rows and columns and constraints on maximum possible value of scores 
  # to avoid calculations of too small values
  # 
 inv_rescal_calc_chnk<-function(rch,cch){
    t0=proc.time()
    cIx=cSorted[(1+(cch-1)*cchLen):min(cch*cchLen,N)]
    rIx=rSorted[(1+(rch-1)*rchLen):min(rch*rchLen,N)]
    
    tpv=AR[rIx,,drop=FALSE]%*%t(A[cIx,,drop=FALSE])
    tq0=proc.time()
    qntl=1-pcnt*1.0/prod(dim(tpv))
    # print(paste(qntl,dim(tpv)))
    if(qntl>0){
      thr=quantile(tpv,qntl)
    }else{
      thr=min(tpv)
    }
    tq1=proc.time()
    thr=max(minThr,thr)
    ik=which(tpv>=thr,arr.ind=TRUE)
    if(nrow(ik)>0){
        Val=tpv[ik]
        ikv=cbind(S=rIx[ik[,1]],O=cIx[ik[,2]],val=Val,rch,cch)  
        t1=proc.time()
        print(sprintf("p=%d, rch=%d, cch=%d, minThr=%f, qntltime=%.2f, ttime=%.2f",p,rch,cch,minThr,(tq1-tq0)[3],(t1-t0)[3]))
        
        return(list(ikv=ikv,minThr=thr))
    }else{
           print(sprintf("p=%d, rch=%d, cch=%d, no values",p,rch,cch))
           # minThr=min(ikv[,'val'])
           return(list(ikv=NULL,minThr=thr))
    }
  }
 
    #require(parallel)
	#require(doParallel)

    if(OS_WIN){
               cluster <- parallel::makeCluster(ncores,outfile=sprintf("%s_recon_1prd_topn_par_nc%d.log",dsname,ncores))
        }else{ 
               cluster <- parallel::makeForkCluster(ncores,outfile=sprintf("%s_recon_1prd_topn_par_nc%d.log",dsname,ncores))#not in windows, copy only changed
    }
	registerDoParallel(cluster)
  
        print("calc max rsA:")
        rsA=rowSums(abs(A))
        mxrsA=max(rsA)
    
    N=nrow(A)
    ts0=proc.time()
    t1=proc.time()
    print(sprintf('calc RAt, p=%d',p))
    RAt=as.matrix(R[[p]])%*%t(A)
    print(sprintf('calc rowMax..RAt, p=%d',p))
    mxRAtc=apply(abs(RAt),2,max)##
    rm(RAt)
    # mxmxRAtc=max(mxRAtc)*mxrsA;  if(mxmxRAtc < pve){  print(mxmxRAtc) }
    print(sprintf('calc AR, p=%d',p))
    AR = A%*% as.matrix(R[[p]])
    print(sprintf('calc rowMax..AR, p=%d',p))
    mxARr=apply(abs(AR),1,max)
    t2=proc.time()
    print(t2-t1)      
      
    cSorted=(1:N)[order(-mxRAtc)]
    rSorted=(1:N)[order(-mxARr)]
    cDone=rep(FALSE,N)
    rDone=rep(FALSE,N)
    
    ikv=NULL    
    minThr=0
    # nonTestedCThr=max(mxRAtc[!cDone])*mxrsA
    # nonTestedRThr=nonTestedCThr=1
    # Iter=1
    # grpCh=expand.grid(rch=1:mxrIter,cch=1:min(mxrIter,mxcIter))
    #p1 top 5M rows, top 30k cols mxrIter=20, mxcIter=60, Iter =ifelse(ix<mxrIter^2,floor(sqrt(ix)),(ix-mxrIter)/mxrIter 1:mxrIter+incr 1 every mxrIter
    Tbl=cbind(Iter=1,rch=1,cch=1)
    mxrIter=min(mxrIter,ceiling(N/rchLen))
    mxcIter=min(mxcIter,ceiling(N/cchLen))
    for(Iter in 2:mxrIter){#mxrIter<mxcIter
      for(rch in 1:(Iter-1) ) Tbl=rbind(Tbl,cbind(Iter,rch=rch,cch=Iter))
      for(cch in 1:Iter ) Tbl=rbind(Tbl,cbind(Iter,rch=Iter,cch=cch))
    }
    # browser()
    Iter=mxrIter
    if(mxrIter<mxcIter){
        for (cch in (Iter+1):mxcIter){
          Iter=Iter+1
          Tbl=rbind(Tbl,cbind(Iter,rch=1:mxrIter,cch=cch))
        }
    }        
    grpLen=min(nrow(Tbl),grpLen)
    grpCnt=ceiling(nrow(Tbl)/grpLen)
   
    # chkIx=(1:grpCnt) * grpLen
    # chkIter=ifelse(Tbl[chkIx,'Iter'] < mxrIter,Tbl[chkIx,'Iter']-1,Tbl[chkIx,'Iter'])#Tbl[,'cch'])
    for(g in 1:grpCnt ){
       print(paste("Group: ",g))
       t0=proc.time()
       opts <- list(preschedule=FALSE)#not using roundrobin instead when ended pick one from queue 
       # tmp=foreach(i=1+(g-1)*grpLen:min(grpLen*g,nrow(Tbl)),.packages='Matrix',.options.multicore=opts) %dopar%{#Lnx
       tmp=foreach(i=(1+(g-1)*grpLen):min(grpLen*g,nrow(Tbl)),.packages='Matrix',.export=c("A","R"),.options.multicore=opts) %dopar%{
       # for(i in 1+(g-1)*grpLen:min(grpLen*g,nrow(Tbl))) {
         tmpR=inv_rescal_calc_chnk(rch=Tbl[i,'rch'],cch=Tbl[i,'cch'])
       }
          
          # browser()
          # minThr,ikv
          tge0=proc.time()
  
          for(i in 1:length(tmp)){
            if(is.null(tmp[[i]]$ikv)) next;
            ikv=rbind(ikv,tmp[[i]]$ikv[tmp[[i]]$ikv[,'val']>minThr,,drop=FALSE])
            if(nrow(ikv) > pcnt){
              print(paste("resorting ... ",i))
              ik_order = order(ikv[,'val'],decreasing=TRUE)
              ikv = ikv[ik_order,,drop=FALSE][1:pcnt,,drop=FALSE]
            }
            minThr=min(ikv[,'val'])
           }        
          
          tge1=proc.time()
         Iter=Tbl[min(grpLen*g,nrow(Tbl)),'Iter']
         if(Iter< mxrIter) Iter=Iter-1
         print(paste('Iterations completely done:',Iter))
        # itercIx=cSorted[(1+(Iter-1)*cchLen):min(Iter*cchLen,N)]
        itercIx=cSorted[1:min(Iter*cchLen,N)]#All previous
        cDone[itercIx]=TRUE
        nonTestedCThr=max(mxRAtc[!cDone])*mxrsA
        # iterrIx=rSorted[(1+(Iter-1)*rchLen):min(Iter*rchLen,N)] 
        # IterR= Iter
        if(Iter<= mxrIter) {     
          iterrIx=rSorted[1:min(Iter*rchLen,N)]  #All previous      
          rDone[iterrIx]=TRUE
          nonTestedRThr=max(mxARr[!rDone])*mxrsA
        }
        tge2=proc.time()
        print(sprintf('g=%d Iter=%d, minThr=%f, nonTestedCThr=%f, nonTestedRThr=%f, tcalc=%.1f,tge1=%.1f, tge2=%.1f',g,Iter,minThr,
            nonTestedCThr,nonTestedRThr,(tge0-t0)[3],(tge1-tge0)[3],(tge2-tge1)[3]))
        if(minThr >= nonTestedCThr && minThr >= nonTestedRThr ){
         print("Remaining value are too small, will be ignored")
         break;
        }
   }
   ##
   
    stopCluster(cluster)
     return(list(ikv=ikv,minThr=minThr,Iter=Iter))
   }

