#26/5/2020
 
# Different methods to quantize RESCAL back by getting top scores, 
#    suitable for graphs of less than 5Million triples and 300K entities
# calculate all possible values then choose top values according to parameter

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

#7/1/2020 : wrapper
 RescalReconstructBack<-function(A,R,trpo=NULL,otnsr=NULL,chkZR=FALSE,prd_cnt=NULL,sf=1,
                                 verbose=1,ncores=3,OS_WIN=FALSE,pve=1e-10,grpLen=NULL,ChnkLen=1000,generateLog=FALSE){
    #require(Matrix)
    if(is.null(otnsr) && is.null(trpo)) stop("must provide trpo or otnsr parameters.")
    
        if(is.null(otnsr) && !is.null(trpo)){
            otnsr=getTensor(trpo)
        }
        if(!is.null(otnsr) && is.null(trpo)){
              trpo=tnsr2trp(otnsr)
        }       
        
       if(is.null(prd_cnt)){
           pcnt=table(trpo[,2])
           prd_cnt=rep(0,length(otnsr$P))
           prd_cnt[match(names(pcnt),otnsr$P)]=pcnt
       }
       
       if(chkZR){
            zr=NULL
            for(i in 1:length(R)){
              zr=rbind(zr,data.frame(stringsAsFactors=FALSE,i,prd=otnsr$P[i],prdcnt=prd_cnt[i],mx=max(R[[i]]),mn=min(R[[i]])))
            }
            # prd_cnt1=prd_cnt
            flg=zr[,'mn']==zr[,'mx']
            if(sum(flg)>0){
                if(verbose>0) print("Excluding  following prds:")
                if(verbose>0) print(zr[zr[,'mn']==zr[,'mx'],])
                prd_cnt[zr[zr[,'mn']==zr[,'mx'],'i']]=0 #Exclude fixed prds
            }
       }
       
        t1=proc.time()
        tmp_res=inv_rescal_sf_prd_chnkgrp(R,A,prd_cnt,scale_fact=sf,verbose=verbose,ncores=ncores,ChnkLen=ChnkLen,
                                          grpLen=grpLen,OS_WIN=OS_WIN,pve=pve,generateLog=generateLog)
        t2=proc.time()  
        if(verbose>1) print(t2-t1)

        ijk=tmp_res[[1]]
        trp1=cbind(otnsr$SO[ijk[,1]],otnsr$P[ijk[,2]],otnsr$SO[ijk[,3]])
        flag= paste(unlist(trpo[,1]),unlist(trpo[,2]),unlist(trpo[,3])) %in% paste(unlist(trp1[,1]),unlist(trp1[,2]),unlist(trp1[,3]))
        if(verbose>1) print(table( flag))#Tp
    # trp[flag,]
    return(list(tmp_res,rectrp=trp1,TP=flag))
}

##-----------------------------------------------------------
#3/7/2019
inv_rescal_1prd_thr<-function(R,A,prd,threshold,verbose=1,ncores=3){
# find triples that scores more than
    AT=Matrix::t(A)
    
    n=nrow(A)
    AR=A%*% R[[prd]]

	#require(parallel)
	#require(doParallel)
	cluster <- makeCluster(ncores,outfile="")
	registerDoParallel(cluster)
    
    get1obj<-function(k){
      aiT=AT[,k,drop=FALSE]
      # print(k)
      Vk=as.vector(AR%*%aiT)
      if(sum(Vk>=threshold)>0){
         ijkv=cbind(S=which(Vk>=threshold),O=k,val=Vk[Vk>=threshold])
         return(ijkv)
      }
    }
    ##
    k=NULL#Used to suppress NOTE while checking the package
    tmp=foreach(k =1:n,.packages='Matrix', .combine="rbind") %dopar%{
    # for(k in 1:n){
        get1obj(k)
        # tmp=get1obj(k)
    }
    stopCluster(cluster)
    return(list(ik=tmp[,1:2],val=tmp[,'val']))
}
##-----------------------------------------------------------

inv_rescal_sf_prd_chnkgrp <- function(R,A,cnt_prd,scale_fact=1,verbose=1,ncores=3,ChnkLen=1000,
                                      saveRes=FALSE,OS_WIN=FALSE,pve=1e-10,
      grpLen=NULL, dsname='', rmx=NULL,readjustChnkLen=TRUE,TotalChnkSize=5e8,chTrpCntTol=1000,generateLog=FALSE){
# #NB:to avoid fluctuations in partitions take all required from each part then choose top
#pve: positive e: insure that score is not 0 or less
### NB: much memory required
## grpLen: groups of iterations 
## rmx optional parameter, gives rowMax(abs(A))
   ## TotalChnkSize: max number of elemnts in one chank >=prod(dim(Ach))
   # chTrpCntTol: tolerance in number of triples returned in one chunk, will be eliminated at end of group
     
 getTrpChnk<-function(p_ch){
# a=1;b=0.1;lea=10;leab=100;cnt=30; mid=b+1.0*cnt*(a-b)/(leab-lea)
    quantile_bisect<-function(a,b,cnt,mxIter=30,cntTol=10){#a>b
     lea=sum(t1n_ch_v>=a)
     leb=sum(t1n_ch_v>=b)
      plea=pleb=-1#previous values
     iter=0
     #linear approx:
     # mid=b+1.0*cnt*(a-b)/(leab-lea) 
     while(cnt>lea && cnt<leb && a>b && iter<=mxIter && (!(plea==lea && pleb==leb))){
         #mid=a-1.0*cnt*(a-b)/(leb-lea)
         mid=b+(a-b)/2
         lemid= sum(t1n_ch_v>=mid)
         if((lemid-cnt)< cntTol && lemid>=cnt) {
             return (mid)
             }
         if(lemid<cnt){a=mid;plea=lea;lea=lemid
         }else{b=mid;pleb=leb; leb=lemid}
         iter=iter+1
         print(sprintf('i=%d, ch=%d, a=%f, b=%f, lea=%d, leb=%d, cnt=%d',iter,p_ch,a,b,lea,leb,cnt))
    }
    if(cnt<=lea) return(a)
    return(b)
    } 
    ##
    quantile_topn<-function(cnt){#if there are only few required take it as number of iterations [<30]
         if(cnt>1){
             for(i in 1:(cnt-1)){
               ix=which.max(t1n_ch_v)
               t1n_ch_v[ix]=-1
             }
        }
        return(max(t1n_ch_v))
    }
   
        t1=proc.time()
            # A_ch= A[(1+(ch-1)*ChnkLen):min(ch*ChnkLen,N),,drop=FALSE]
            chIX=grpIX[(1+(p_ch-1)*ChnkLen):min(p_ch*ChnkLen,length(grpIX))]
            A_ch= A[chIX,,drop=FALSE]
            t1n_ch=A_ch%*%RAtMul 
            # thr=quantile(t1n_ch@x,1-qntl)
            t2=proc.time()
            cntCand=sum(t1n_ch>=minThr)#No of candidates
            if(cntCand==0){#cntCand
                 print(sprintf("s=%d,g=%d, p_ch=%d, cntCand==0, qnt will not be calc, skipped",s,g,p_ch))
                 Res=rep(0,3)
            }else{
                t3=proc.time()
            if(cntCand < chTrpCntTol){
               tq0=proc.time()
               thr=minThr
            }else{
                mx=max(t1n_ch)
                cnt=scale_fact*cnt_prd[s]-ifelse(is.null(tmpallg),0,length(tmpallg[tmpallg[,'val']>=mx,1])) #Cnt required
                cnt=min(cnt,cntCand)#count possible            
                if(cnt < 1){
                    print(sprintf("s=%d,g=%d, p_ch=%d,  mx=%f, minThr=%f, chnk skipped",s,g,p_ch,mx,minThr))
                   Res=rep(0,3)
                } else{
                    ##How many less than mx
                    tq0=proc.time()
                    # qntl=scale_fact*cnt_prd[s]/(dim(t1n_ch)[1]*dim(t1n_ch)[2])
                    # thr=quantile(t1n_ch,1-qntl)
                    if(cnt==1){
                      thr=mx 
                    }else{
                      t1n_ch_v=t1n_ch[t1n_ch >= minThr]#make the time for iteration less
                      if( cnt < chTrpCntTol ){#will be eliminated at end of group
                        # thr=quantile_topn(cnt)
                        thr=minThr
                      }else{
                        thr=quantile_bisect(a=mx,b=minThr,cnt=cnt,cntTol=scntTol)
                      }
                    }
                   }
                   }
                    tq1=proc.time()
                    thr=max(thr,minThr)
                    # print(sprintf('-------------%f-------[nv:%d]',thr,sum(t1n_ch>=thr)))
                    print(sprintf('--g=%d,  p_ch=%d------------%f-------',g,p_ch,thr))
                    ik=which(t1n_ch>=thr,arr.ind=TRUE)
                   
                    t4=proc.time()
                    print(sprintf("s=%d, g=%d, p_ch=%d, times t21=%.2f, t32=%.2f, tq02=%.2f, tqntl=%.2f t4q1=%.2f",s,g, p_ch,(t2-t1)[3],(t3-t2)[3],
                                    (tq0-t2)[3],(tq1-tq0)[3],(t4-tq1)[3]))
             
                print(sprintf("prd:%d, #triples=%d, #generated triples=%d",s,cnt_prd[s],nrow(ik)))
                Val=as.vector(t1n_ch[ik])
                
                print(paste("count triples:",nrow(ik)))
                if(nrow(ik)==0){
                   Res=rep(0,3)
                }else{       
                   Res=as.vector(t(cbind(chIX[ik[,1]],which(flgMul)[ik[,2]],Val)) )       
                }
        }
        return(Res)
    }
    #########################
    #require(parallel)
	#require(doParallel)
	# cluster <- makeCluster(ncores,outfile="sf_prd.log")
    if(OS_WIN){
               cluster <- parallel::makeCluster(ncores,outfile=ifelse(generateLog,
                                                                      sprintf("%s_sf_prd_chgrp_nc%d.log",dsname,ncores),''))
        }else{ 
               #not in windows, copy only changed
               cluster <- parallel::makeForkCluster(ncores,outfile=ifelse(generateLog,
                                                                  sprintf("%s_sf_prd_fork_chgrp_nc%d.log",dsname,ncores),''))
    }
	registerDoParallel(cluster)
    N=nrow(A)
    r=ncol(A)
    mxA=max(abs(A))
    cntChnks=ceiling(N/ChnkLen)
    if(is.null(grpLen)) grpLen=5*ncores;
    cntGrp=ceiling(cntChnks/grpLen)
    if(is.null(rmx)) rmx=apply(abs(A),1,max)
    
    rsA=rowSums(abs(A))
    mxrsA=max(rsA)
    rm(rsA)
    print(sprintf("NEnt=%d, ChkLen=%d, cntChnks=%d, cntGrp=%d, grpLen=%d",N, ChnkLen, cntChnks, cntGrp, grpLen))
    act_thr=NULL
    all_res=NULL
    
    for(s in which(cnt_prd>0)){
      if(verbose>0) print(sprintf('-------------=======s=%d==========---------',s))
      ts0 = proc.time()
      # RAt=R[[s]]%*%Matrix::t(A)
      RAt = as.matrix(R[[s]])%*%t(A)
      if(verbose>1) print('calc rowMax..RAt')
      mxRAtr = apply(abs(RAt),2,max)##
      # mxRAtr=colSums(abs(RAt))/r##using mean is more precise than max
      # sum_mxrRAt=sum(apply(abs(RAt),1,max))
      sumRAtc = colSums(abs(RAt))##
      mxSumRAtr = max(sumRAtc)
      # readjusted=FALSE
      minThr = pve     
      # scntTol=max(1,min(100,cnt_prd[s]/10000))
      scntTol = max(10,min(1000,cnt_prd[s]/1000))
      tmpallg = NULL;#data.frame(S=0,P=0,O=0,val=0.0)
      remIX = rmx > (minThr/mxSumRAtr) 
      cntRem = sum(remIX)
      grpIX = which(remIX)[1:min(ChnkLen*grpLen,cntRem)]
      remIX[grpIX] = FALSE
      calcNtrp = 0#Calculated triples so far
      for(g in 1:cntGrp){
          if(verbose>1) print(paste("Group:",g))
          tg0 = proc.time()
          flgMul = mxRAtr > (minThr/mxrsA)
          actNCol=sum(flgMul)
          newChnkLen=floor(TotalChnkSize/actNCol)
          if(verbose>1) print(sprintf('s=%d, g=%d, actNCol=%d, newChnkLen=%d, mxrsA=%f',s,g,actNCol,newChnkLen,mxrsA))
          if(ChnkLen < newChnkLen && readjustChnkLen ){
            ChnkLen=min(cntRem,newChnkLen)
            remIX[grpIX]=TRUE#temporarly reset to TRUE
            grpIX=which(remIX)[1:min(ChnkLen*grpLen,cntRem)]
            remIX[grpIX]=FALSE
            if(verbose>1) print(sprintf('s=%d, g=%d, Chunk size readjusted to:%d',s,g,ChnkLen))
            # readjusted=TRUE#change once
          }
          if(actNCol == 0) break;
          opts <- list(preschedule=FALSE)#not using roundrobin instead when ended pick one from queue
          RAtMul=RAt[,flgMul]
          ch=NULL#Used to suppress NOTE while checking the package
          tmp=foreach(ch =1:ceiling(length(grpIX)/ChnkLen),.packages='Matrix',.options.multicore=opts) %dopar%{#, ,.combine="append"
          # for(ch in 1:ceiling(length(grpIX)/ChnkLen)) {#, ,.combine="append"
               getTrpChnk(p_ch=ch)
          }#dopar
         tg1=proc.time()
         tge0=proc.time()
       tmp=unlist(tmp)
        ix=(1:(length(tmp)/3))-1
         tmp1=data.frame(stringsAsFactors=FALSE,S=as.integer(tmp[3*ix+1]),P=s,O=as.integer(tmp[3*ix+2]),val=as.numeric(tmp[3*ix+3]))
        if(length(tmp1)>0 && any(tmp1[,'S']>0)){
         tmpallg=rbind(tmpallg,tmp1[tmp1[,'S']>0,])
         if(cnt_prd[s]<=nrow(tmpallg)){
            grp_order=order(tmpallg[,'val'],decreasing=TRUE)
            tmpallg=(tmpallg[grp_order,,drop=FALSE])[1:cnt_prd[s],,drop=FALSE]
            calcNtrp=nrow(tmpallg)
            minThr=min(tmpallg[,'val'])                   
            remIX= rmx > (minThr/mxSumRAtr) & remIX# sum_mxrRAt
         }
         # print(minThr)
         }else{ if(verbose>1) print('no new candidate triples found.')}
         ## new grpIX
         
         cntRem=sum(remIX)
         if(cntRem==0 ){
           if(g < cntGrp)   print('remaining values in A are too small, will be skipped ..')
            break;
         }
            grpIX=which(remIX)[1:min(ChnkLen*grpLen,cntRem)]
            remIX[grpIX]=FALSE
            tge1=proc.time()
            if(verbose>1) print(sprintf("=========Group %d, ttime=%.2f, tge=%.2f, cntrows:%d==================",g,(tg1-tg0)[3],(tge1-tge0)[3],nrow(tmp1)))
            if(verbose>1) print(sprintf('remaining no of grp:%d',ceiling(cntRem/(ChnkLen*grpLen))))
        }#grp g
        preccnt=ifelse(is.null(tmpallg),0,nrow(tmpallg))
        if( cnt_prd[s]>preccnt){
            if(verbose>1) print(sprintf("Warning slice:%d Res does not have enough triples (req:%d,calc:%d).",s,cnt_prd[s],preccnt))
        }else{
            if(verbose>1) print(sprintf("**** slice:%d  (req:%d, calc:%d).",s,cnt_prd[s],preccnt))
        }
        all_res=rbind(all_res,tmpallg)
        act_thr=c(act_thr,minThr)
        ts1=proc.time()
        if(verbose>0) print(sprintf("s=%d,   cnt=%d, tmp_cnt:%d, thr:%f ttime=%.2f",s,cnt_prd[s],preccnt,minThr,(ts1-ts0)[3]))
        
        if(saveRes){
          tmpRes=list(ijk=all_res[,1:3],val=all_res[,'val'],act_thr=act_thr)#all_res[,'thr']
          save(file=paste0("inv_rescal_sf_prd_chnk_",s,".RData"),tmpRes)
        }
        # gc()
    }#end s:slice
    stopCluster(cluster)
   return(list(ijk=all_res[,1:3],val=all_res[,'val'],act_thr=act_thr))#ch_thr=all_res[,'thr'],
}
# ---------------------------------------
