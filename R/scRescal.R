# 6/2/2018
# From Yago factorization paper, RESCAL, Nickle code,M Hoffman code

#  Copyright (C) 2018 Abdelmoneim Amer Desouki, 
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

	# Initalize
	# Loop
	# 	updateA
	# 	updateR
	# 	measureProgress
	
	# Settings
	# lambdaA=0.3
	# lambdaR=0.01
	# epsilon=1e-3
	# maxIter=100
    
# ---------- initialize ------------------------------------------------
#Sclable RESCAL
scRescal <- function(X, rnk,ainit = 'nvecs',verbose=2,Ainit=NULL,Rinit=NULL,lambdaA=0,	lambdaR=0, lambdaV=0,
	epsilon=1e-3,maxIter=100, minIter=1, P = list(),orthogonalize=FALSE,func_compute_fval='compute_fit',retAllFact=FALSE,useQR=FALSE,
    ncores=0,OS_WIN=FALSE,savePath='',oneCluster=TRUE,useXc=FALSE,saveARinit=FALSE,saveIter=0,dsname='',maxNrows=50000,
	generateLog=FALSE){
    #retAllFact :flag to return intermediate values of A & R
    #useQR: was suggested by Nickel in Factorizing Yago and implemented by Michail Huffman; found to be converging more slowly.
    #ncores: number of cores used to run in parallel, 0 means no paralllelism
    #useXc: compact the sparse tensor: require more space but much faster.
    #saveIter: iterations in which A&R to be saved, default 0 :none
    #saveARinit: option to save init A and R
    #maxNrows: used as limit to decide if the compact form of the sparse matrix is to be dense matrix (#rows <maxNrows) or to be sparse
    # used in updateA to consider the predicate having rows (in compact form) more than the given number as Big and hence return dense matrix.
    #OS_WIN : True when the operating system is windows, used to allow using Forkin when running in parallel
    
# to be able to continue from a given state Ainit, Rinit can be used
# Parameters
    # ----------
    # X : list
        # List of frontal slices X_k of the tensor X.
        # The shape of each X_k is ('N', 'N').
        # X_k's are expected to be Matrix::sparseMatrix
    # rnk : int
        # Rank of the factorization
    # lmbdaA : float, optional
        # Regularization parameter for A factor matrix. 0 by default
    # lmbdaR : float, optional
        # Regularization parameter for R_k factor matrices. 0 by default
    # lmbdaV : float, optional
        # Regularization parameter for V_l factor matrices. 0 by default
    # init : string, optional
        # Initialization method of the factor matrices. 'nvecs' (default)
        # initializes A based on the eigenvectors of X. 'random' initializes
        # the factor matrices randomly.
    # compute_fit : boolean, optional
        # If true, compute the fit of the factorization compared to X.
        # For large sparse tensors this should be turned of. None by default.
    # maxIter : int, optional
        # Maximium number of iterations of the ALS algorithm. 500 by default.
    # conv : float, optional
        # Stop when residual of factorization is less than conv. 1e-5 by default

    # Returns
    # -------
    # A : ndarray
        # array of shape ('N', 'rank') corresponding to the factor matrix A
    # R : list
        # list of 'M' arrays of shape ('rank', 'rank') corresponding to the
        # factor matrices R_k
    # f : float
        # function value of the factorization
    # itr : int
        # number of iterations until convergence
    # exectimes : ndarray
        # execution times to compute the updates in each iteration
    allA=list()
    allR=list()
	n=nrow(X[[1]])
    if(ncores>0 ){
        #require(parallel)
        #require(doParallel)
    if( oneCluster){
        if(OS_WIN){
               cluster <- parallel::makeCluster(ncores,outfile=ifelse(generateLog,
                                                            sprintf("%s_rescal_par_nc%d.log",dsname,ncores),''))
        }else{ 
                                                  #not in windows, copy only changed
               cluster <- parallel::makeForkCluster(ncores,outfile=ifelse(generateLog,
                                                  sprintf("%s_rescal_par_fork_nc%d.log",dsname,ncores),''))
        }
        doParallel::registerDoParallel(cluster)
        
        if(OS_WIN) clusterExport(cl=cluster, list("sum2sprsMat"), envir=environment())
        }else{cluster=NULL}
    }
	if(!is.null(Ainit)){
        if(ncol(Ainit)!=rnk){
            stop(sprintf("Ainit number of columns (%d) must equal rank (%d)",ncol(Ainit),rnk))
        }
        if(nrow(Ainit)!=n){
            stop(sprintf("Ainit number of rows (%d) must equal n (%d)",nrow(Ainit),n))
        }
		A=Ainit;
        rm(Ainit)
        gc()
	}else{
		if(verbose>0) print('Initializing A')
		if (ainit == 'random'){
			if(verbose>0) print("Initializing A with random values...")
			A <- matrix(runif(n*rnk), ncol=rnk,byrow=TRUE)
			if(verbose>2){print(A)}
			}
		else {
			if (ainit == 'nvecs'){#fails for large matrices
				S = Matrix::sparseMatrix(x=0,i=1,j=1,dims=c(n, n)) #csr_matrix((n, n), dtype=dtype)
                k = length(X)
				
                S=tensor2mat(X,binary=FALSE,symmetrize=TRUE)
				if(verbose>0) print("Calculating  eigen vectors...")
				# tt=eigen(S)#NB: need only rnk vectors but no option to specify it in R eigen
                tt=rARPACK::eigs_sym(S,k=rnk,which="LM",opts = list(retvec = TRUE))
				A = tt$vectors[,1:rnk]
                if(verbose>1 && rnk<=10 && n<=100){
                    print("initial A ................")
                    print(A)
                }
			}
			else{
				stop(sprintf('Unknown init option:%s, should be random or nvecs' ,ainit))
			}
		}
	}
    # --------- use Xc : compact representation of sparse matrices (i.e. slices)
    if(useXc){
        Xc=get_Xic(X,spsThr=0.01,maxMem=1000)
    }
    # ---------------------
    # ------- initialize R and Z ---------------------------------------------
    if(verbose>0) print("initialize R and Z...")
	if(!is.null(Rinit)){
		R=Rinit
	}else{
		t1=proc.time()
		if(useQR){
          if(ncores==0){
                R = updateR_qr(X, A, lambdaR,verbose)
            }else{
                R = updateR_qr_par(X, A, lambdaR,verbose=verbose,ncores=ncores)
            }
        }else{
            if(ncores==0){
                R = updateR(X, A, lambdaR,verbose)
            }else{
                 if(!useXc){
                       R = updateR_par(X, A, lambdaR,verbose,ncores=ncores,clstr=cluster,OS_WIN=OS_WIN)
                 }else{
                      R = updateR_Xc_par(X, Xc, A, lambdaR,verbose,ncores=ncores,clstr=cluster,OS_WIN=OS_WIN)
                  }
                if(saveARinit){
                       save(file=sprintf('%s%s_rescal_init_%s_AR.RData',savePath,dsname,format(Sys.time(),'%Y%m%d%H%M')),R,A)
                }
            }
        }
        
		tUR= proc.time()
		if(verbose>2){print(R)}
		if(verbose>1) print(sprintf('updateR in:%.3f sec',(tUR-t1)[3]))
	}
    Z = updateZ(A, P, lambdaV)

	# precompute norms of X
    normX = lapply(X,function(M){(norm(as.matrix(M@x),'f')^2)})#[sum(M.data ** 2) for M in X]
    sumNorm=0
    for(p in 1:length(X)){
        sumNorm = sumNorm + sum(X[[p]]@x^2)#RDF tensor [sum(M.data ** 2) for M in X]
    }
    if(verbose>1) print(sprintf("sumNorm:%f",sumNorm))
    #  ------ compute factorization ------------------------------------------
    if(verbose>0) print("compute factorization...")
	fit = fitchange = fitold = 0
    exectimes = NULL
	all_err=NULL #maintain a list of errors in every iteration 
    #########Stats to improve performance
    nnzcc=nnzrc=rep(0,length(X))
                for(p in 1:length(X)){
                    if(verbose>2) print(p)
                    Xp=methods::as(X[[p]],"TsparseMatrix")
                    nnzrc[p]=length(unique(Xp@i[Xp@x==1]))
                    nnzcc[p]=length(unique(Xp@j[Xp@x==1]))        
                }
                
    for (itr in 1:maxIter){
		if(verbose>0) print(sprintf("-----------------------------iteration: %d ----------------------------------",itr))
        tic = proc.time()
        fitold = fit
		Aold=A
        gc()
        if(ncores==0){
                A = updateA(X, A, R, P, Z, lambdaA, orthogonalize)
            }else{
                #
                if(!useXc){
                    A = updateA_par(X, A, R, P, Z, lambdaA, orthogonalize,ncores=ncores,verbose=verbose,clstr=cluster,OS_WIN=OS_WIN,maxNrows=maxNrows)
                }else{
                    A = updateA_Xc_par(X, Xc,A, R, P, Z, lambdaA, orthogonalize,ncores=ncores,verbose=verbose,clstr=cluster,OS_WIN=OS_WIN,maxNrows=maxNrows)
                }
        }
        
	    if(verbose>2){print(A)}
		tUA= proc.time()
		if(verbose>1) print(sprintf('updateA in:%.3f sec',(tUA-tic)[3]))
		if(verbose>1) print("Update R...")
        Rold=R
        if(useQR){
            if(ncores==0){
                R = updateR_qr(X, A, lambdaR,verbose)
            }else{
                R = updateR_qr_par(X, A, lambdaR,verbose,ncores=ncores)
            }
        }else{
            if(ncores==0){
                R = updateR(X, A, lambdaR,verbose)
            }else{
                if(!useXc){
                   R = updateR_par(X, A, lambdaR,verbose,ncores=ncores,clstr=cluster,OS_WIN=OS_WIN)
                }else{
                   R = updateR_Xc_par(X, Xc, A, lambdaR,verbose,ncores=ncores,clstr=cluster,OS_WIN=OS_WIN)
                   }
            }
        }
        
		tUR= proc.time()
		if(verbose>1) print(sprintf('updateR in:%.3f sec',(tUR-tUA)[3]))
			if(verbose>2){print(R)}
        Z = updateZ(A, P, lambdaV)

        # compute fit value!!!
		# if(itr%%fent==0){
			if (!is.null(func_compute_fval )){
                if(func_compute_fval=='measureProgress'){
                    fit=measureProgress(Aold,Rold,A,R, ErrThr=epsilon)
                    fitchange=fit
                }else{
                  fit=compute_fit(X, A, R ,normX=sumNorm)#, lambdaA, lambdaR,
                  fitchange = abs(fitold - fit)
                }
                    # print(fit)
			          if(verbose>1) print(sprintf('Calc fit in:%.3f sec',(proc.time()-tUR)[3]))
			}else{
				fit = Inf
			}

        toc = proc.time()
        exectimes=c(exectimes,toc[3] - tic[3])
		all_err=rbind(all_err,cbind(itr,fit=fit,delta=fitchange,Ex_time=exectimes[itr]))
        if(retAllFact){
            allA[[itr]]=A
            allR[[itr]]=R
        } 
        
        if(verbose>0) print(sprintf('[%3d] fit: %0.5f | delta: %7.1e | secs: %.5f',
                     itr, fit, fitchange, exectimes[itr])
        )
        
        if(itr %in% saveIter) save(file=sprintf('%s%s_rescal_i%d_%s_AR.RData',savePath,dsname,itr,format(Sys.time(),'%Y%m%d%H%M')),A,R)
          
        if (itr >= minIter && fitchange < epsilon){
            break;
		}
	}
    
    if(ncores>0 && oneCluster){
        stopCluster(cluster)
    }
    return (list(A=A, R=R, all_err=all_err, nitr=itr , times=as.vector(exectimes),allA=allA,allR=allR))
}

##--------------------------update A--------------------------------------------------
updateA <- function(X, A, R, P, Z, lambdaA, orthogonalize){
    # """Update step for A"""
    n=nrow(A); rnk = ncol(A);
    F = Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(n, rnk))#zeros((n, rnk), dtype=A.dtype)
    E = Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(rnk, rnk))#zeros((rnk, rnk), dtype=A.dtype)

    AtA = Matrix::t(A)%*% A#dot(A.T, A) r x r
	
    for (i in 1:length(X)){
        F = F + X[[i]]%*% (A %*% Matrix::t(R[[i]])) + Matrix::t(X[[i]]) %*% ((A %*% R[[i]]))
        E = E + (R[[i]] %*% (AtA %*% Matrix::t(R[[i]]))) + (Matrix::t(R[[i]]) %*% (AtA %*% R[[i]]))
    }
    # regularization:zeros
    I = lambdaA * diag(rnk)
# save(file=sprintf('%srescal_updateA_EuF_i%d_%s.RData',get("savePath", parent.frame()),get("itr", parent.frame()),format(Sys.time(),'%Y%m%d%H%M')),F,E)
    
    # attributes
	if(length(Z)>0){
		for( i in 1:length(Z)){
			F = F + (P[[i]] %*% (Matrix::t(Z[[i]])))
			E = E + (Z[[i]] %*% (Matrix::t(Z[[i]])))#dot(Z[i], Z[i].T)
		}
	}
    # finally compute update for A  O(r^3), solve eqns: E A = F, for A
    A = Matrix::t(Matrix::solve(I + Matrix::t(E), Matrix::t(F)))
    if(orthogonalize) {
        stop("orth(A) not implemented.")
	 }else{
        return(A)
     }
}

# ---------------------------------------------------------------------------
# Tensor_norm
# ---------------------------------------------------------------------------
Tensor_error <- function(X, A, R){
	err= 0 
	for (i in 1:length(X)){
        t0=proc.time()
        ARAt = A %*% (R[[i]]%*% Matrix::t(A))
		tmp=X[[i]] - ARAt
        err = err + (Matrix::norm(tmp,'f') ^ 2)  
        t1=proc.time()
		print(sprintf("i=%d, Err=%f, t=%.3f",i,err,(t1-t0)[3]))
        rm(tmp)#memory issues
		gc()
	}
	return (err)
}

# ------------------ Update R ------------------------------------------------
 updateR <- function(X, A, lambdaR,verbose=1){
    rnk = ncol(A)

	if(verbose>1) print('calculating SVD of A...')
	t1=proc.time()
	stmp=svd(x=A)#default,nu=min(dim(A)),nv=min(dim(A))
	U=stmp$u
	S=matrix(stmp$d,nrow=1)
	Vt=Matrix::t(stmp$v)
	t2=proc.time()
    if(verbose>1){
        print(sprintf('time calc. SVD of A:%.3f',(t2-t1)[3]))
        print('Calculating Shat...')
    }
    Shat = kronecker(S, S)#kronecker product: element-wise product d x d
    Shat = matrix(Shat / (Shat^2 + lambdaR),ncol=rnk, nrow=rnk)

	if(verbose>1) print('Updating R...');
    R = list()
    for(i in 1:length(X)){
        # Rn = Shat * dot(U.T, X[i].dot(U))
        Rn = Shat * (Matrix::t(U) %*% (X[[i]] %*% U))
        Rn = Matrix::t(Vt) %*% (Rn %*% Vt)
        R[[i]]=Rn
	}
    return (R)
}
##--------------------------update A--------------------------------------------------
updateA_par <- function(X, A, R, P, Z, lambdaA, orthogonalize,verbose=1,ncores,clstr=NULL,OS_WIN=FALSE,maxNrows=200000){#grpLen=40,
    # """Update step for A"""
    n=nrow(A); rnk = ncol(A);
    # F = Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(n, rnk))#zeros((n, rnk), dtype=A.dtype)
    # E = Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(rnk, rnk))#zeros((rnk, rnk), dtype=A.dtype)
    calc_F1<-function(Xi,Ri,p_i){#15 sec , 200+ secs
       t0=proc.time()
       tmp = Xi %*% A 
       t1=proc.time()
          F1i = tmp %*% Ri
          
       t2=proc.time()
       rm(tmp)
       
       if(verbose>1) print(sprintf("Fun calc_F1, i=%d, p=%d, ntrp=%d , t_tmp: %.f secs, t_F1i:%.2f secs, sz(F1)=%.fMB...",i,p_i,sum(Xi@x==1),(t1-t0)[3],(t2-t1)[3],object.size(F1i)/(2.0^20)))
       return(F1i)
    }
    
	F1ViaSum<-function(Xi,Ri,p_i){
       t0=proc.time()
        Xi=methods::as(Xi,"TsparseMatrix")

        ii=Xi@i[Xi@x==1]+1
        iis=sort(unique(ii))
        jj=Xi@j[Xi@x==1]+1
        print(paste("i=",p_i,length(ii),length(iis)))
        # FL  = Xi[iis,,drop=FALSE] %*% A
        tmpl=lapply(iis,function(x){colSums(A[jj[ii==x],,drop=FALSE])})
        length(tmpl) 
        tmpL1=unlist(tmpl)
        #is(tmpl[[1]])
        FL=matrix(tmpL1,ncol=ncol(A),byrow=TRUE)
        t1=proc.time()
        F1tmp = FL %*% Ri
        # if(p_i==8)print('Ok 1')
        spsRatio=(sum(F1tmp!=0.0)/prod(dim(A)))
        if(verbose>1) print(sprintf(" i=%d , p=%d, sprsRatio=%f",i,p_i,spsRatio ))
        if(spsRatio < 0.1){#sparsity limit
            ix=which(F1tmp!=0,arr.ind=TRUE)
          if(length(ix)==0){
              F1i=sparseMatrix(x=0,i=1,j=1,dims=c(nrow(Xi),ncol(Ri)))
              if(verbose>1) print(sprintf("Fun F1ViaSum, i=%d, p=%d, ntrp=%d , t_tmp: %.f secs, Fli=0 nonzerovalues...",i,p_i,sum(Xi@x==1),(t1-t0)[3]))
          }else{
            # if()
            # val=as.vector(t(F1tmp))#byrow
            val=F1tmp[ix]
            F1i=sparseMatrix(x=val,i=iis[ix[,1]],j=ix[,2],dims=c(nrow(Xi),ncol(Ri)))
            rm(ix)
            rm(val)
            }
        }else{
            F1i=matrix(0,nrow=nrow(Xi),ncol=ncol(Ri))
            F1i[iis,]=F1tmp#no use of ,drop=FALSE
        }
        t2=proc.time()
       rm(F1tmp)
       rm(tmpL1)
       rm(FL)
       # print(gc())
       if(verbose>1) print(sprintf("Fun F1ViaSum, i=%d, p=%d, ntrp=%d , t_tmp: %.f secs, t_F1i:%.2f secs, sz(F1i)=%.fMB...",i,p_i,sum(Xi@x==1),(t1-t0)[3],(t2-t1)[3],object.size(F1i)/(2.0^20)))
      
        return(F1i)
    }
    # # regularization:zeros
    if(verbose>1) print("Calc F ... ")
    # n=nrow(X[[1]])
    grpL=which(get("nnzrc", parent.frame())<=maxNrows)#tested to be less than 10 sec
    grpB=which(get("nnzrc", parent.frame())>maxNrows)
   
    if(is.null(clstr)){
        if(OS_WIN){
               cluster <- parallel::makeCluster(ncores,outfile=sprintf("rescal_par_nc%d.log",ncores))
        }else{ 
               cluster <- parallel::makeForkCluster(ncores,outfile=sprintf("rescal_updateA_par_fork_nc%d.log",ncores))#not in windows, copy only changed
        }
        doParallel::registerDoParallel(cluster)
       }
        FgrpB=NULL
    if(length(grpB)>0){
    grpLen=3*ncores
    print(sprintf("Processing predicates having high number of rows t(Ri)[%d]..: grpLen %d",length(grpB),grpLen))
    for(g in 1:ceiling(length(grpB)/grpLen)){
        tf0=proc.time()
        
        print(paste("Group: ",g))
        FgrpB1g=foreach(i =((g-1)*grpLen+1):min(length(grpB),(g*grpLen)),.packages='Matrix', .combine="+") %dopar%{
            p=grpB[i]
            if(verbose>1) print(sprintf("Calc F: i=%d, p=%d, mem:%s, %s",i,p,pryr::mem_used(),format(Sys.time(),"%H:%M:%S")))        
            Xpi=X[[p]]
            Fi =calc_F1(Xpi,t(R[[p]]),p)
          }
        tf1=proc.time()
        if(is.null(FgrpB)){
            FgrpB=FgrpB1g
        }else{
            FgrpB=FgrpB+FgrpB1g
        }

        rm(FgrpB1g)
        if(verbose>1) 
           print(sprintf("Calc FgrpB grp=%d in %.f secs, sz(F1)=%.fMB, mem_used=%.2fGB...",g,(tf1-tf0)[3],object.size(FgrpB)/(2.0^20),pryr::mem_used()/(2.0^30)))
    }#g  t(Ri)
    
     print(sprintf("Processing predicates having high number of rows F2: t(Xi)[%d]..: grpLen %d",length(grpB),grpLen))
    for(g in 1:ceiling(length(grpB)/grpLen)){
        tf0=proc.time()
        
        print(paste("Group: ",g))
        # F1=foreach(i =1:length(X),.packages='Matrix', .combine="+") %dopar%{
        FgrpB1g=foreach(i =((g-1)*grpLen+1):min(length(grpB),(g*grpLen)),.packages='Matrix', .combine="+") %dopar%{
            p=grpB[i]
            if(verbose>1) print(sprintf("Calc F t(Xi): i=%d, p=%d, mem:%s, %s",i,p,pryr::mem_used(),format(Sys.time(),"%H:%M:%S")))        
            Xpi=X[[p]]
            Fi =calc_F1(Matrix::t(Xpi),R[[p]],p)
          }
        tf1=proc.time()
        if(is.null(FgrpB)){
            FgrpB=FgrpB1g
        }else{
            FgrpB=FgrpB+FgrpB1g
        }
        rm(FgrpB1g)
        if(verbose>1) 
           print(sprintf("Calc FgrpB grp=%d in %.f secs, sz(F1)=%.fMB, mem_used=%.2fGB...",g,(tf1-tf0)[3],object.size(FgrpB)/(2.0^20),pryr::mem_used()/(2.0^30)))
    }#g  t(Xi)
   
    }
   
    #grpL: low no of rows, no need to grouping
    tl0=proc.time()
    FgrpL=NULL
    grpLen=4*ncores #
    if(length(grpL)>0){
      if(verbose>1)print(sprintf("Processing predicates having LOW number of rows (<%d)[cnt:%d]..: grpLen %d",maxNrows,length(grpL),grpLen))
        for(g in 1:ceiling(length(grpL)/grpLen)){
          if(verbose>1) print(paste("grpL Group: ",g))
        # FgrpL=foreach(i =1:length(grpL),.packages='Matrix', .combine="+") %dopar%{
        tf0=proc.time()
        FgrpL1g=foreach(i =((g-1)*grpLen+1):min(length(grpL),(g*grpLen)),.packages='Matrix', .combine="+") %dopar%{
                p=grpL[i]
                if(verbose>1) print(sprintf("Calc F: i=%d, p=%d, mem:%s, %s",i,p,pryr::mem_used(),format(Sys.time(),"%H:%M:%S")))        
                Xpi=X[[p]]
                Fi =F1ViaSum(Xpi,t(R[[p]]),p)+F1ViaSum(Matrix::t(Xpi),R[[p]],p)
              }
              if(is.null(FgrpL)){
                FgrpL=FgrpL1g
            }else{
                FgrpL=FgrpL+FgrpL1g
            }
            # rm(Fi)
            rm(FgrpL1g)
            # print(gc())
            tf1=proc.time()
            if(verbose>1) 
               # print(sprintf("Calc F1 grp=%d in %.f secs, sz(F1)=%.fMB, mem_used=%.2fGB...",g,(tf1-tf0)[3],object.size(F1)/(2.0^20),pryr::mem_used()/(2.0^30)))
               print(sprintf("Calc FgrpL grp=%d in %.f secs, sz(FgrpL)=%.fMB, mem_used=%.2fGB...",g,(tf1-tf0)[3],object.size(FgrpL)/(2.0^20),pryr::mem_used()/(2.0^30)))
          }#
            tl1=proc.time()
        
            if(verbose>1) print(sprintf('UpdateA, time in calc grpL[%d]: low no of rows %.2f, size of FgrpL:%.fMB',length(grpL),(tl1-tl0)[3],object.size(grpL)/(2^20)))
        }
    # save(file=sprintf('%srescal_updateA_F1_i%d_%s.RData',get("savePath", parent.frame()),get("itr", parent.frame()),format(Sys.time(),'%Y%m%d%H%M')),F1)
    ###################################################
    if(!is.null(FgrpB)){
       if(!is.null(FgrpL)){
           mF=FgrpB + FgrpL
       }else{  mF = FgrpB }
    }else{  mF = FgrpL}
    
    rm(FgrpB)
    rm(FgrpL)
     
    if(verbose>1) print("Calc AtA...")
    t0=proc.time()
    AtA = Matrix::t(A)%*% A#dot(A.T, A) r x r
    t1=proc.time()
    if(verbose>1) print(sprintf("Calc AtA in %.f secs...",(t1-t0)[3]))
     te0=proc.time()
       
    if(verbose>1) print("Calc E ... ")      
    E=foreach(i =1:length(X),.packages='Matrix', .combine="+") %dopar%{
        if(verbose>1) print(sprintf("Calc E: i=%d, mem:%s, %s",i,pryr::mem_used(),format(Sys.time(),"%H:%M:%S")))   
        Ei= (R[[i]] %*% (AtA %*% Matrix::t(R[[i]]))) + (Matrix::t(R[[i]]) %*% (AtA %*% R[[i]]))
      }
   te1=proc.time()
   if(verbose>1) 
       print(sprintf("updateA: Calc E in %.f secs, sz(E)=%.fMB, mem_used=%.2fGB...",(te1-te0)[3],object.size(E)/(2.0^20),pryr::mem_used()/(2.0^30)))
   
    I = lambdaA * diag(rnk)##lambdA's are zeros

    # attributes
	if(length(Z)>0){
		for( i in 1:length(Z)){
			mF = mF + (P[[i]] %*% (Matrix::t(Z[[i]])))
			E = E + (Z[[i]] %*% (Matrix::t(Z[[i]])))#dot(Z[i], Z[i].T)
		}
	}
    # finally compute update for A  O(r^3), solve eqns: E A = F, for A
    if(verbose>1) print("UpdateA_par: Solving for A ...");
    ts0=proc.time()
    A = Matrix::t(Matrix::solve(I + Matrix::t(E), Matrix::t(mF)))
    ts1=proc.time()
    if(verbose>1) print(sprintf("UpdateA: time in solving for A:%.2f secs",(ts1-ts0)[3]))
    if(is.null(clstr)){
        stopCluster(cluster)
    }
    if(orthogonalize) {
        stop("orth(A) not implemented.")
	 }else{
        return(as.matrix(A))
     }
}

##--------------------------Calc Xc--------------------------------------------------
get_Xic<-function(X,spsThr=0.01,maxMem=1000){#500 MB
# Calculate compact form of slices, improves the speed of single iterations by factor 4 for large tensors.
# spsThr: threshold of sparsity, slices with sparsity>spsThr will be dense
# maxMem: maximum memory allowed for one slice as compact matrix,
    allXic=list()
    for(p in 1:length(X)){
        Xi=X[[p]]
        ii=Xi@i[Xi@x==1]+1
        iis=sort(unique(ii))
        jj=Xi@j[Xi@x==1]+1
        jjs=sort(unique(jj))
        
        Xicsps=sum(Xi@x==1)/(length(iis)*1.0*length(jjs))
        mem=length(iis)*10.0*length(jjs)/(2^20)
        if(Xicsps >= spsThr & mem<=maxMem){
           Xic=matrix(0,nrow=length(iis),ncol=length(jjs))
           Xic[cbind(match(ii,iis),match(jj,jjs))]=1
        }else{
            Xic=Matrix::sparseMatrix(i=match(ii,iis),j=match(jj,jjs),x=1,dims=c(length(iis),length(jjs)))
        }
        allXic[[p]]=list(Xic=Xic,iis=iis,jjs=jjs)
    }
    
   print(sprintf("Size of all Xic=%.2fMB",object.size(allXic)/(2^20)))
   return(allXic)
}
##---------------------------------------------------------------------------------
sum2sprsMat<-function(m1,m2){
        m1=methods::as(m1,'TsparseMatrix')
        m2=methods::as(m2,'TsparseMatrix')
        alli=append(m1@i[m1@x!=0],m2@i[m2@x!=0])+1
        allj=append(m1@j[m1@x!=0],m2@j[m2@x!=0])+1
        allv=append(m1@x[m1@x!=0],m2@x[m2@x!=0])        

        m12=Matrix::sparseMatrix(i=alli,j=allj,x=allv,dims=c(nrow(m1), ncol(m1))) 
        m12=methods::as(m12,'TsparseMatrix')
      return(m12)
}

sumSprsMat_par<-function(){#list of matrices in matLst variable
  # `%dopar%` <- foreach::`%dopar%`
    mL=get("matLst", parent.frame())
    if(length(mL)==1) return(mL[[1]]);
    while(1==1){
      # print(length(mL))
      i=NULL#Used to suppress NOTE while checking the package
       mL1=foreach::foreach(i=1:floor(length(mL)/2),.packages='Matrix') %dopar%{
          print(i)
          # m1=mL[[2*i-1]]   
          # m2=mL[[2*i]]
          # m12=m1+m2
          sum2sprsMat(mL[[2*i-1]],mL[[2*i]])
        }
     
       if(length(mL)%%2==1) {
            mL= append(mL1,mL[[length(mL)]])
       }else{
         mL=mL1
       }
       
        if(length(mL)<2) break;
    }
 return(mL[[1]])
}
##---------------------------------------------------------------------------------
updateA_Xc_par <- function(X, Xc, A, R, P, Z, lambdaA, orthogonalize,verbose=1,ncores,clstr=NULL,OS_WIN=FALSE,maxNrows=50000){#grpLen=40,
    # """Update step for A"""
    n=nrow(A); rnk = ncol(A);
    
	F1Xc<-function(Xic,iis,jjs,Ri,p_i,retDenseMat=FALSE){
        t0=proc.time()
        Ac=A[jjs,,drop=FALSE]
        FL=Xic%*%Ac
        t1=proc.time()
        #print(t1-t0)#28-38
        FL=as.matrix(FL)
        F1tmp = FL %*% Ri
        t2=proc.time()
    spsRatio=(sum(F1tmp!=0.0)/prod(dim(A)))
        if(spsRatio < 0.01 && !retDenseMat){#sparsity limit
            ix=which(F1tmp!=0,arr.ind=TRUE)
          if(length(ix)==0){
              F1i=sparseMatrix(x=0,i=1,j=1,dims=c(nE,ncol(Ri)))
              print(sprintf("Fun F1Xc, i=%d, p=%d, ntrp=%d , t_tmp: %.f secs, Fli=0 nonzerovalues...",p,p,sum(Xic),(t2-t1)[3]))
          }else{
            # if()
            # val=as.vector(t(F1tmp))#byrow
            val=F1tmp[ix]
            F1i=sparseMatrix(x=val,i=iis[ix[,1]],j=ix[,2],dims=c(nE,ncol(Ri)))
            rm(ix)
            rm(val)
            }
        }else{
            F1i=matrix(0,nrow=nE,ncol=ncol(Ri))
            F1tmp=as.matrix(F1tmp)
            F1i[iis,]=F1tmp#no use of ,drop=FALSE
        }
      t3=proc.time() 
       print(sprintf("p=%d, nrows=%d, spsRatio=%.2f timeFL=%.2f, timeF1i=%.2f",p_i,length(iis),spsRatio,(t3-t1)[3],(t1-t0)[3]))
        return(F1i)
    }
    # # regularization:zeros
    if(verbose>1) print("Calc F ... ")
    # n=nrow(X[[1]])
    grpL=which(get("nnzrc", parent.frame())<=maxNrows)#tested to be less than 10 sec
    grpB=which(get("nnzrc", parent.frame())>maxNrows)
    nE=nrow(X[[1]])
    tmpo=get("nnzrc", parent.frame())*1.0*get("nnzcc", parent.frame())
    tmpoL=order(tmpo[grpL])#Let the bigger matrices come last
    tmpoB=order(tmpo[grpB])
    if(is.null(clstr)){
        if(OS_WIN){
               cluster <- parallel::makeCluster(ncores,outfile=sprintf("rescal_par_nc%d.log",ncores))
        }else{ 
               cluster <- parallel::makeForkCluster(ncores,outfile=sprintf("rescal_updateA_par_fork_nc%d.log",ncores))#not in windows, copy only changed
        }
        doParallel::registerDoParallel(cluster)
        
        if(OS_WIN) clusterExport(cl=cluster, list("sum2sprsMat"), envir=environment())
       }
        FgrpB=NULL
    if(length(grpB)>0){
    grpLen=3*ncores
    print(sprintf("Processing predicates having high number of rows t(Ri)[%d]..: grpLen %d",length(grpB),grpLen))
    for(g in 1:ceiling(length(grpB)/grpLen)){
        tf0=proc.time()
        
        print(paste("Group: ",g))
        FgrpB1g=foreach(i =((g-1)*grpLen+1):min(length(grpB),(g*grpLen)),.packages='Matrix', .combine="+") %dopar%{
            p=grpB[tmpoB][i]
            if(verbose>1) print(sprintf("Calc F: i=%d, p=%d, mem:%s, %s",i,p,pryr::mem_used(),format(Sys.time(),"%H:%M:%S")))        
            
            Fi =F1Xc(Xc[[p]][['Xic']],Xc[[p]][['iis']],Xc[[p]][['jjs']],t(R[[p]]),p,retDenseMat=TRUE)
          }
        tf1=proc.time()
        if(is.null(FgrpB)){
            FgrpB=FgrpB1g
        }else{
            FgrpB=FgrpB+FgrpB1g
        }
        rm(FgrpB1g)
        if(verbose>1) 
           print(sprintf("Calc FgrpB grp=%d in %.f secs, sz(F1)=%.fMB, mem_used=%.2fGB...",g,(tf1-tf0)[3],object.size(FgrpB)/(2.0^20),pryr::mem_used()/(2.0^30)))
    }#g  t(Ri)
    
     print(sprintf("Processing predicates having high number of rows F2: t(Xi)[%d]..: grpLen %d",length(grpB),grpLen))
    for(g in 1:ceiling(length(grpB)/grpLen)){
        tf0=proc.time()
        
        print(paste("Group: ",g))
        FgrpB1g=foreach(i =((g-1)*grpLen+1):min(length(grpB),(g*grpLen)),.packages='Matrix', .combine="+") %dopar%{
            p=grpB[tmpoB][i]
            if(verbose>1) print(sprintf("Calc F t(Xi): i=%d, p=%d, mem:%s, %s",i,p,pryr::mem_used(),format(Sys.time(),"%H:%M:%S")))        
            Fi =F1Xc(Matrix::t(Xc[[p]][['Xic']]),Xc[[p]][['jjs']],Xc[[p]][['iis']],R[[p]],p,retDenseMat=TRUE)
          }
        tf1=proc.time()
        if(is.null(FgrpB)){
            FgrpB=FgrpB1g
        }else{
            FgrpB=FgrpB+FgrpB1g
        }
        rm(FgrpB1g)
        if(verbose>1) 
           print(sprintf("Calc FgrpB grp=%d in %.f secs, sz(F1)=%.fMB, mem_used=%.2fGB...",g,(tf1-tf0)[3],object.size(FgrpB)/(2.0^20),pryr::mem_used()/(2.0^30)))
    }#g  t(Xi)
   
    }
   
    #grpL: low no of rows, no need to grouping
    tl0=proc.time()
    FgrpL=NULL
    grpLen=3*ncores #
    if(length(grpL)>0){
        print(sprintf("Processing predicates having LOW number of rows (<%d)[cnt:%d]..: grpLen %d",maxNrows,length(grpL),grpLen))
        for(g in 1:ceiling(length(grpL)/grpLen)){
          if(verbose>1) print(paste("grpL Group: ",g))
        tf0=proc.time()
        matLst=foreach(i =((g-1)*grpLen+1):min(length(grpL),(g*grpLen)),.packages='Matrix') %dopar%{
                p=grpL[tmpoL][i]
                if(verbose>1) print(sprintf("Calc F: i=%d, p=%d, mem:%s, %s",i,p,pryr::mem_used(),format(Sys.time(),"%H:%M:%S")))        
                Xpi=X[[p]]
                Fi=sum2sprsMat(F1Xc(Xc[[p]][['Xic']],Xc[[p]][['iis']],Xc[[p]][['jjs']],t(R[[p]]),p),
                               F1Xc(Matrix::t(Xc[[p]][['Xic']]),Xc[[p]][['jjs']],Xc[[p]][['iis']],R[[p]],p) )
        }
        
        tf1=proc.time()
        if(is.null(FgrpL)){
            FgrpL=sumSprsMat_par()
        }else{
            FgrpL=FgrpL+sumSprsMat_par()
        }
        tf2=proc.time()
            if(verbose>1) 
               print(sprintf("Calc FgrpL grp=%d in %.2f secs, sumSprs time=%.2f, sz(FgrpL)=%.fMB, mem_used=%.2fGB...",g,(tf2-tf0)[3],(tf2-tf1)[3],object.size(FgrpL)/(2.0^20),pryr::mem_used()/(2.0^30)))
          }#
            tl1=proc.time()
        
            if(verbose>1) print(sprintf('UpdateA, time in calc grpL[%d]: low no of rows %.2f, size of FgrpL:%.fMB',length(grpL),(tl1-tl0)[3],object.size(FgrpL)/(2^20)))
        }
    # save(file=sprintf('%srescal_updateA_F1_i%d_%s.RData',get("savePath", parent.frame()),get("itr", parent.frame()),format(Sys.time(),'%Y%m%d%H%M')),F1)
    ###################################################
    if(!is.null(FgrpB)){
       if(!is.null(FgrpL)){
           mF=FgrpB + FgrpL
       }else{  mF = FgrpB }
    }else{  mF = FgrpL}
    
    
    # save(file=sprintf('%srescal_updateA_par_mF_i%d_%s.RData',get("savePath", parent.frame()),get("itr", parent.frame()),format(Sys.time(),'%Y%m%d%H%M')),mF)
    rm(FgrpB)
    rm(FgrpL)
     
    if(verbose>1) print("Calc AtA...")
    t0=proc.time()
    AtA = Matrix::t(A)%*% A#dot(A.T, A) r x r
    t1=proc.time()
    if(verbose>1) print(sprintf("Calc AtA in %.f secs...",(t1-t0)[3]))
     te0=proc.time()
        
    if(verbose>1) print("Calc E ... ")      
    E=foreach(i =1:length(X),.packages='Matrix', .combine="+") %dopar%{
        if(verbose>1) print(sprintf("Calc E: i=%d, mem:%s, %s",i,pryr::mem_used(),format(Sys.time(),"%H:%M:%S")))   
        Ei= (R[[i]] %*% (AtA %*% Matrix::t(R[[i]]))) + (Matrix::t(R[[i]]) %*% (AtA %*% R[[i]]))
      }
   te1=proc.time()
   if(verbose>1) 
       print(sprintf("updateA: Calc E in %.f secs, sz(E)=%.fMB, mem_used=%.2fGB...",(te1-te0)[3],object.size(E)/(2.0^20),pryr::mem_used()/(2.0^30)))
           
    ####
    I = lambdaA * diag(rnk)##lambdA's are zeros

    # attributes
	if(length(Z)>0){
		for( i in 1:length(Z)){
			mF = mF + (P[[i]] %*% (Matrix::t(Z[[i]])))
			E = E + (Z[[i]] %*% (Matrix::t(Z[[i]])))#dot(Z[i], Z[i].T)
		}
	}
    # finally compute update for A  O(r^3), solve eqns: E A = F, for A
    if(verbose>1) print("UpdateA_par: Solving for A ...");
    ts0=proc.time()
    A = Matrix::t(Matrix::solve(I + Matrix::t(E), Matrix::t(mF)))
    ts1=proc.time()
    if(verbose>1) print(sprintf("UpdateA: time in solving for A:%.2f secs",(ts1-ts0)[3]))
    if(is.null(clstr)){
        stopCluster(cluster)
    }
    if(orthogonalize) {
        stop("orth(A) not implemented.")
        # return(orth(A));
	 }else{
        return(as.matrix(A))
     }
}
# ------------------ Update R ------------------------------------------------
 updateR_par <- function(X, A, lambdaR,verbose=1,ncores,clstr=NULL,OS_WIN=FALSE){
    rnk = ncol(A)#A.shape[1]

	if(verbose>1) print('calculating SVD of A...')
	t1=proc.time()
	stmp=svd(x=A)#default,nu=min(dim(A)),nv=min(dim(A))
	U=stmp$u
	S=matrix(stmp$d,nrow=1)
	Vt=Matrix::t(stmp$v)
	t2=proc.time()
    if(verbose>1){
        print(sprintf('time calc. SVD of A:%.3f',(t2-t1)[3]))
        print('Calculating Shat...')
    }
    Shat = kronecker(S, S)#kronecker product: element-wise product
    Shat = matrix(Shat / (Shat^2 + lambdaR),ncol=rnk, nrow=rnk)

	if(verbose>1) print(sprintf('Updating R par ncores:%d ...sz(U)=%.2f,sz(Vt)=%.2f, sz(Shat)=%.2f',ncores,
                                          object.size(U)/(1024*1024),object.size(Vt)/(1024*1024) ,object.size(Shat)/(1024*1024) ));
    Ut=Matrix::t(U)
    if(is.null(clstr)){
        if(OS_WIN){
               cluster <- parallel::makeCluster(ncores,outfile=sprintf("rescal_par_nc%d.log",ncores))
        }else{ 
               cluster <- parallel::makeForkCluster(ncores,outfile=sprintf("rescal_updateR_par_fork_nc%d.log",ncores))#not in windows, copy only changed
        }
        doParallel::registerDoParallel(cluster)
    }
    #idea: better multiply: mm=cbind(r=X[[i]]@i,v=U[X[[i]]@j+1,]), by(mm[,'v'],mm[,'r'],sum),
    i=NULL#Used to suppress NOTE while checking the package
    tmpR=foreach(i =1:length(X),.packages='Matrix') %dopar%{
        print(sprintf('-------- %d ------',i))
        Rn = Shat * (Ut %*% (X[[i]] %*% U))
        print(sprintf("sz(Rn)=%.2f",object.size(Rn)/(1024*1024)))
        Rn = as.matrix(Matrix::t(Vt) %*% (Rn %*% Vt))
        # R[[i]]=Rn 
	}
    if(is.null(clstr)){
        stopCluster(cluster)
    }
    return (tmpR)
}
# ------------------ Update R Xc------------------------------------------------
 updateR_Xc_par <- function(X, Xc, A, lambdaR,verbose=1,ncores,clstr=NULL,OS_WIN=FALSE){
    rnk = ncol(A)#A.shape[1]
    calcUtXU<-function(p_i){
        t0=proc.time()
        tmp=Xc[[p_i]][['Xic']] %*% U[Xc[[p_i]][['jjs']],,drop=FALSE]
        t1=proc.time()
       tmp2=(Ut[,Xc[[p_i]][['iis']],drop=FALSE] %*% tmp)
       t2=proc.time()
       rm(tmp)
       print(sprintf('p=%d, nrows=%d, ttmp=%.2f , ttmp2=%.2f',p_i,length(Xc[[p_i]][['iis']]),(t1-t0)[3],(t2-t1)[3]))
       return(tmp2)
    }
	if(verbose>1) print('calculating SVD of A...')
	t1=proc.time()
	stmp=svd(x=A)#default,nu=min(dim(A)),nv=min(dim(A))
	U=stmp$u
	S=matrix(stmp$d,nrow=1)
	Vt=Matrix::t(stmp$v)
	t2=proc.time()
    if(verbose>1){
        print(sprintf('time calc. SVD of A:%.3f',(t2-t1)[3]))
        print('Calculating Shat...')
    }
    Shat = kronecker(S, S)#kronecker product: element-wise product
    Shat = matrix(Shat / (Shat^2 + lambdaR),ncol=rnk, nrow=rnk)

	if(verbose>1) print(sprintf('Updating R par ncores:%d ...sz(U)=%.2f,sz(Vt)=%.2f, sz(Shat)=%.2f',ncores,
                                          object.size(U)/(1024*1024),object.size(Vt)/(1024*1024) ,object.size(Shat)/(1024*1024) ));
    # R = list()
    # for(i in 1:length(X)){
    Ut=Matrix::t(U)
    if(is.null(clstr)){
        if(OS_WIN){
               cluster <- parallel::makeCluster(ncores,outfile=sprintf("rescal_par_nc%d.log",ncores))
        }else{ 
               cluster <- parallel::makeForkCluster(ncores,outfile=sprintf("rescal_updateR_par_fork_nc%d.log",ncores))#not in windows, copy only changed
        }
        doParallel::registerDoParallel(cluster)
    }
    #idea: better multiply: mm=cbind(r=X[[i]]@i,v=U[X[[i]]@j+1,]), by(mm[,'v'],mm[,'r'],sum),
    i=NULL#Used to suppress NOTE while checking the package
    tmpR=foreach(i =1:length(X),.packages='Matrix') %dopar%{
        print(sprintf('-------- %d ------',i))
        # Rn = Shat * (Ut %*% (X[[i]] %*% U))
        Rn = Shat * calcUtXU(i)
        print(sprintf("sz(Rn)=%.2f",object.size(Rn)/(1024*1024)))
        Rn = as.matrix(Matrix::t(Vt) %*% (Rn %*% Vt))
        # R[[i]]=Rn 
	}
    if(is.null(clstr)){
        stopCluster(cluster)
    }
    return (tmpR)
}
# ------------------ Update R Unreg using QR------------------------------------------------
 updateR_qr <- function(X, A, R,verbose=0){
	# computes the updates of the Core tensor
	unregularizedRUpdate<-function(AL,BL){
		#Here we solve ASA.T = X, using that A is upper triangular
		# x <- backsolve (R, b) solves R x = b where R is upper triangular
		M1 = backsolve(AL, BL)#, overwrite_b=True, check_finite=False
		return (Matrix::t(backsolve(AL, Matrix::t(M1))))
	}
	# Q, A_hat = np.linalg.qr(A, mode='reduced')
	qrR<-qr(A, LAPACK = TRUE)
	Q <- qr.Q(qrR)
	A_hat <- qr.R(qrR)
	R=list()
	for (i in 1:length(X)){ 
		X_hat = Matrix::t(Q) %*% (X[[i]]%*%Q)
		# R[[i]] = unregularizedRUpdate(A_hat, X_hat)
        M1 = backsolve(A_hat, X_hat)
        R[[i]] = Matrix::t(backsolve(A_hat, Matrix::t(M1)))
	}
	return (R)
}
# ------------------ Update R Unreg using QR  parallel------------------------------------------------
 updateR_qr_par <- function(X, A, R,verbose=0,ncores){
	# computes the updates of the Core tensor
	# unregularizedRUpdate
    
	# Q, A_hat = np.linalg.qr(A, mode='reduced')
    t1=proc.time()
    if(verbose>1) print(sprintf('Updating R QR par ncores:%d ...',ncores));
	qrR<-qr(A, LAPACK = TRUE)
	Qm <- qr.Q(qrR)
	A_hat <- qr.R(qrR)
    t2=proc.time()
    if(verbose>1) print(sprintf('time in qr & A_hat:%.2f secs',(t2-t1)[3] ))
	# R=list()
	# for (i in 1:length(X)){ 
    i=NULL#Used to suppress NOTE while checking the package
    tmpR=foreach(i =1:length(X),.packages='Matrix') %dopar%{
		X_hat = Matrix::t(Qm) %*% (X[[i]]%*%Qm)
		# R[[i]] = unregularizedRUpdate(A_hat, X_hat)
        M1 = backsolve(A_hat, X_hat)
        Rn = Matrix::t(backsolve(A_hat, Matrix::t(M1)))        
	}
    # browser()
	return (tmpR)
}
# ------------------ Update Z ------------------------------------------------
 updateZ<-function(A, P, lmbdaZ){
    Z = list()
    if (length(P) == 0)
        return (Z);
    #_log.debug('Updating Z (Norm EQ, %d)' % len(P))
    # pinvAt = inv(dot(A.T, A) + lmbdaZ * eye(A.shape[1], dtype=A.dtype))
    # pinvAt = t(pinvAt %*% t(A))
    # for (i in 1:length(P)){
        # if issparse(P[[i]]){
            # Zn = P[i].tocoo().T.tocsr().dot(pinvAt).T
		# }else{
            # Zn = (pinvAt %*% P[[i]])
		# }
        # Z[[i]]=Zn
	# }
    return (Z)
}
##---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
compute_fit <-function(X, A, R,normX=NULL){#unused Params:, P, Z, lmbdaA, lmbdaR, lmbdaZ,
    # """Compute fit for full slices"""
    f = 0
    # precompute norms of X
    if(is.null(normX)){
        sumNorm=0
        for(p in 1:length(X)){
            sumNorm = sumNorm + sum(X[[p]]@x^2)#sum(X[[p]]@x==1)#RDF tensor [sum(M.data ** 2) for M in X]
        }
    }else{
        sumNorm = normX
       }

   for (i in 1:length(X)){
	    # print(sprintf('slice:%d',i))
		t1=proc.time()
        ARAt = A %*% (R[[i]]%*% Matrix::t(A))
		t2=proc.time()
		# print(sprintf('Time to calc ARAt:%f',(t2-t1)[3]))
		tmp=X[[i]] - ARAt
        rm(ARAt)#memory issues        
		t3=proc.time()
		# print(sprintf('Time to calc difference:%f',(t3-t2)[3]))
        f = f + (Matrix::norm(tmp,'f') ^ 2)
        rm(tmp)#memory issues        
   }
    if(get("verbose", parent.frame())>1) print(sprintf("f=%f",f))
    return (1 - f / sumNorm)
}


# ---------------------------------------------------------------------------

compute_fit_par <-function(X, A, R,normX=NULL,ncores,ChnkLen=10000){#unused Params:, P, Z, lmbdaA, lmbdaR, lmbdaZ,
    # """Compute fit for full slices"""
    #require(parallel)
	#require(doParallel)
	cluster <- makeCluster(ncores,outfile="comp_fit_par.log")
	doParallel::registerDoParallel(cluster)


    f = 0
    # precompute norms of X
    if(is.null(normX)){
        sumNorm=0
        for(p in 1:length(X)){
            sumNorm = sumNorm + sum(X[[p]]@x^2)#sum(X[[p]]@x==1)#RDF tensor [sum(M.data ** 2) for M in X]
        }
    }else{
        sumNorm = normX
       }
   A=as.matrix(A)
   for (i in 1:length(X)){
	    print(sprintf('slice:%d',i))
        RAt=(R[[i]]%*% Matrix::t(A))
        ones_ix=cbind(X[[i]]@i[X[[i]]@x==1],X[[i]]@j[X[[i]]@x==1])+1
        # for(ch in 1:ceiling(dim(A)[1]/ChnkLen)){
            f_ch=0
            ch=NULL#Used to suppress NOTE while checking the package
        tmp=foreach::foreach(ch =1:ceiling(dim(A)[1]/ChnkLen),.packages='Matrix', .combine="sum") %dopar%{
            print(paste("ch:",ch))
            st_ix=(1+(ch-1)*ChnkLen)
            end_ix=min(ch*ChnkLen,dim(A)[1])
            # ix=(1+(ch-1)*ChnkLen):min(ch*ChnkLen,dim(A)[1])
            ix=st_ix:end_ix
            t1=proc.time()
            ARAt_ch = A[ix,,drop=FALSE] %*% RAt
            t2=proc.time()
            print(sprintf('Time to calc ARAt:%f',(t2-t1)[3]))
            ones_ix_ch=ones_ix[ones_ix[,1]>=st_ix & ones_ix[,1]<=end_ix,,drop=FALSE]
            ones_ix_ch[,1]=ones_ix_ch[,1]-st_ix+1
            if(nrow(ones_ix_ch)>0){
               ones_tmp=2*sum(ARAt_ch[ones_ix_ch])
            }else{
                ones_tmp=0
            }
            totalErr=(Matrix::norm(ARAt_ch,'f') ^ 2)-ones_tmp#+cnt_ones
            rm(ARAt_ch)#memory issues        
            t3=proc.time()
            print(sprintf('Time to calc difference:%f',(t3-t2)[3]))          
            f_ch = totalErr
            # rm(tmp)#memory issues        
        }
        f_prd=tmp+sum(X[[i]]@x==1)
        print(paste("===============prd:",i," totalErr:",f_prd))
        f = f +  f_prd
    }
   print(sprintf("f=%f",f))
    return (1 - f / sumNorm)
}

# ---------------------------------------------------------------------------

compute_fval<-function(X, A, R, lambdaA, lambdaR, normX){
    print("Compute fit for full slices")
    if(lambdaA!=0){
		f = lambdaA * norm(A,"f") ^ 2;
	}else{ f=0;}
	# browser()
    for (i in 1:length(X)){
	    print(sprintf('slice:%d',i))
		t1=proc.time()
        ARAt = A %*% (R[[i]]%*% Matrix::t(A))
		t2=proc.time()
		print(sprintf('Time to calc ARAt:%f',(t2-t1)[3]))
		tmp=X[[i]] - ARAt
		t3=proc.time()
		print(sprintf('Time to calc difference:%f',(t3-t2)[3]))
        f = f + (Matrix::norm(tmp,'f') ^ 2) / normX[[i]] 
        # f = f + 0.5*(Matrix::norm(tmp,'f') ^ 2)  # Yago paper
		t4=proc.time()
		rm(tmp)#memory issues
		gc(T)
		print(sprintf('Time to calc norm:%f',(t4-t3)[3]))
		if(lambdaR!=0){
		 f = f + lambdaR * norm(R[[i]],'f') ^ 2
		}
	}
    return (f)
}
# ---------------------------------------------------------------------------
## 27/12/2019: avoid calculating difference
## error_in_ones sum((1-ARAt[@i,@j])^2)
## totalErr=norm(ARAt)-sum((ARAt[@i,@j])^2)+error_in_ones
##---------------------------------------------------------------------------
measureProgress <- function(Aold,Rold,Anew,Rnew,ErrThr=Inf){
     err=0
	 for (i in 1:length(Rnew)){
		err = err + norm(as.matrix(Rnew[[i]]-Rold[[i]]),'f')
        # if(err > ErrThr) return(err)
	 }
     if(err>ErrThr) return(err)
	 err=err+norm(as.matrix(Anew-Aold),'f')
	 return(err)
}


