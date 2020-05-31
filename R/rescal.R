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
rescal <- function(X, rnk,ainit = 'nvecs',verbose=2,Ainit=NULL,Rinit=NULL,lambdaA=0,	lambdaR=0, lambdaV=0,
	epsilon=1e-3,maxIter=100, minIter=1, P = list(),orthogonalize=FALSE,func_compute_fval='compute_fit'){
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

	n=nrow(X[[1]])
	if(!is.null(Ainit)){
        if(ncol(Ainit)!=rnk){
            stop(sprintf("Ainit number of columns (%d) must equal rank (%d)",ncol(Ainit),rnk))
        }
        if(nrow(Ainit)!=n){
            stop(sprintf("Ainit number of rows (%d) must equal n (%d)",nrow(Ainit),n))
        }
		A=Ainit;
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
				for (i in 1:k){
					S = S + X[[i]]
					S = S + Matrix::t(X[[i]])
				}
				if(verbose>0) print("Calculating  eigen vectors...")
				tt=eigen(S)#NB: need only rnk vectors but no option to specify it in R eigen
				A = tt$vectors[,1:rnk]
                if(verbose>2 && rnk<=10){
                    print("initial A ................")
                    print(A)
                }
			}
			else{
				stop(sprintf('Unknown init option:%s, should be random or nvecs' ,ainit))
			}
		}
	}
    # ------- initialize R and Z ---------------------------------------------
    if(verbose>0)print("initialize R and Z...")
	if(!is.null(Rinit)){
		R=Rinit
	}else{
		t1=proc.time()
		R = updateR(X, A, lambdaR,0)
        # print(R[[1]])
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
    if(verbose>1)print(sprintf("sumNorm:%f",sumNorm))
    #  ------ compute factorization ------------------------------------------
    if(verbose>1) print("compute factorization...")
	fit = fitchange = fitold = 0
    exectimes = NULL
	all_err=NULL #maintain a list of errors in every iteration 
    for (itr in 1:maxIter){
		 if(verbose>0) print(sprintf("-----------------------------iteration: %d ----------------------------------",itr))
        tic = proc.time()
        fitold = fit
		Aold=A
        A = updateA(X, A, R, P, Z, lambdaA, orthogonalize)
			if(verbose>2){print(A)}
		tUA= proc.time()
		if(verbose>1) print(sprintf('updateA in:%.3f sec',(tUA-tic)[3]))
		if(verbose>1) print("Update R...")
        Rold=R
		R = updateR(X, A, lambdaR,verbose)
		tUR= proc.time()
		if(verbose>1) print(sprintf('updateR in:%.3f sec',(tUR-tUA)[3]))
			if(verbose>2){print(R)}
        Z = updateZ(A, P, lambdaV)

        # compute fit value!!!
		# if(itr%%fent==0){
			if (!is.null(func_compute_fval )){
                if(func_compute_fval=='measureProgress'){
                    fit=measureProgress(Aold,Rold,A,R)
                }else{
                  fit=compute_fit(X, A, R ,normX=sumNorm)#, lambdaA, lambdaR,
                }
                    # print(fit)
                    if(verbose>1) print(sprintf('Calc fit in:%.3f sec',(proc.time()-tUR)[3]))
                    fitchange = abs(fitold - fit)
			}else{
				fit = Inf
			}
		# }

        toc = proc.time()
        exectimes=c(exectimes,toc[3] - tic[3])
		all_err=rbind(all_err,cbind(itr,fit=fit,delta=fitchange,Ex_time=exectimes[itr]))
        if(verbose>0) print(sprintf('[%3d] fit: %0.5f | delta: %7.1e | secs: %.5f', itr, fit, fitchange, exectimes[itr])
        )
        if (itr > minIter && fitchange < epsilon){
            break;
		}
	}
    return (list(A=A, R=R, all_err, nitr=itr , times=as.vector(exectimes)))
}

##--------------------------update A--------------------------------------------------
updateA <- function(X, A, R, P, Z, lambdaA, orthogonalize){
    # """Update step for A"""
    n=nrow(A); rnk = ncol(A);
    F = Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(n, rnk))#zeros((n, rnk), dtype=A.dtype)
    E = Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(rnk, rnk))#zeros((rnk, rnk), dtype=A.dtype)

    AtA = Matrix::t(A)%*% A#dot(A.T, A)
	# browser()
    for (i in 1:length(X)){
        F = F + X[[i]]%*% (A %*% Matrix::t(R[[i]])) + Matrix::t(X[[i]]) %*% ((A %*% R[[i]]))
        E = E + (R[[i]] %*% (AtA %*% Matrix::t(R[[i]]))) + (Matrix::t(R[[i]]) %*% (AtA %*% R[[i]]))
    }
    # regularization:zeros
    I = lambdaA * diag(rnk)

    # attributes
	if(length(Z)>0){
		for( i in 1:length(Z)){
			F = F + (P[[i]] %*% (Matrix::t(Z[[i]])))
			E = E + (Z[[i]] %*% (Matrix::t(Z[[i]])))#dot(Z[i], Z[i].T)
		}
	}
    # finally compute update for A  O(r^3)
    A = Matrix::t(Matrix::solve(I + Matrix::t(E), Matrix::t(F)))
    if(orthogonalize) {
        stop("orth(A) not implemented.")
        # return(orth(A));
	 }else{
        return(A)
     }
}

##---------------------------------------------------------------------------

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
		if(get("verbose", parent.frame())>1)print(sprintf("i=%d, Err=%f, t=%.3f",i,err,(t1-t0)[3]))
        rm(tmp)#memory issues
		gc()
	}
	return (err)
}

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
		t3=proc.time()
		# print(sprintf('Time to calc difference:%f',(t3-t2)[3]))
        f = f + (Matrix::norm(tmp,'f') ^ 2) 
   }
   if(get("verbose", parent.frame())>1) print(sprintf("f=%f",f))
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

measureProgress <- function(Aold,Rold,Anew,Rnew){
	 err=norm(as.matrix(Anew-Aold),'f')
	 for (i in 1:length(Rnew)){
		err = err + norm(as.matrix(Rnew[[i]]-Rold[[i]]),'f')
	 }
	 return(err)
}


# ------------------ Update R ------------------------------------------------
 updateR <- function(X, A, lambdaR,verbose=1){
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
# ------------------ Update R Unreg using QR------------------------------------------------
 RUpdate <- function(X, A, R){
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
		R[[i]] = unregularizedRUpdate(A_hat, X_hat)
	}
	return (R)
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
