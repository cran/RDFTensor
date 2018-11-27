# 23/5/2018
# translated from serial_parCube.m
# This function is translated from MATLAB code of MATLAB Tensor Toolbox

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

# I tried to keep the same names of variables (but F and options are reseved in R)

serial_parCube <- function(X,R,sample_factor,repetitions,opts=NULL){
# %Vagelis Papalexakis, 2012 - 2014
# %School of Computer Science, Carnegie Mellon University
# %Serial 3-mode ParCube algorithm for memory resident tensors
# %
# %The algorithm first appeared in:
# %   Papalexakis, Evangelos E., Christos Faloutsos, and Nicholas D. Sidiropoulos. 
# %   "Parcube: Sparse parallelizable tensor decompositions." Machine Learning and Knowledge Discovery in Databases. 
# %   Springer Berlin Heidelberg, 2012. 521-536.
# %
# %INPUTS:
# %   X: input tensor is "Tensor Toolbox" format, list(subs,vals,size)
# %   R: rank of the decomposition
# %   sample_factor: [s1 s2 s3] such that each sampled tensor is of size
# %       [I/s1 J/s2 K/s3]
# %   repetitions: number of sampling repetitions
# %   opts: structure that stores options of the algorithm. For default,
# %       leave blank or use 'default_parcube_options()'.
# %           opts.p: percentage of common indices
# %           opts.nonneg: nonnegativity constraint enforced (binary)
# %           opts.loss: loss function (opts: 'fro' for Frobenius norm
# %               'kl' for KL-divergence)
# %           opts.weights: function of calculating sampling weights
# %               (opts: 'sum_abs' for sum of absolute values or 'sum_squares'
# %               for sum of squares)

# %           opts.normalize: normalize the factor matrices to unit norm per column (binary);
# %           opts.tolerance: the numerical tolerance of the algorithm
# %               (everything smaller than that is considered zero)
# %           opts.internal_tolerance: the tolerance of the solvers
# %               used interally
# %           opts.num_restarts: number of repetitions of each decomposition of each sample
# %
# %OUTPUT:
# %   Fac_out: structure that holds a cell array with the factor matrices
# %       (Fac_out.U) and the scaling lambda (Fac_out.lambda)
# %
if (is.null(opts)){
   opts = default_parcube_options(); 
}

nonneg = opts$nonneg;
p = opts$p;
s = X$size; I = s[1]; J = s[2]; K = s[3];
A = Matrix::spMatrix(x=0,i=1,j=1,nrow=I,ncol=R);
B = Matrix::spMatrix(x=0,i=1,j=1,nrow=J,ncol=R);
C = Matrix::spMatrix(x=0,i=1,j=1,nrow=K,ncol=R);
lambda=matrix(0,R,repetitions)
for (i in 1:repetitions){
    if (i == 1){
        tmp = parCube_core(X,sample_factor,list(),opts);
		print('parCube_core')
		print(tmp)
		Xs=tmp[[1]]; idx_i=tmp[[2]]; idx_j=tmp[[3]]; idx_k=tmp[[4]]
        fixed_i = idx_i[1:ceiling(p*length(idx_i))];
        fixed_j = idx_j[1:ceiling(p*length(idx_j))];
        fixed_k = idx_k[1:ceiling(p*length(idx_k))];
		fixed=list()
        fixed[[1]] = fixed_i; fixed[[2]] = fixed_j; fixed[[3]]=fixed_k;
    }else{
		tmp = parCube_core(X,sample_factor,fixed,opts);
		Xs=tmp[[1]]; idx_i=tmp[[2]]; idx_j=tmp[[3]]; idx_k=tmp[[4]]
        # [Xs idx_i idx_j idx_k] = parCube_core(X,sample_factor,fixed,opts);
    }
    
    if (Xs$nnz > 0){
		normXs = spt_norm(Xs);

		all_loss = rep(0,opts$num_restarts);
		all_factors = list();
		for (num_restart in 1:opts$num_restarts){
			if (opts$loss== 'fro'){
				if(nonneg=='1'){
					opt_nmu=list();
					opt_nmu[['tol']] = opts$internal_tolerance;
					print('calc nmu...')
					tmp = cp_nmu(Xs,R,opt_nmu);
					factors=tmp[[1]]
				}else{
					# stop('cp_als: not yet implemented')
                    print('calc als...')
					tmp_factors = cp_als(Xs,R,list(tol=opts$internal_tolerance));
                    factors=tmp_factors[[1]]
				}
			}else{
				  if(opts$loss=='kl'){
					# stop('cp_apr: not yet implemented')
                    print('calc apr ...')
					tmp_factors = cp_apr(Xs,R,list(tol=opts$internal_tolerance));
                    factors=tmp_factors[[1]]
                    }
			}
			 all_loss[num_restart] = normXs^2 + ktensor_norm(factors)^2 - 2 * tt_innerprod(factors,Xs);
			 all_factors[[num_restart]] = factors;
		}
		min_idx = which.min(all_loss);
		factors = all_factors[[min_idx]];
		As=factors$u[[1]];Bs=factors$u[[2]];Cs=factors$u[[3]];#%lambda = factors.lambda;As = As*diag(lambda);
		
		lambda[,i] = factors$lambda;
		As = As%*%diag(lambda[,i]^(1/3)); Bs = Bs%*%diag(lambda[,i]^(1/3)); Cs = Cs%*%diag(lambda[,i]^(1/3));  
    }else{
       As = 0; Bs = 0; Cs = 0; 
    }
    
    Atemp = Matrix::spMatrix(x=0,i=1,j=1,nrow=I,ncol=R);
	Btemp = Matrix::spMatrix(x=0,i=1,j=1,nrow=J,ncol=R);
	Ctemp = Matrix::spMatrix(x=0,i=1,j=1,nrow=K,ncol=R);
    # Atemp = sparse(I,R); Btemp = sparse(J,R);Ctemp = sparse(K,R);
    Atemp[idx_i,] = As; Btemp[idx_j,] = Bs; Ctemp[idx_k,] = Cs;
    allCorrespondences = rep(0,R)#zeros(R,1);
    
     # %now, normalize the common part of every column
	 eps=2.2e-16
    for (r in 1:R){
       norm_a = norm(Atemp[fixed_i,r,drop=FALSE],'f');Atemp[,r] = Atemp[,r]/(norm_a+eps); lambda[r,i] = norm_a;
       norm_b = norm(Btemp[fixed_j,r,drop=FALSE],'f');Btemp[,r] = Btemp[,r]/(norm_b+eps); lambda[r,i] = lambda[r,i]*norm_b;
       norm_c = norm(Ctemp[fixed_k,r,drop=FALSE],'f');Ctemp[,r] = Ctemp[,r]/(norm_c+eps); lambda[r,i] = lambda[r,i]*norm_c;
    }
    
    # %Do the merge
    if (i ==1){
        A = Atemp; B = Btemp; C = Ctemp;
    }else{
        valA=valB=valC=rep(0,R);
        for (f1 in 1:R){
            for (f2 in 1:R){
                valA[f2] = Matrix::t(Atemp[fixed_i,f1]) %*% A[fixed_i,f2];
                valB[f2] = Matrix::t(Btemp[fixed_j,f1]) %*% B[fixed_j,f2];
                valC[f2] = Matrix::t(Ctemp[fixed_k,f1]) %*% C[fixed_k,f2];
            }
            idx = which.max(valA);
            mask2 = which(A[,idx]==0);
            A[mask2,idx] = A[mask2,idx] + Atemp[mask2,f1];#%Update ONLY the zero values
            allCorrespondences[f1] = idx;
            
# % % % %      idx = max(valB);
            mask2 = which(B[,idx]==0);
            B[mask2,idx] = B[mask2,idx] + Btemp[mask2,f1];#%Update ONLY the zero values
            
# % % % %             [junk idx] = max(valC);
            mask2 = which(C[,idx]==0);
            C[mask2,idx] = C[mask2,idx] + Ctemp[mask2,f1];#%Update ONLY the zero values
        }
    lambda[,i] = lambda[allCorrespondences,i];
    }
}
print(lambda)
newlambda = rowMeans(lambda)#mean(lambda,2);
print(newlambda)
A = as(A%*%diag(newlambda),'sparseMatrix');
lambda = rep(1,R)#ones(R,1);

if (opts$normalize == 1){
    for (i in 1:R){
        normA = norm(A[,i,drop=FALSE],'f');
        normB = norm(B[,i,drop=FALSE],'f');
        normC = norm(C[,i,drop=FALSE],'f');
        lambda[i] = normA*normB*normC;
        A[,i] = A[,i]/normA;
        B[,i] = B[,i]/normB;
        C[,i] = C[,i]/normC;
    }
    tmp = sort(lambda,decreasing=TRUE,index.return=TRUE);
	lambda=tmp$x;lam_idx=tmp$ix
    A = A[,lam_idx];
    B = B[,lam_idx];
    C = C[,lam_idx];
}


	A[abs(A)<opts$tolerance] = 0;
	B[abs(B)<opts$tolerance] = 0;
	C[abs(C)<opts$tolerance] = 0;
	Fac_out = list(lambda=lambda,u=list(A=A, B=B, C=C));

return(Fac_out)
}

###---------------------------------------------------------------##

 parCube_core <- function(X,sample_factor,fixed_set,opts){
 # [Xs idx_i idx_j idx_k Xma Xmb Xmc] 
# %Vagelis Papalexakis, 2012
# %School of Computer Science, Carnegie Mellon University
# %Core sampling function for the ParCube algorithm, for memory resident
# %tensors.
	if (length(sample_factor)>1){
		s1 = sample_factor[1]; 
		s2 = sample_factor[2];
		s3 = sample_factor[3];
	}else{
		s1 = sample_factor; 
		s2 = sample_factor;
		s3 = sample_factor;
	}

	if (length(fixed_set)!=0){
	   fixed_i = fixed_set[[1]];
	   fixed_j = fixed_set[[2]];
	   fixed_k = fixed_set[[3]];
	}else{
		fixed_i = list();
		fixed_j = list();
		fixed_k = list();
	}

	if (opts$weights== 'sum_abs'){
		negative_values = which(X$vals<0);#not needed for RDFTensor
		if (length(negative_values)!=0){
			stop("Tensor contains negative values, not yet implemented.")
			# X_abs = double(X);
			# X_abs = abs(X_abs);
			# X_abs = sptensor(X_abs);
		}else{
		   X_abs = X; 
		}
	}else{
		 if(opts$weights== 'sum_squares'){
		# X_abs = sptensor(tensor(X).^2);
			# % ## all 1's   RDFTensor
			X_abs=X
		}
	}
	s = X$size; I = s[1]; J = s[2]; K = s[3];

	Xma = spt_collapse(X_abs,2:3);#returns a vector
	Xmb = spt_collapse(X_abs,c(1,3));
	Xmc = spt_collapse(X_abs,1:2);

	idx_i = which(Xma!=0);
	idx_i = randsample_weighted(idx_i,I/s1,Xma[idx_i]);#full: not needed, already dense vector
	if (length(fixed_i)!=0){
		idx_i = union(idx_i,fixed_i);
	}

	idx_j = which(Xmb!=0);
	idx_j = randsample_weighted(idx_j,J/s2,Xmb[idx_j]);
	if (length(fixed_j)!=0){
		idx_j = union(idx_j,fixed_j);
	}

	idx_k = which(Xmc!=0);
	idx_k = randsample_weighted(idx_k,K/s3,Xmc[idx_k]);
	if (length(fixed_k)!=0){
		idx_k = union(idx_k,fixed_k);
	}

	# Xs = X(idx_i,idx_j,idx_k);
	newsubs_i=match(X$subs[,1], idx_i); 
	newsubs_j=match(X$subs[,2], idx_j); 
	newsubs_k=match(X$subs[,3], idx_k); 
	# & X$subs[,2] %in% idx_j & X$subs[,3] %in% idx_k
	newsubs=cbind(newsubs_i,newsubs_j,newsubs_k)[(!is.na(newsubs_i)) & (!is.na(newsubs_j)) & (!is.na(newsubs_k)),]
	# newvals=X$vals[X$subs[,1] %in% idx_i & X$subs[,2] %in% idx_j & X$subs[,3] %in% idx_k]
	newvals=X$vals[(!is.na(newsubs_i)) & (!is.na(newsubs_j)) & (!is.na(newsubs_k))]
	newsiz=c(length(idx_i),length(idx_j),length(idx_k))
	Xs=list(subs=newsubs,vals=newvals,size=newsiz,nnz=length(newvals))
	# print(Xs)
	print("Size of Xs:")
	print(newsiz)
	return(list(Xs=Xs,idx_i=idx_i,idx_j=idx_j,idx_k=idx_k))
}

###
default_parcube_options<-function(){
	# %%Vagelis Papalexakis, 2012 - 2014
	# %School of Computer Science, Carnegie Mellon University
	# %
	# %Default options for the ParCube algorithm
	# %
	# %The algorithm first appeared in:
	# %   Papalexakis, Evangelos E., Christos Faloutsos, and Nicholas D. Sidiropoulos. 
	# %   "Parcube: Sparse parallelizable tensor decompositions." Machine Learning and Knowledge Discovery in Databases. 
	# %   Springer Berlin Heidelberg, 2012. 521-536.
	# %
	opts=list()
	opts[['p']] = 0.35;
	opts[['nonneg']] = 1;
	opts[['loss']] = 'fro';
	# % opts[['loss']] = 'kl';
	# % opts[['weights']] = 'sum_abs';
	opts[['weights']] = 'sum_squares';
	opts[['MAXMEM']] = 8192;
	opts[['parallel_merge']] = 0;
	opts[['normalize']] = 1;
	opts[['tolerance']] = 10^-5;#10^-8
	opts[['internal_tolerance']] = 10^-4;#10^-6
	opts[['num_restarts']] = 1;
	
	return(opts)
}

spt_norm<-function(X){
	return(sqrt(sum(X$vals^2)))
}

###--------------------------------------------------##

 spt_collapse <- function(tn,dims=NULL,fun='sum'){
# %COLLAPSE Collapse sparse tensor along specified dimensions.
# %
# %  S = COLLAPSE(T,DIMS) sums the entries of T along all dimensions
# %  specified in DIMS. If DIMS is negative, then T is summed across
# %  all dimensions *not* specified by -DIMS.
# %
# %  S = COLLAPSE(T) is shorthand for S = COLLAPSE(T,1:ndims(T)).
# %
# %  S = COLLAPSE(T,DIMS,FUN) accumulates the entries of T using the
# %  accumulation function @FUN.
# %
# %  Examples
# %  subs = [1 1 1; 1 1 3; 2 2 4; 4 4 4]
# %  vals = [10.5; 1.5; 2.5; 3.5]
# %  X = sptensor(subs,vals,[4 4 4]);
# %  Y = collapse(X,[2 3]) %<-- sum of entries in each mode-1 slice
# %  Y = collapse(ones(X),[1 2]) %<-- nnz in each mode-3 slide
# %  Y = collapse(ones(X),[1 2],@max) %<-- 1 if mode-3 has any entry
# %  Y = collapse(ones(X),-3,@max); %<-- equivalent
# %
# %  See also SPTENSOR, SPTENSOR/SCALE.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.


	if (is.null(dims)){
		dims = 1:ndims(tn);
	}

	dims = tt_dimscheck(dims,ndims(tn));
	remdims = (1:ndims(tn))[-dims];

	# % Check for the case where we accumulate over *all* dimensions
	if (length(remdims)==0){
		s = eval(parse(text=paste(fun,'(tn$vals)')));
		return(s);
	}

	# % Calculate the size of the result
	newsiz = tn$size[remdims];
	newsubs = tn$subs[,remdims];

	# % Check for the case where the result is just a dense vector
	if (length(remdims) == 1){
		# if (tn$nnz > 0){
		if (nrow(tn$subs) > 0){
			# s = accumarray(t.subs(:,remdims), t.vals, [newsiz 1], fun);
			s =  eval(parse(text=paste('pracma::accumarray(newsubs,tn$vals,sz=newsiz,func=',fun,')')));
		}else{
			s = rep(0,newsiz)#zeros(newsiz,1);
		}
		return(s);
	}

	# % Create the result
	if (nrow(tn$subs) > 0){
		s =  eval(parse(text=paste('pracma::accumarray(newsubs,tn$vals,sz=newsiz,func=',fun,')')));
	  # s = list(subs=newsubs, vals=tn$vals, newsiz, fun);
	}else{
	  s = list(subs=cbind(1,1,1)[1==2,],vals=numeric(0),size=newsiz,nnz=0);#sptensor
	}
	return(s);
}

###--------------------------------------###

randsample_weighted<-function(p,n,w){
# %Vagelis Papalexakis, 2012
# %School of Computer Science, Carnegie Mellon University
n = ceiling(n);
nnz_w = sum(w!=0);
n = min(n, nnz_w);#%this guard makes sure that we always have data to sample from
s = rep(0,n);

	for (i in 1:n){
			# if(strfind(version('-release'),'2009'))
			# %http://www.mathworks.com/matlabcentral/newsreader/view_thread/255206
			# if size(w,1) == 1
				# w = w';
			# end
			# pp = w/sum(w);
			# temp = [0;cumsum(pp)];
			# e = min(temp,1);
			# e(end) = 1;
			# pp = diff(e);
			# w = pp;
		# end
			s[i] = sample(p,1,replace=TRUE,prob=w);
			w[p==s[i]]=0;
	}
	return(s)
}