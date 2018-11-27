# 21/5/2018
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

cp_nmu<-function(X,R,opts=list()){
# %CP_NMU Compute nonnegative CP with multiplicative updates.
# translated from cp_nmu.m : MATLAB Tensor Toolbox
# X is assumed to be RDFTensor: Boolean sparse, [subs,vals,size]#list of frontal slices
# R: rank
# %   P = CP_NMU(X,R) computes an estimate of the best rank-R PARAFAC
# %   model of a tensor X with nonnegative constraints on the factors.
# %   This version uses the Lee & Seung multiplicative updates from
# %   their NMF algorithm.  The input X can be a tensor, sptensor,
# %   ktensor, or ttensor. The result P is a ktensor.
# %
# %   P = CP_NMU(X,R,OPTS) specify options:
# %   OPTS.tol: Tolerance on difference in fit {1.0e-4}
# %   OPTS.maxiters: Maximum number of iterations {50}
# %   OPTS.dimorder: Order to loop through dimensions {1:ndims(A)}
# %   OPTS.init: Initial guess [{'random'}|'nvecs'|cell array] nvecs not avail.
# %   OPTS.printitn: Print fit every n iterations {1}
# %
# %   [P,U0] = CP_NMU(...) also returns the initial guess.
# %
# %   Examples:
# %   X = sptenrand([5 4 3], 10);
# %   P = cp_nmu(X,2);
# %   P = cp_nmu(X,2,struct('dimorder',[3 2 1]));
# %   P = cp_nmu(X,2,struct('dimorder',[3 2 1],'init','nvecs'));
# %   U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of P
# %   P = cp_nmu(X,2,struct('dimorder',[3 2 1],'init',{U0}));
# %
# %   See also KTENSOR, TENSOR, SPTENSOR, TTENSOR.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.

# %% Extract number of dimensions and norm of X.
# N = ndims(X);
N=3#list(n=nrow(X[[1]]),m=ncol(X[[1]]),l=length(X))
# normX = norm(X);
# normX = lapply(X,function(M){(norm(as.matrix(M@x),'f')^2)})
# normX=sqrt(sum(unlist(lapply(X,function(M) {sum(M)}))))
normX=sqrt(sum(X$vals^2));#!! no need for ^2 in RDFTensor
sizeX=X$size;#c(n=nrow(X[[1]]),l=length(X),m=ncol(X[[1]]))

# %% Set algorithm parameters from input or by using defaults
fitchangetol = ifelse(is.null(opts[['tol']]),1e-4,opts[['tol']]);
maxiters = ifelse(is.null(opts[['maxiters']]),100,opts[['maxiters']]);
if(is.null(opts[['dimorder']])){
	dimorder = 1:N;
}else{
	dimorder=opts[['dimorder']];
}

init = ifelse(is.null(opts[['init']]),'random',opts[['init']]);
printitn = ifelse(is.null(opts[['printitn']]),1,opts[['printitn']]);
epsilon = 1e-12;  #% Small number to protect against round-off error

# %% Error checking 
# % Error checking on maxiters
if (maxiters < 0)
    stop('OPTS.maxiters must be positive');

print(dimorder)
# % Error checking on dimorder
if(!identical(1:N,sort(dimorder)))
    stop('OPTS.dimorder must include all elements from 1 to ndims(X)');


# %% Set up and error checking on initial guess for U.
# if iscell(init)
    # Uinit = init;
    # if numel(Uinit) ~= N
        # error('OPTS.init does not have %d cells',N);
    # end
    # for n = dimorder(1:end);
        # if ~isequal(size(Uinit{n}),[size(X,n) R])
            # error('OPTS.init{%d} is the wrong size',n);
        # end
    # end
# else
if (is(init)[1]=="list"){
    Uinit = init;
    if(length(Uinit)!=N){
        stop(sprintf('OPTS.init does not have %d cells',N))
    }    
}else{
   if (init=='random'){
        Uinit = list(N);#cell(N,1);
        for (n in  dimorder){
            Uinit[[n]] = matrix(runif(sizeX[n]*R), ncol=R,byrow=TRUE)#rand(size(X,n),R) + 0.1;
        }
    }else{
			if (init=='nvecs' || init=='eigs'){ 
				Uinit = list(N);#cell(N,1);
				for (n in dimorder){
					k = min(R,sizeX[n]-2);
					print(sprintf('  Computing %d leading e-vectors for factor %d.',k,n));
					# Uinit[[n]] = abs(nvecs(X,n,k));##!!!Check
					print("Calculating  eigen vectors...")
					# tt=eigen(S)#NB: need only rnk vectors but no option to specify it in R eigen
					# A = tt$vectors[,1:R]
					if (k < R){#difficult to be achieved
					  Uinit[[n]] = rBind(Uinit[[n]],matrix(runif(sizeX[n]*(R-k)), ncol=R-k,byrow=TRUE))#[Uinit{n} rand(size(X,n),R-k)]; 
					}
				}
			}else{
				stop('The selected initialization method is not supported');

				}
		}
}

# %% Set up for iterations - initializing U and the fit.
U = Uinit;
fit = 0;
stats=NULL

	if (printitn>0)#verbose
	  print('Nonnegative PARAFAC:');


	# %% Main Loop: Iterate until convergence
	for (iter in 1:maxiters){
		t0=proc.time()
		fitold = fit;

		# % Iterate over all N modes of the tensor
		for (n in dimorder){

			# % Compute the matrix of coefficients for linear system
			Y = matrix(1,R,R)#ones(R,R);
			# for i = [1:n-1,n+1:N]
			for(i in (1:N)[-n]){
				# Y = Y .* (U{i}'*U{i});
				Y = Y * (Matrix::t(U[[i]])%*% U[[i]]);
			}
			Y = U[[n]] %*% Y;

			# % Initialize matrix of unknowns
			Unew = U[[n]];
			# print(sprintf("n:%d, Unew:",n))
			# print(Unew)
			# % Calculate Unew = X_(n) * khatrirao(all U except n, 'r').
			if (printitn>0) print(sprintf("n:%d, calc spt_mttkrp...",n))
			t1=proc.time()
			tmp = spt_mttkrp(X,U,n) + epsilon;
			t2=proc.time()
			# print(tmp)
			# % Update unknowns
			Unew = Unew * tmp;
			Unew = Unew / (Y + epsilon);

			U[[n]] = Unew;
		}

		# P = ktensor(U); assume lamda's as ones
		P=list(lambda=rep(1,R),u=U)
		#print(sprintf("normX:%f, norm:%f, innprod:%f",normX, ktensor_norm(P)^2,tt_innerprod(P,X)))
		normresidual = sqrt( normX^2 + ktensor_norm(P)^2 - 2 * tt_innerprod(P,X) );
		fit = 1 - (normresidual / normX); #%fraction explained by model
		fitchange = abs(fitold - fit);
		t3=proc.time()
		if ((iter%%printitn)==0){#mod(iter,printitn)==0
		  print(sprintf(' Iter %2d: fit = %e fitdelta = %7.1e, itime=%.2f', iter, fit, fitchange,itime=(t3-t0)[3]));
		}
		
		stats=rbind(stats,cbind(iter=iter,fit=fit,fitchange=fitchange,itime=(t3-t0)[3],mttkrp_time=(t2-t1)[3]))
		# % Check for convergence
		if ((iter > 1) && (fitchange < fitchangetol)){
			break;
		}

	}

# %% Clean up final result
# % Arrange the final tensor so that the columns are normalized.
	if (printitn>0) print("arranging results...")
	P = arrange(P);

	if (printitn>0){
	  normresidual = sqrt( normX^2 + ktensor_norm(P)^2 - 2 * tt_innerprod(P,X) );
	  fit = 1 - (normresidual / normX); #%fraction explained by model
	  print(sprintf(' Final fit = %e ', fit));
	}

	return(list(P=P,Uinit=Uinit,stats=stats,fit=fit));
}

###--------------------------------------------------------------------------------##

spt_mttkrp<-function(X,U,n){
# %MTTKRP Matricized tensor times Khatri-Rao product for sparse tensor.
# %
# %   V = MTTKRP(X,U,n) efficiently calculates the matrix product of the
# %   n-mode matricization of X with the Khatri-Rao product of all
# %   entries in U, a cell array of matrices, except the nth.  How to
# %   most efficiently do this computation depends on the type of tensor
# %   involved.
# %
# %   See also SPTENSOR, TENSOR/MTTKRP, SPTENSOR/TTV
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.

# % In the sparse case, it is most efficient to do a series of TTV operations
# % rather than forming the Khatri-Rao product.

	N = ndims(X);

	if (n==1){
		R = ncol(U[[2]])#size(U{2},2);
	}else{
		R = ncol(U[[1]])#size(U{1},2);
	}
# if(n==2) browser()
	V = matrix(0,nrow=X$size[n],ncol=R);
	for(r in 1:R){
		# % Set up cell array with appropriate vectors for ttv multiplication
		Z = list(N);
		for (i in (1:N)[-n]){#[1:n-1,n+1:N]
			Z[[i]] = U[[i]][,r]#(:,r);
		}
		# % Perform ttv multiplication
		V[,r] = as.numeric(spt_ttv(X, Z, -n));
	}

return(V)
}

##---------------------------------------------------------------##

kt_mttkrp<- function(X,U,n){
# %MTTKRP Matricized tensor times Khatri-Rao product for ktensor.
# %
# %   V = MTTKRP(X,U,n) efficiently calculates the matrix product of the
# %   n-mode matricization of X with the Khatri-Rao product of all
# %   entries in U, a cell array of matrices, except the nth.  How to
# %   most efficiently do this computation depends on the type of tensor
# %   involved.
# %
# %   See also KTENSOR, KTENSOR/TTV
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.

	N = ndims(X);

	if (n==1){
		R = ncol(U[[2]])#size(U{2},2);
	}else{
		R = ncol(U[[1]])#size(U{1},2);
	}

	# % Compute matrix of weights
	# W = repmat(X.lambda,1,R);
	W = kronecker(matrix(1,1,R),X$lambda)
	# for i = [1:n-1,n+1:N]
	  # W = W .* (X.u{i}' * U{i});
		for(i in (1:N)[-n]){
            W = W * (Matrix::t(X$u[[i]])%*% U[[i]]);
        }

	# % Find each column of answer by multiplying columns of X.u{n} with weights 
	V = X$u[[n]] %*% W;
	return(V)
}

###-------------------------------------------------------------------------------------------##

arrange<-function(X,foo=NULL){
# %ARRANGE Arranges the rank-1 components of a ktensor.
# %
# %   ARRANGE(X) normalizes the columns of the factor matrices and then sorts
# %   the ktensor components by magnitude, greatest to least.
# %
# %   ARRANGE(X,N) absorbs the weights into the Nth factor matrix instead of
# %   lambda. 
# %
# %   ARRANGE(X,P) rearranges the components of X according to the
# %   permutation P. P should be a permutation of 1 to NCOMPOMENTS(X). 
# %
# %   See also KTENSOR, NCOMPONENTS, NORMALIZE.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.


# %% Just rearrange and return if second argument is a permutation
if (!is.null(foo) && length(foo)>1){
    X[['lambda']] = X$lambda[foo];
    for (i in 1 : 3){#ndims(X)
        X$u[[i]] = X$u[[i]][,foo]#(:,foo);
    }
    return(X);
}

# %% Ensure that matrices are normalized

X = normalize(X);

# %% Sort
# [X.lambda, idx] = sort(X.lambda, 1, 'descend');
	tmp = sort(X$lambda, decreasing=TRUE, index.return=TRUE);
	idx=tmp$ix
	X$lambda=tmp$x;
	for(i in 1 : ndims(X)){
		X$u[[i]] = X$u[[i]][,idx];
	}

# %% Absorb the weight into one factor, if requested
# if (!is.null(foo)){
    # r = length(X$lambda);
    # % X.u{end} = full(X.u{end} * spdiags(X.lambda,0,r,r));
    # X.u[[length(X.u)]] = full(X.u{end} * spdiags(X.lambda,0,r,r));
    # X$lambda = rep(1,length(X$lambda))#ones(size(X.lambda));
# }

	return(X)
}
###------------------------------------------------------------------##
ktensor_norm<-function(A){
# %NORM Frobenius norm of a ktensor(list(lambda,u)).
# %
# %   NORM(T) returns the Frobenius norm of a ktensor.
# %
# %   See also KTENSOR.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.

# % Retrieve the factors of A
	U = A[['u']];#A.u

	# % Compute the matrix of correlation coefficients
	coefMatrix = A[['lambda']] %*% Matrix::t(as.matrix(A[['lambda']]));
	for (i in 1:ndims(A)){
	  coefMatrix = coefMatrix * (Matrix::t(U[[i]])%*%U[[i]]);
	}

	nrm = sqrt(abs(sum(coefMatrix)));

	return(nrm);
}

###------------------------------------##
ndims<-function(X){#RDF Tensor
	return(3)
}

###--------------------------------------------------------------##

# function X = normalize(X,N,normtype,pmode)
normalize <-function(X,N=-1,normtype='f',pmode=NULL){ 
#normtype 2 is Fobinious norm
# no mode and N=-1
# %NORMALIZE Normalizes the columns of the factor matrices.
# %
# %   NORMALIZE(X) normalizes the columns of each factor matrix using the
# %   vector 2-norm, absorbing the excess weight into lambda. Also ensures
# %   that lambda is positive.   
# %
# %   NORMALIZE(X,N) absorbs the weights into the Nth factor matrix instead
# %   of lambda. (All the lambda values are 1.)
# %
# %   NORMALIZE(X,0) equally divides the weights across the factor matrices.
# %   (All the lambda values are 1.)
# %
# %   NORMALIZE(X,[]) is equivalent to NORMALIZE(X). 
# %
# %   NORMALIZE(X,'sort') is the same as the above except it sorts the
# %   components by lambda value, from greatest to least. 
# %
# %   NORMALIZE(X,V,1) normalizes using the vector one norm (sum(abs(x))
# %   rather than the two norm (sqrt(sum(x.^2))), where V can be any of the
# %   second arguments decribed above.
# %
# %   NORMALIZE(X,[],1,I) just normalizes the I-th factor using whatever norm
# %   is specified by the 3rd argument (1 or 2).
# %
# %   See also KTENSOR, ARRANGE, REDISTRIBUTE, TOCELL.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.


	if (N=='sort')    N = -2;

	if(!is.null(pmode)){
		for (r in 1:length(X$lambda)){
			# tmp = norm(X.u[[pmode]][,r],normtype);
			if(normtype=='1'){
				tmp = sum(X$u[[pmode]][,r])
			}else{
				 tmp = norm(as.matrix(X$u[[pmode]][,r]),normtype);
			}
			if (tmp > 0){
				X$u[[pmode]][,r] = X$u[[pmode]][,r] / tmp;
			}
     		X$lambda[r] = X$lambda[r] * tmp; 
		}
		return(X);
	}

# %% Ensure that matrices are normalized
	for (r in 1:length(X$lambda)){
		for (n in 1:ndims(X)){
			if(normtype=='1'){
				tmp = sum(X$u[[n]][,r])#norm(as.matrix(X$u[[n]][,r]),normtype);
			}else{
				 tmp = norm(as.matrix(X$u[[n]][,r]),normtype);
			}
				
			if (tmp > 0){
				X$u[[n]][,r] = X$u[[n]][,r] / tmp;
			}
     		X$lambda[r] = X$lambda[r] * tmp;        
		}
	}

# %% Check that all the lambda values are positive
	idx = which(X$lambda < 0);
	X$u[[1]][,idx] = -1 * X$u[[1]][,idx];
	X$lambda[idx] = -1 * X$lambda[idx];

# %% Absorb the weight into one factor, if requested
# if (N == 0)
    # D = diag(nthroot(X.lambda,ndims(X)));
    # X.u = cellfun(@(x) x*D, X.u, 'UniformOutput', false);
    # X.lambda = ones(size(X.lambda));
# else
if (N > 0){
    X$u[[N]] = X$u[[N]] %*% diag(X$lambda);
    X$lambda = rep(1,length(X$lambda))#ones(size(X.lambda));
}
# elseif (N == -2)
    # if ncomponents(X) > 1
        # [~,p] = sort(X.lambda,'descend');
        # X = arrange(X,p);
    # end
# end
 return(X)
}

###--------------------------------------------------------------------------##

kt_innerprod<-function(X,Y){
# %INNERPROD Efficient inner product with a ktensor.
# %
# %   R = INNERPROD(X,Y) efficiently computes the inner product between
# %   two tensors X and Y.  If Y is a ktensor, the inner product is
# %   computed using inner products of the factor matrices, X{i}'*Y{i}.
# %   Otherwise, the inner product is computed using ttv with all of
# %   the columns of X's factor matrices, X{i}.
# %
# %   See also KTENSOR, KTENSOR/TTV
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.

 M = X$lambda %*% t(Y$lambda);
    for( n in 1:ndims(X)){
        M = M * (Matrix::t(X$u[[n]]) %*% Y$u[[n]]);
    }
    res = sum(M);
	
   return(res)
}

tt_innerprod<-function(X,Y){
#X ktensor, Y tensor   
  # case {'tensor','sptensor','ttensor'}
    R = length(X$lambda);
    vecs = list(ndims(X));
    res = 0;
    for (r in 1:R){
      for( n in 1:ndims(X)){
        vecs[[n]] = X$u[[n]][,r]#(:,r);
      }
      res = res + X$lambda[[r]] * spt_ttv(Y,vecs);
    }
    
  # otherwise
    # disp(['Inner product not available for class ' class(Y)]);
# end

return(res);
}

###------------------------##
ttv<-function(a,v,dims=NULL){
# %TTV Tensor times vector for ktensor.
# a tensor
# %
# %   Y = TTV(X,A,N) computes the product of Kruskal tensor X with a
# %   (column) vector A.  The integer N specifies the dimension in X
# %   along which A is multiplied.  If size(A) = [I,1], then X must have
# %   size(X,N) = I.  Note that ndims(Y) = ndims(X) - 1 because the N-th
# %   dimension is removed.
# %
# %   Y = TTV(X,{A1,A2,...}) computes the product of tensor X with a
# %   sequence of vectors in the cell array.  The products are computed
# %   sequentially along all dimensions (or modes) of X. The cell array
# %   contains ndims(X) vectors.
# %
# %   Y = TTV(X,{A1,A2,...},DIMS) computes the sequence tensor-vector
# %   products along the dimensions specified by DIMS.
# %
# %   See also TENSOR/TTV, KTENSOR, KTENSOR/TTM.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.


# %%%%%%%%%%%%%%%%%%%%%%
# %%% ERROR CHECKING %%%
# %%%%%%%%%%%%%%%%%%%%%%

# % Check the number of arguments
# if (nargin < 2)
    # error('TTV requires at least two arguments.');
# end

	# % Check for 3rd argument
	if (is.null(dims)){
		dims = list();
	}

	# % Check that 2nd argument is cell array. If not, recall with v as a
	# % cell array with one element.
	if (mode(v)!='list'){
		cc = ttv(a,list(v),dims);
		return(cc);
	}

	# % Get sorted dims and index for multiplicands
	tmp = tt_dimscheck(dims,ndims(a),length(v));       
	dims=tmp[[1]];vidx=tmp[[2]];
	# % Check that each multiplicand is the right size.
	# for i = 1:numel(dims)
		# if ~isequal(size(v{vidx(i)}),[size(a,dims(i)) 1])
			# error('Multiplicand is wrong size');
		# end
	# end

	# % Figure out which dimensions will be left when we're done
	
	remdims = 1:ndims(a)
	remdims=remdims[-dims]

	# % Collapse dimensions that are being multiplied out
	newlambda = a$lambda;
	for(i in 1:length(dims)){ 
			newlambda = newlambda * ( Matrix::t(a$u[[dims[i]]]) %*% v[[vidx[i]]] );
	}

	# % Create final result
	if (length(remdims)==0){
		 cc = sum(newlambda);
	}else{
		ul=list()
		for(i in remdims) ul[[i]]=a$u[[i]]
		cc = list(newlambda,ul);#ktensor(newlambda,a.u{remdims});
	}

	return(cc)
}

tt_dimscheck<-function(dims,N,M=NULL){
# %TT_DIMSCHECK Used to preprocess dimensions tensor dimensions.
# %
# %   NEWDIMS = TT_DIMCHECK(DIMS,N) checks that the specified dimensions
# %   are valid for a tensor of order N. If DIMS is empty, then
# %   NEWDIMS=1:N. If DIMS is negative, then NEWDIMS is everything
# %   but the dimensions specified by -DIMS. Finally, NEWDIMS is
# %   returned in sorted order.
# %
# %   [NEWDIMS,IDX] = TT_DIMCHECK(DIMS,N,M) does all of the above but
# %   also returns an index for M muliplicands. 
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.

# % Fix empty case
if (length(dims)==0){
    dims = 1:N;
}

# % Fix "minus" case
if (max(dims) < 0){
    # % Check that every member of dims is in 1:N
    # tf = ismember(-dims,1:N);
    if (!all(-dims %in% 1:N)){#min(tf) == 0
        stop('Invalid dimensions specified');
    }
    dims = (1:N)[dims];
}

# % Check that every member of dims is in 1:N
# tf = ismember(dims,1:N);
if (!all(dims %in% 1:N)){
    stop('Invalid dimensions specified');
}

# % Save the number of dimensions in dims
P = length(dims);

# % Reorder dims from smallest to largest (this matters in particular
# % for the vector multiplicand case, where the order affects the
# % result)
tmp = sort(dims,index.return=TRUE);
sdims=tmp$x; sidx=tmp$ix;

if(is.null(M)) return(sdims);
# if (nargout == 2)
    # % Can't have more multiplicands them dimensions
    if (M > N){
        stop('Cannot have more multiplcands than dimensions');
    }
    
    # % Check that the number of mutliplicands must either be
    # % full-dimensional (i.e., M==N) or equal to the number of specified
    # % dimensions (i.e., M==P).
    if ((M != N) && (M != P)){
        stop('Invalid number of multiplicands');
    }
    
    # % Check sizes to determine how to index multiplicands
    if (P == M){
        # % Case 1: Number of items in dims and number of multiplicands
        # % are equal; therefore, index in order of how sdims was sorted.
        vidx = sidx;   
    }else{
        # % Case 2: Number of multiplicands is equal to the number of
        # % dimensions in the tensor; therefore, index multiplicands by
        # % dimensions specified in dims argument.
        vidx = sdims; # % index multiplicands by (sorted) dimension
    }
	
	return(list(sdims,vidx))
}

spt_ttv<-function(a,v,dims=NULL){
# %TTV Sparse tensor times vector.
# %
# %   Y = TTV(X,V,N) computes the product of a sparse tensor X with a
# %   (column) vector V.  The integer N specifies the dimension in X
# %   along which V is multiplied.  If size(V) = [I,1], then X must have
# %   size(X,N) = I.  Note that ndims(Y) = ndims(X) - 1 because the N-th
# %   dimension is removed. 
# %
# %   Y = TTV(X,U) computes the product of a sparse tensor X with a
# %   sequence of vectors in the cell array U.  The products are
# %   computed sequentially along all dimensions (or modes) of X. The
# %   cell array U contains ndims(X) vectors.
# %
# %   Y = TTV(X,U,DIMS) computes the sequence tensor-vector products
# %   along the dimensions specified by DIMS.
# %
# %   In all cases, the result Y is a sparse tensor if it has 50% or
# %   fewer nonzeros; otherwise ther result is returned as a dense
# %   tensor.
# %
# %   See also SPTENSOR, SPTENSOR/TTM, TENSOR, TENSOR/TTV.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.

# % Check the number of arguments
# if (nargin < 2)
    # error('TTV requires at least two arguments.');
# end

# % Check for 3rd argument
	if (is.null(dims)){
		dims = list();
	}

# % Check that 2nd argument is cell array. If not, recall with v as a
# % cell array with one element.
	if (mode(v)!='list'){
		cc = spt_ttv(a,list(v),dims);
		return(cc);
	}

# % Get sorted dims and index for multiplicands
tmp = tt_dimscheck(dims,ndims(a),length(v));       
	dims=tmp[[1]];vidx=tmp[[2]];
	remdims = 1:ndims(a)
	remdims=remdims[-dims]

# % Check that each multiplicand is the right size.
# for i = 1:numel(dims)
    # if ~isequal(size(v{vidx(i)}),[size(a,dims(i)) 1])
        # error('Multiplicand is wrong size');
    # end
# end

# % Multiply each value by the appropriate elements of the
# % appropriate vector
newvals = a$vals;
subs = a$subs;
for (n in 1:length(dims)){
     idx = subs[,dims[n]];# % extract indices for dimension n
     w = v[[vidx[n]]];        #% extract nth vector
     bigw = w[idx];         #% stretch out the vector
     newvals = newvals * bigw;
}

# % Case 0: If all dimensions were used, then just return the sum
if (length(remdims)==0){
		 cc = sum(newvals);
		 return(cc)
}

# % Otherwise, figure out the subscripts and accumuate the results.
newsubs = a$subs[,remdims];
newsiz = a$size[remdims];
# % Case I: Result is a vector
if (length(remdims) == 1){
    cc = pracma::accumarray(newsubs,newvals,sz=newsiz);
    # if (nnz(cc) <= 0.5 * newsiz){
        # cc = sptensor((1:newsiz)',c,newsiz);
        # cc = list(subs=(1:newsiz),vals=cc,size=newsiz);
    # }else{
        # cc = tensor(c,newsiz);
    # }
 return(cc);
}

# % Case II: Result is a multiway array
# c = sptensor(newsubs, newvals, newsiz);

# % Convert to a dense tensor if more than 50% of the result is nonzero.
# if nnz(c) > 0.5 * prod(c.size)
    # c = tensor(c);
# end

return(NULL)
}


 spt_ones<-function(X){
 return(list(subs=X$subs,vals=rep(1,length(X$vals)),size=X$size,nnz=X$nnz))
 }
 
redistribute<-function(X,pmode){
# %REDISTRIBUTE Distribute lambda values to a specified mode. 
# %
# %   K = REDISTRIBUTE(K,N) absorbs the weights from the lambda vector
# %   into mode N. Set the lambda vector to all ones.
# %
# %   See also KTENSOR, NORMALIZE.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.

    lambda=X$lambda
	u=X$u
	for(r1 in 1:length(lambda)){
		u[[pmode]][,r1] = u[[pmode]][,r1] * lambda[r1];
		lambda[r1] = 1;
	}
	return(list(lambda=lambda,u=u))
}
