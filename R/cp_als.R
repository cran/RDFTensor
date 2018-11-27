#17/9/2018
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

cp_als<-function(X,R, opts=list()){
# Return: [P,Uinit,output] =
# %CP_ALS Compute a CP decomposition of any type of tensor.
# %
# %   P = CP_ALS(X,R) computes an estimate of the best rank-R
# %   CP model of a tensor X using an alternating least-squares
# %   algorithm.  The input X can be a tensor, sptensor, ktensor, or
# %   ttensor. The result P is a ktensor.
# %
# %   P = CP_ALS(X,R,'param',value,...) specifies optional parameters and
# %   values. Valid parameters and their default values are:
# %      'tol' - Tolerance on difference in fit {1.0e-4}
# %      'maxiters' - Maximum number of iterations {50}
# %      'dimorder' - Order to loop through dimensions {1:ndims(A)}
# %      'init' - Initial guess [{'random'}|'nvecs'|cell array]
# %      'printitn' - Print fit every n iterations; 0 for no printing {1}
# %
# %   [P,U0] = CP_ALS(...) also returns the initial guess.
# %
# %   [P,U0,out] = CP_ALS(...) also returns additional output that contains
# %   the input parameters.
# %
# %   Note: The "fit" is defined as 1 - norm(X-full(P))/norm(X) and is
# %   loosely the proportion of the data described by the CP model, i.e., a
# %   fit of 1 is perfect.
# %
# %   NOTE: Updated in various minor ways per work of Phan Anh Huy. See Anh
# %   Huy Phan, Petr Tichavsky, Andrzej Cichocki, On Fast Computation of
# %   Gradients for CANDECOMP/PARAFAC Algorithms, arXiv:1204.1586, 2012.
# %
# %   Examples:
# %   X = sptenrand([5 4 3], 10);
# %   P = cp_als(X,2);
# %   P = cp_als(X,2,'dimorder',[3 2 1]);
# %   P = cp_als(X,2,'dimorder',[3 2 1],'init','nvecs');
# %   U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of P
# %   [P,U0,out] = cp_als(X,2,'dimorder',[3 2 1],'init',U0);
# %   P = cp_als(X,2,out.params); %<-- Same params as previous run
# %
# %   See also KTENSOR, TENSOR, SPTENSOR, TTENSOR.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.

# % This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.
# % http://www.sandia.gov/~tgkolda/TensorToolbox.
# % Copyright (2015) Sandia Corporation. Under the terms of Contract
# % DE-AC04-94AL85000, there is a non-exclusive license for use of this
# % work by or on behalf of the U.S. Government. Export of this data may
# % require a license from the United States Government.
# % The full license terms can be found in the file LICENSE.txt



# %% Extract number of dimensions and norm of X.
N = 3;# 3-way tensors: RDFTensor
normX = normX=sqrt(sum(X$vals^2));;
sizeX=X$size;

# %% Copy from params object
# fitchangetol = params.Results.tol;
fitchangetol = ifelse(is.null(opts[['tol']]),1e-4,opts[['tol']]);
# maxiters = params.Results.maxiters;
maxiters = ifelse(is.null(opts[['maxiters']]),50,opts[['maxiters']]);
# dimorder = params.Results.dimorder;
if(is.null(opts[['dimorder']])){
	dimorder = 1:N;
}else{
	dimorder=opts[['dimorder']];
}

# init = params.Results.init;
init = ifelse(is.null(opts[['init']]),'random',opts[['init']]);
# printitn = params.Results.printitn;
printitn = ifelse(is.null(opts[['printitn']]),1,opts[['printitn']]);

# %% Error checking 

# %% Set up and error checking on initial guess for U.
if (is(init)[1]=="list"){
    Uinit = init;
    if(length(Uinit)!=N){
        stop(sprintf('OPTS.init does not have %d cells',N))
    }    
}else{
    if (init=='random'){
    # % Observe that we don't need to calculate an initial guess for the
        # % first index in dimorder because that will be solved for in the first
        # % inner iteration.    
    Uinit = list(N);#cell(N,1);
        for (n in  dimorder){#calculate for first also
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

if (printitn>0)#verbose
  print('CP_ALS:');


# %% Main Loop: Iterate until convergence

# if (isa(X,'sptensor') || isa(X,'tensor')) && (exist('cpals_core','file') == 3)
 
    # %fprintf('Using C++ code\n');
    # [lambda,U] = cpals_core(X, Uinit, fitchangetol, maxiters, dimorder);
    # P = ktensor(lambda,U);
    
# else
    # UtU = zeros(R,R,N);
    UtU=list(N)
    for(n in 1:N){
        UtU[[n]]=matrix(0,R,R)
    }
    for(n in 1:N){
        if (!is.null(U[[n]])){#~isempty(U{n})
            UtU[[n]] = Matrix::t(U[[n]]) %*% U[[n]];
        }
    }
    
    for(iter in 1:maxiters){
        t0=proc.time()
        fitold = fit;
        
        # % Iterate over all N modes of the tensor
        for (n in dimorder){
            
            # % Calculate Unew = X_(n) * khatrirao(all U except n, 'r').
            if(printitn>0) print( "INFO:Calculate X_(n) * khatrirao(all U except n,)...");
            Unew = spt_mttkrp(X,U,n);
            # print(n)
            # % Compute the matrix of coefficients for linear system
            tmp_dd=(1:N)[-n]
            Y=UtU[[tmp_dd[1]]]*UtU[[tmp_dd[2]]]
            # Y = prod(UtU(:,:,[1:n-1 n+1:N]),3);
            (Unew = Unew %*% solve(Y));
            # if issparse(Unew)
                # Unew = full(Unew);  # % for the case R=1
            # end
# browser()                        
            # % Normalize each vector to prevent singularities in coefmatrix
            if(iter == 1){
                # lambda = sqrt(sum(Unew.^2,1))'; %2-norm
                lambda = sqrt(colSums(Unew^2));# %2-norm
            }else{
                # lambda = max( max(abs(Unew),[],1), 1 )'; %max-norm
                lambda = apply(Unew,2,max); #%max-norm
            }         
            
            # Unew = bsxfun(@rdivide, Unew, lambda');
            tmpUnew = Matrix::t(apply(Unew,1,function(x){x/lambda}))

            U[[n]] = tmpUnew;
            # UtU(:,:,n) = U{n}'*U{n};
            UtU[[n]] = Matrix::t(U[[n]]) %*% U[[n]];
        }
        
        # P = ktensor(lambda,U);
		P=list(lambda=lambda,u=U)

        if (normX == 0){
            # fit = norm(P)^2 - 2 * innerprod(X,P);
            fit = ktensor_norm(P)^2 - 2 * tt_innerprod(P,X);
        }else{
            # normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
            normresidual = sqrt( normX^2 + ktensor_norm(P)^2 - 2 * tt_innerprod(P,X) );
            fit = 1 - (normresidual / normX); #%fraction explained by model
        }
        fitchange = abs(fitold - fit);
        
        # % Check for convergence
        if ((iter > 1) && (fitchange < fitchangetol)){
            flag = 0;
        }else{
            flag = 1;
        }
        t3=proc.time()
        if ((printitn >0 && (iter%%printitn)==0) || ((printitn>0) && (flag==0))){
            # print(sprintf(' Iter %2d: f = %e f-delta = %7.1e\n', iter, fit, fitchange);
            print(sprintf(' Iter %2d: fit = %e fitdelta = %7.1e, itime=%.2f', iter, fit, fitchange,itime=(t3-t0)[3]));
        }
        
        # % Check for convergence
        if (flag == 0)
            break;
                
    }  



# %% Clean up final result
# % Arrange the final tensor so that the columns are normalized.
if (printitn>0) print("arranging results...")
P = arrange(P);
# % Fix the signs
P = fixsigns(P);

if (printitn>0){
    if (normX == 0){
        # fit = norm(P)^2 - 2 * innerprod(X,P);
        fit = ktensor_norm(P)^2 - 2 * tt_innerprod(P,X);
    }else{
        # normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
          normresidual = sqrt( normX^2 + ktensor_norm(P)^2 - 2 * tt_innerprod(P,X) );
        fit = 1 - (normresidual / normX); #%fraction explained by model
    }
  print(sprintf(' Final f = %e ', fit));
}

    return(list(P=P,Uinit=Uinit,iters=iter,fit=fit))
}

fixsigns<-function(K){
# %FIXSIGNS Fix sign ambiguity of a ktensor.
# %
# %   K = FIXSIGNS(K) makes it so that the largest magnitude entries for
# %   each vector in each factor of K are positive, provided that the
# %   sign on *pairs* of vectors in a rank-1 component can be flipped.
# %
# %   K = FIXSIGNS(K,K0) returns a version of K where some of the signs of
# %   the columns of the factor matrices have been flipped to better align
# %   with K0. 
# %
# %   See also KTENSOR and KTENSOR/ARRANGE.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.

# % This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.
# % http://www.sandia.gov/~tgkolda/TensorToolbox.
# % Copyright (2015) Sandia Corporation. Under the terms of Contract
# % DE-AC04-94AL85000, there is a non-exclusive license for use of this
# % work by or on behalf of the U.S. Government. Export of this data may
# % require a license from the United States Government.
# % The full license terms can be found in the file LICENSE.txt


# if nargin == 1
    # K = fixsigns_oneargin(K);
# else 
    # K = fixsigns_twoargin(K, K0);
# end


# %%
# function K = fixsigns_oneargin(K)
    R = length(K$lambda);
    for (r in 1 : R){    
        sgn=numeric(length(K$u))
        for (n in 1:length(K$u)){
            # [val(n),idx(n)] = max(abs(K.u{n}(:,r)));
            idx= which.max(abs(K$u[[n]][,r]));            
            # sgn(n) = sign(K.u{n}(idx(n),r));
            sgn[n] = sign(K$u[[n]][idx,r]);
        }

        negidx = which(sgn == -1);
        nflip = 2 * floor(length(negidx)/2);
        if(nflip>0){
            for( i in 1:nflip){
                n = negidx[i];
                K$u[[n]][,r] =  -1*K$u[[n]][,r];
            }
        }
    }
  return(K)
}