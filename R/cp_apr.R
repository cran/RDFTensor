# 28/5/2018
# Trnaslated from matlab file cp_apr.m  MATLAB Tensor Toolbox.

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

 cp_apr <- function(X, R, opts=list()){
 # Output: [M, Minit, output] =
# %CP_APR Compute nonnegative CP with alternating Poisson regression.
# %
# %   M = CP_APR(X, R) computes an estimate of the best rank-R CP model of a
# %   nonnegative tensor X using an alternating Poisson regression. This is
# %   most appropriate for sparse count data (i.e., nonnegative integer
# %   values) because it uses Kullback-Liebler divergence.  The input X can
# %   be sptensor( list(subs,vals,size,nnz)). The result M is a ktensor.  Input data must be
# %   nonnegative, and the computed ktensor factors are all nonnegative.   
# %
# %   Different algorithm variants are available (selected by the 'alg'
# %   parameter):
# %     'pqnr' - row subproblems by projected quasi-Newton (default)
# %     'pdnr' - row subproblems by projected damped Hessian
# %     'mu'   - multiplicative update (default in version 2.5)
# %
# %   M = CP_APR(X, R, 'param', value, ...) specifies optional parameters and
# %   values. Some parameters work in all situations, others apply only for
# %   a particular choice of algorithm.
# %
# %   Valid parameters and their default values are:
# %      'alg'           - Algorithm ['mu'|'pdnr'|'pqnr'] {'pqnr'}
# %      'stoptol'       - Tolerance on the overall KKT violation {1.0e-4}
# %      'stoptime'      - Maximum number of seconds to run {1e6}
# %      'maxiters'      - Maximum number of iterations {1000}
# %      'init'          - Initial guess [{'random'}|ktensor]
# %      'maxinneriters' - Maximum inner iterations per outer iteration {10}
# %      'epsDivZero'    - Safeguard against divide by zero {1.0e-10}
# %      'printitn'      - Print every n outer iterations; 0 for none {1}
# %      'printinneritn' - Print every n inner iterations {0}
# %
# %   Additional input parameters for algorithm 'mu':
# %      'kappa'         - Offset to fix complementary slackness {100}
# %      'kappatol'      - Tolerance on complementary slackness {1.0e-10}
# %
# %   Additional input parameters for algorithm 'pdnr':
# %      'epsActive'     - Bertsekas tolerance for active set {1.0e-8}
# %      'mu0'           - Initial damping parameter {1.0e-5}
# %      'precompinds'   - Precompute sparse tensor indices {TRUE}
# %      'inexact'       - Compute inexact Newton steps {TRUE}
# %
# %   Additional input parameters for algorithm 'pqnr':
# %      'epsActive'     - Bertsekas tolerance for active set {1.0e-8}
# %      'lbfgsMem'      - Number vector pairs to store for L-BFGS {3}
# %      'precompinds'   - Precompute sparse tensor indices {TRUE}
# %
# %   [M,M0] = CP_APR(...) also returns the initial guess.
# %
# %   [M,M0,out] = CP_APR(...) also returns additional output:
# %      out.kktViolations - maximum KKT violation per iteration
# %      out.nInnerIters   - number of inner iterations per outer iteration
# %      out.obj           - final negative log-likelihood objective
# %      out.ttlTime       - time algorithm took to converge or reach max
# %      out.times         - cumulative time through each outer iteration
# %    If algorithm is 'mu':
# %      out.nViolations   - number of factor matrices needing complementary
# %                          slackness adjustment per iteration
# %    If algorithm is 'pdnr' or 'pqnr':
# %      out.nZeros        - number of zero factor entries per iteration
# %
# %   REFERENCES: 
# %   * E. C. Chi and T. G. Kolda. On Tensors, Sparsity, and Nonnegative
# %     Factorizations, SIAM J. Matrix Analysis,  33(4):1272-1299, Dec. 2012,
# %     http://dx.doi.org/10.1137/110859063  
# %   * S. Hansen, T. Plantenga and T. G. Kolda, Newton-Based Optimization
# %     for Kullback-Leibler Nonnegative Tensor Factorizations, 
# %     Optimization Methods and Software, 2015, 
# %     http://dx.doi.org/10.1080/10556788.2015.1009977
# %
# %   See also CP_ALS, KTENSOR, TENSOR, SPTENSOR.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.


	# %% Set the algorithm choice and initial guess from input or defaults.
	# params = inputParser;
	# params.addParameter('alg', 'pqnr', @(x) (ismember(x,{'mu','pdnr','pqnr'})) );
	opts[['alg']]= ifelse(is.null(opts[['alg']]),'pqnr',opts[['alg']]);
	if(!(opts[['alg']] %in% c('mu','pdnr','pqnr'))) stop(sprintf('invalid alg:%s, possible values are mu,pdnr,pqnr',opts[['alg']]))
	# params.addParameter('init','random', @(x) (isa(x,'ktensor') || ismember(x,{'random'})) );
	if(is.null(opts[['init']])){
			opts[['init']] = 'random';
	}
	opts[['KeepUnmatched']] = TRUE;
	# params.parse(varargin{:});

	alg   = opts[['alg']];
	Minit = opts[['init']];

	# % Extract the number of modes in tensor X.
	N = ndims(X);

	if (R <= 0){
		stop('Number of components(Rank) requested must be positive');
	}

	# %% Check that the data is nonnegative.
	# tmp = which(X < 0.0);
	if (any(X$vals < 0.0)){
		stop('Data tensor must be nonnegative for Poisson-based factorization');
	}

# %% Set up an initial guess for the factor matrices.
	if (length(opts[['init']])!=1 ){#not 'random'
		# % User provided an initial ktensor; validate it.

		if (ncol(Minit$u[[1]]) != R){
			stop('Initial guess does not have the right number of components(columns)');
		}

		for (n in 1:N){
			if (nrow(Minit$u[[n]]) != X$size[n]){
				stop(sprintf('Mode %d of the initial guess is the wrong size',n));
			}
			# if (min(min(Minit.U{n})) < 0.0){
				# stop(sprintf('Initial guess has negative element in mode %d',n));
			# }
		}
		if (any(Minit$lambda < 0.0)){
			stop('Initial guess has a negative ktensor weight');
		}

	}else{
		# % Choose random values for each element in the range (0,1).
		Uinit=list(N);#cell(N,1);
		for (n in  1:N){
				Uinit[[n]] = matrix(runif(X$size[n]*R)+0.1, ncol=R,byrow=TRUE)#rand(size(X,n),R) + 0.1;
		}

		Minit = list(lambda=rep(1,R),u=Uinit);#ktensor(F);
	}


# %% Call a solver based on the choice of algorithm parameter, passing
# %  all the other input parameters.
	outparams=list()
	if (alg=='mu'){
		tmp = tt_cp_apr_mu (X, R, Minit, opts);#[['Unmatched']]
		M=tmp$M; output=tmp$output;	
		 outparams[['alg']] = 'mu';
	}else{
	  if (alg=='pdnr'){
		 tmp = tt_cp_apr_pdnr (X, R, Minit, opts);
		  M=tmp$M; output=tmp$output;	
		 outparams[['alg']] = 'pdnr';
	  }else{
		  if (alg=='pqnr'){
		  tmp = tt_cp_apr_pqnr (X, R, Minit, opts);#[['Unmatched']]
		  M=tmp$M; output=tmp$output;		  
		  outparams[['alg']] = 'pqnr';
		 }
	  }
	 }

	 return(list(M=M,Minit=Minit,output=output,outparams=outparams))
	}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %  Main algorithm PQNR
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  tt_cp_apr_pqnr <- function(X, R, Minit, opts){
# %TT_CP_APR_PQNR Compute nonnegative CP with alternating Poisson regression.
# %
# %   tt_cp_apr_pqnr(X, R, ...) computes an estimate of the best rank-R
# %   CP model of a tensor X using an alternating Poisson regression.
# %   The algorithm solves "row subproblems" in each alternating subproblem,
# %   using a quasi-Newton Hessian approximation.
# %   The function is typically called by cp_apr.
# %
# %   The model is solved by nonlinear optimization, and the code literally
# %   minimizes the negative of log-likelihood.  However, printouts to the
# %   console reverse the sign to show maximization of log-likelihood.
# %
# %   The function call can specify optional parameters and values.
# %   Valid parameters and their default values are:
# %      'stoptol'       - Tolerance on the overall KKT violation {1.0e-4}
# %      'stoptime'      - Maximum number of seconds to run {1e6}
# %      'maxiters'      - Maximum number of iterations {1000}
# %      'maxinneriters' - Maximum inner iterations per outer iteration {10}
# %      'epsDivZero'    - Safeguard against divide by zero {1.0e-10}
# %      'printitn'      - Print every n outer iterations; 0 for no printing {1}
# %      'printinneritn' - Print every n inner iterations {0}
# %      'epsActive'     - Bertsekas tolerance for active set {1.0e-8}
# %      'lbfgsMem'      - Number vector pairs to store for L-BFGS {3}
# %      'precompinds'   - Precompute sparse tensor indices to run faster {TRUE}
# %
# %   Return values are:
# %      M                 - ktensor model with R components
# %      out.fnEvals       - number of row obj fn evaluations per outer iteration
# %      out.kktViolations - maximum KKT violation per iteration
# %      out.nInnerIters   - number of inner iterations per outer iteration
# %      out.nZeros        - number of factor elements equal to zero per iteration
# %      out.obj           - final log-likelihood objective
# %                          (minimization objective is actually -1 times this)
# %      out.ttlTime       - time algorithm took to converge or reach max
# %      out.times         - cumulative time through each outer iteration
# %
# %   REFERENCE: Samantha Hansen, Todd Plantenga, Tamara G. Kolda.
# %   Newton-Based Optimization for Nonnegative Tensor Factorizations,
# %   arXiv:1304.4964 [math.NA], April 2013,
# %   URL: http://arxiv.org/abs/1304.4964. Submitted for publication.
# %
# %   See also CP_APR, KTENSOR, TENSOR, SPTENSOR.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.


# %% Set algorithm parameters from input or by using defaults.
# params = inputParser;
# params.addParamValue('epsActive', 1e-8, @isscalar);
epsActSet= ifelse(is.null(opts[['epsActive']]),1e-8,opts[['epsActive']]);
epsDivZero= ifelse(is.null(opts[['epsDivZero']]),1e-10,opts[['epsDivZero']]);
nSizeLBFGS= ifelse(is.null(opts[['lbfgsMem']]),3,opts[['lbfgsMem']]);
maxInnerIters= ifelse(is.null(opts[['maxinneriters']]),10,opts[['maxinneriters']]);
maxOuterIters= ifelse(is.null(opts[['maxiters']]),200,opts[['maxiters']]);
precomputeSparseIndices= ifelse(is.null(opts[['precompinds']]),TRUE,opts[['precompinds']]);
printInnerItn= ifelse(is.null(opts[['printinneritn']]),0,opts[['printinneritn']]);
printOuterItn= ifelse(is.null(opts[['printitn']]),1,opts[['printitn']]);
stoptime= ifelse(is.null(opts[['stoptime']]),1e6,opts[['stoptime']]);
stoptol= ifelse(is.null(opts[['stoptol']]),1e-4,opts[['stoptol']]);

# % Extract the number of modes in tensor X.
N = ndims(X);

# % If the initial guess has any rows of all zero elements, then modify
# % so the row subproblem is not taking log(0).  Values will be restored to
# % zero later if the unfolded X for the row has no nonzeros.
for (n in 1:N){
  rs = rowSums(Minit$u[[n]]);
  # tmpIx = find(rowsum == 0);
  # if (any(rs==0)){
    Minit$u[[n]][rs==0,1] = 1.0e-8;
}

# % Start with the initial guess, normalized using the vector L1 norm.
M = normalize(Minit,normtype='1');

# % Sparse tensor flag affects how Pi and Phi are computed.
# if isa(X,'sptensor')
    isSparse = TRUE;
# else
    # isSparse = FALSE;
# end

	# % Initialize output arrays.
	fnEvals = rep(0,maxOuterIters)#zeros(maxOuterIters,1);
	kktViolations = rep(-1,maxOuterIters)#-ones(maxOuterIters,1);
	nInnerIters = rep(0,maxOuterIters);
	nzeros = rep(0,maxOuterIters);
	times = rep(0,maxOuterIters);

	if (printOuterItn > 0)
		print(sprintf('CP_PQNR (alternating Poisson regression using quasi-Newton)'));

	dispLineWarn = (printInnerItn > 0);

	# % Start the wall clock timer.
	tic=proc.time();


	if (isSparse && precomputeSparseIndices){
		# % Precompute sparse index sets for all the row subproblems.
		# % Takes more memory but can cut execution time significantly in some cases.
		if (printOuterItn > 0)
			print('  Precomputing sparse index sets...');
		
		sparseIx = list(N);
		for (n in 1:N){
			num_rows = nrow(M$u[[n]]);
			sparseIx[[n]] = list(num_rows,1);
			for (jj in 1:num_rows){
				sparseIx[[n]][[jj]] = which(X$subs[,n] == jj);
			}
		}
		if (printOuterItn > 0)
			print('done');
		
	}


# %% Main Loop: Iterate until convergence or a max threshold is reached.
  lno=0;
  for( iter in 1:maxOuterIters){
    isConverged = TRUE;  
    kktModeViolations = rep(0,N);
    countInnerIters = rep(0,N);

    # % Alternate thru each factor matrix, A_1, A_2, ... , A_N.
    for (n in 1:N){

        # % Shift the weight from lambda to mode n.
        M = redistribute(M,n);

        # % Calculate Khatri-Rhao product of all matrices but the n-th.
        # if (isSparse == FALSE)
            # % Data is not a sparse tensor.
            # Pi = tt_calcpi_prowsubprob(X, isSparse, M, R, n, N, []);
            # X_mat = double(tenmat(X,n));
        # end

        num_rows = nrow(M$u[[n]]);
        isRowNOTconverged = rep(0,num_rows);

        # % Loop over the row subproblems in mode n.
        for (jj in 1:num_rows){

            # % Get data values for row jj of matricized mode n.
            if (isSparse){
                # % Data is a sparse tensor.
                if (precomputeSparseIndices == FALSE){
                    sparse_indices = which(X$subs[,n] == jj);
                }else{
                    sparse_indices = sparseIx[[n]][[jj]];
                }
                if (length(sparse_indices)==0){
                    # % The row jj of matricized tensor X in mode n is empty.
                    M$u[[n]][jj,] = 0;
                    next;
                }
                x_row = X$vals[sparse_indices];

                # % Calculate just the columns of Pi needed for this row.
                Pi = tt_calcpi_prowsubprob(X, isSparse, M,R, n, N, sparse_indices);
            }else{
                stop("dense tensors not implemented.")
                # x_row = X_mat[jj,];
            }

            # % Get current values of the row subproblem variables.
            m_row = M$u[[n]][jj,];

            # % Initialize L-BFGS storage for the row subproblem.
            delm = matrix(0,R, nSizeLBFGS);
            delg = matrix(0,R, nSizeLBFGS);
            rho = rep(0,nSizeLBFGS);
            lbfgsPos = 1;
            m_rowOLD = list();
            gradOLD = list();

            # % Iteratively solve the row subproblem with projected qNewton steps.
            for (i in 1:maxInnerIters){
                # % Calculate the gradient.
                tmp = calc_grad(isSparse, Pi, epsDivZero, x_row, m_row);
				gradM=tmp[[1]]; phi_row=tmp[[2]]
                if (i == 1){
                    # % Original cp_aprPQN_row code (and plb_row) does a gradient
                    # % step to prime the L-BFGS approximation.  However, it means
                    # % a row subproblem that already converged wastes time
                    # % doing a gradient step before checking KKT conditions.
                    # % TODO: fix in a future release.
                    m_rowOLD = m_row;
                    gradOLD = gradM;
                    tmp = tt_linesearch_prowsubprob(-t(gradM), t(gradM), m_rowOLD, 1, 1/2, 10, 1.0e-4,
                                                    isSparse, x_row, Pi, 
                                                    phi_row, dispLineWarn);
					
					m_row=tmp[[1]];f=tmp[[2]]; f_unit=tmp[[3]]; f_new=tmp[[4]]; num_evals=tmp[[5]] ;
                    fnEvals[iter] = fnEvals[iter] + num_evals;
                    tmp = calc_grad(isSparse, Pi, epsDivZero,x_row, m_row);
					gradM=tmp[[1]]; phi_row=tmp[[2]]					
                }

                # % Compute the row subproblem kkt_violation.
                # % Experiments in the original paper used this:
                # %kkt_violation = norm(abs(min(m_row,gradM')),2);
                # % Now we use | KKT |_inf:
				min_v=ifelse(m_row< Matrix::t(gradM),m_row,Matrix::t(gradM))
                kkt_violation = max(abs(min_v));
				
                # % Report largest row subproblem initial violation.
                if ((i == 1) && (kkt_violation > kktModeViolations[n])){
                     kktModeViolations[n] = kkt_violation;
                }

                if (printInnerItn!=0 && (i %% printInnerItn) == 0){
                    print(sprintf('    Mode = %1d, Row = %d, InnerIt = %d , RowKKT = %.2e%s', n, jj, i,kkt_violation,
					              ifelse(i==1,'',sprintf(', RowObj = %.4e', -f_new))));
                }

                # % Check for row subproblem convergence.
                if (kkt_violation < stoptol){
                    break;
                }else{
                    # % Not converged, so m_row will be modified.
                    isRowNOTconverged[jj] = 1;
                }

                # % Update the L-BFGS approximation.
                tmp_delm = m_row - m_rowOLD;
                tmp_delg = gradM - gradOLD;
                tmp_rho = 1 / (tmp_delm %*% tmp_delg);
				# browser()
				lno=lno+1
				print(sprintf('-lno=%d-i=%d--jj:%d----tmp_rho:%f,lbfgsPos:%d,nSizeLBFGS:%d, kkt_violation:%f',lno,i,jj,tmp_rho,lbfgsPos,nSizeLBFGS,
				          kkt_violation));
                if ((tmp_rho > 0.0) && (tmp_rho!=Inf)){
                    delm[,lbfgsPos] = tmp_delm;
                    delg[,lbfgsPos] = tmp_delg;
                    rho[lbfgsPos] = tmp_rho;
                }else{
                    # % Rho is required to be positive; if not, then skip
                    # % the L-BFGS update pair.  The recommended safeguard for
                    # % full BFGS is Powell damping, but not clear how to damp
                    # % in 2-loop L-BFGS.
                    if (dispLineWarn){
                        print(sprintf('WARNING: skipping L-BFGS update, rho would be 1 / %.2e', (tmp_delm %*% tmp_delg)));
                    }
                    # % Roll back lbfgsPos since it will increment later.
                    if (lbfgsPos == 1){
                        if (rho[nSizeLBFGS] > 0){
                            lbfgsPos = nSizeLBFGS;
                        }else{
                            # % Fatal error, should not happen.
                            print('ERROR: L-BFGS first iterate is bad, init is returned.');
                            return(list(M=M,out=NULL));
                        }
                    }else{
                        lbfgsPos = lbfgsPos - 1;
                    }
                }
				
                # % Calculate the search direction.
                search_dir = getSearchDirPqnr(m_row, gradM, epsActSet, 
									delm, delg, rho, lbfgsPos, 
                                              i, dispLineWarn);
                lbfgsPos = (lbfgsPos%% nSizeLBFGS) + 1;

                m_rowOLD = m_row;
                gradOLD = gradM;

                # % Perform a projected linesearch and update variables.
                # % Start from a unit step length, decrease by 1/2, stop with
                # % sufficient decrease of 1.0e-4 or at most 10 steps.
                tmp   = tt_linesearch_prowsubprob(Matrix::t(search_dir), Matrix::t(gradOLD), m_rowOLD, 
                                                1, 1/2, 10, 1.0e-4, 
                                                isSparse, x_row, Pi, 
                                                phi_row, dispLineWarn);
				m_row=tmp[[1]];f=tmp[[2]]; f_unit=tmp[[3]]; f_new=tmp[[4]]; num_evals=tmp[[5]] ;
                fnEvals[iter] = fnEvals[iter] + num_evals;
            }#i
			
            M$u[[n]][jj,] = m_row;
            countInnerIters[n] = countInnerIters[n] + i;

        }

        # % Test if all row subproblems have converged, which means that
        # % no variables in this mode were changed.
        if (any(isRowNOTconverged != 0) )
            isConverged = FALSE;
        

        # % Shift weight from mode n back to lambda.
        M = normalize(M,normtype=1,pmode=n);

        # % Total number of inner iterations for a given outer iteration,
        # % totalled across all modes and all row subproblems in each mode.
        nInnerIters[iter] = nInnerIters[iter] + countInnerIters[n];
    }

    # % Save output items for the outer iteration.
    num_zero = 0;
    for (n in 1:N){
        num_zero = num_zero + sum(M$u[[n]] == 0.0);
    }
    nzeros[iter] = num_zero;
    kktViolations[iter] = max(kktModeViolations); 

    # % Print outer iteration status.
    if ((iter%%printOuterItn) == 0){
        print(sprintf('%4d. Ttl Inner Its: %d, KKT viol = %.2e, obj = %.8e, nz: %d', 
        iter, nInnerIters[iter], kktViolations[iter], tt_loglikelihood(X,M), 
        num_zero));
    }

    times[iter] = (proc.time()-tic)[3];

    # % Check for convergence
    if (isConverged)
        break;
    
    if (times[iter] > stoptime){
        print('Exiting because time limit exceeded');
        break;
    }

}

t_stop = proc.time();

# %% Clean up final result and set output items.
M = normalize(M,'sort',1);
loglike = tt_loglikelihood(X,M);

if (printOuterItn > 0){
    # % For legacy reasons, compute "fit", the fraction explained by the model.
    # % Fit is in the range [0,1], with 1 being the best fit.
    normX = sqrt(sum(X$vals^2));   
    normresidual = sqrt( normX^2 + ktensor_norm(M)^2 - 2 * tt_innerprod(M,X) );
    fit = 1 - (normresidual / normX);

    print(sprintf('==========================================='));
    print(sprintf(' Final log-likelihood = %e ', loglike));
    print(sprintf(' Final least squares fit = %e ', fit));
    print(sprintf(' Final KKT violation = %7.7e', kktViolations[iter]));
    print(sprintf(' Total inner iterations = %d', sum(nInnerIters)));
    print(sprintf(' Total execution time = %.2f secs', (t_stop-tic)[3]));
}

out = list();
# out.params = params.Results;
out[['obj']] = loglike;
out[['kktViolations']] = kktViolations[1:iter];
out[['fnEvals']] = fnEvals[1:iter];
out[['nInnerIters']] = nInnerIters[1:iter];
out[['nZeros']] = nzeros[1:iter];
out[['times']] = times[1:iter];
out[['ttlTime']] = (t_stop-tic)[3];#t_stop;

	return(list(M=M,out=out))
}

# %----------------------------------------------------------------------

 calc_grad<-function(isSparse, Pi, eps_div_zero, x_row, m_row){
# %function grad_row = calc_grad(isSparse, Pi, eps_div_zero, x_row, m_row)
# % Compute the gradient for a PQNR row subproblem.
# %
# %   isSparse     - TRUE if x_row is sparse, FALSE if dense
# %   Pi           - matrix
# %   eps_div_zero - safeguard value to prevent division by zero
# %   x_row        - row vector of data values for the row subproblem
# %   m_row        - vector of variables for the row subproblem
# %
# %   Returns the gradient vector for a row subproblem.

    if (isSparse){
        v = m_row %*% Matrix::t(Pi);
		max_v=ifelse(v > eps_div_zero,v,eps_div_zero)
        w = Matrix::t(x_row) / max_v;
		phi_row = w %*% Pi;
 
    # else
        # v = m_row * Pi';
        # w = x_row ./ max(v, eps_div_zero);
        # phi_row = w * Pi;
    }

    # grad_row = (ones(size(phi_row)) - phi_row)';
    grad_row = t(matrix(1,nrow(phi_row),ncol(phi_row)) - phi_row);

  return(list(grad_row=grad_row, phi_row=phi_row))
}

# %----------------------------------------------------------------------

getSearchDirPqnr<-function (m_row, grad, epsActSet, 
                                 delta_m, delta_g, rho, lbfgs_pos, 
                                 iters, disp_warn){
# % Compute the search direction by projecting with L-BFGS.
# %
# %   m_row     - current variable values
# %   grad      - gradient at m_row
# %   epsActSet - Bertsekas tolerance for active set determination
# %   delta_m   - L-BFGS array of vector variable deltas
# %   delta_g   - L-BFGS array of gradient deltas
# %   lbfgs_pos - pointer into L-BFGS arrays
# %
# %   Returns
# %     d       - search direction based on current L-BFGS and grad
# %
# %   Adapted from MATLAB code of Dongmin Kim and Suvrit Sra written in 2008.
# %   Modified extensively to solve row subproblems and use a better linesearch;
# %   see the reference at the top of this file for details.

    lbfgsSize = ncol(delta_m)#size(delta_m,2);

    # % Determine active and free variables.
    # % If epsActSet is zero, then the following works:
    # %   fixedVars = find((m_row == 0) & (grad' > 0));
    # % For the general case this works but is less clear and assumes m_row > 0:
    # %   fixedVars = find((grad' > 0) & (m_row <= min(epsActSet,grad')));
    projGradStep = (m_row - Matrix::t(grad)) * (m_row - Matrix::t(grad) > 0);
    wk = norm(m_row - projGradStep);
    fixedVars = which((Matrix::t(grad) > 0) & (m_row <= min(epsActSet,wk)));

    d = -grad;
    d[fixedVars] = 0;

    if ((Matrix::t(delta_m[,lbfgs_pos]) %*% delta_g[,lbfgs_pos]) == 0.0){
        # % Cannot proceed with this L-BFGS data; most likely the iteration
        # % has converged, so this is rarely seen.
        if (disp_warn)
            print('WARNING: L-BFGS update is orthogonal, using gradient');
        
        return(d);
    }

    alpha = rep(1,lbfgsSize);
    k = lbfgs_pos;

    # % Perform an L-BFGS two-loop recursion to compute the search direction.

    for (i in 1 : min(iters, lbfgsSize)){
        alpha[k] = rho[k] * Matrix::t(delta_m[, k]) %*% d;
        d = d - alpha[k] * delta_g[, k];#d = d - alpha(k) * delta_g(:, k);
        k = lbfgsSize - ((1 - k)%% lbfgsSize);
    }

    coeff = as.numeric(1 / rho[lbfgs_pos] / (Matrix::t(delta_g[, lbfgs_pos]) %*% delta_g[, lbfgs_pos]));
    d = coeff * d;

    for(i in 1 : min(iters, lbfgsSize)){
        k = (k %% lbfgsSize) + 1;
        b = as.numeric(rho[k] * Matrix::t(delta_g[, k]) %*% d);#rho(k) * delta_g(:, k)' * d;
        d = d + (alpha[k] - b) * delta_m[, k];
    }

    d[fixedVars] = 0;
	return(d);
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %  Shared Internal Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'


 tt_calcpi_prowsubprob<-function(X, isSparse, M, R, n, N, sparse_indices){
# % TT_CALCPI_PROWSUBPROB Compute Pi for a row subproblem.
# %
# %   X              - data sptensor
# %   isSparse       - TRUE if X is sparse, FALSE if dense
# %   M              - current factor matrices
# %   R              - number of columns in each factor matrix
# %   n              - mode
# %   N              - number of modes (equals the number of factor matrices)
# %   sparse_indices - indices of row subproblem nonzero elements
# %
# %   Returns Pi matrix.
# %
# %   Intended for use by CP_PDN and CP_PQN.
# %   Based on calculatePi() in CP_APR, which computes for an entire mode
# %   instead of a single row subproblem.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.#'

    if (isSparse){
        # % X is a sparse tensor.  Compute Pi for the row subproblem specified
        # % by sparse_indices.
        num_row_nnz = length(sparse_indices);

        Pi = matrix(1,num_row_nnz, R);
        # for nn = [1:n-1,n+1:N]
		for(nn in (1:N)[-n]){
            # Pi = M{nn}(X.subs(sparse_indices,nn),:) .* Pi;
            Pi = M$u[[nn]][X$subs[sparse_indices,nn],] * Pi;
        }
    # else
        # % X is a dense tensor.  Compute Pi for all rows in the mode.
        # U = M.U;
        # Pi = khatrirao(U{[1:n-1,n+1:N]},'r');
    }
	return(Pi)
}

# %----------------------------------------------------------------------

# function [m_new, f_old, f_1, f_new, num_evals] ...
   tt_linesearch_prowsubprob<-function(d, grad, m_old, step_len, step_red, 
                                max_steps, suff_decr, isSparse, x_row, Pi, 
                                phi_row, disp_warn){
# % TT_LINESEARCH_PROWSUBPROB Perform a line search on a row subproblem.
# %
# %   d         - search direction
# %   grad      - gradient vector at m_old
# %   m_old     - current variable values
# %   step_len  - initial step length, which is the maximum possible step length
# %   step_red  - step reduction factor (suggest 1/2)
# %   max_steps - maximum number of steps to try (suggest 10)
# %   suff_decr - sufficient decrease for convergence (suggest 1.0e-4)
# %   isSparse  - sparsity flag for computing the objective
# %   x_row     - row subproblem data, for computing the objective
# %   Pi        - Pi matrix, for computing the objective
# %   phi_row   - 1-grad, more accurate if failing over to multiplicative update
# %   disp_warn - TRUE means warning messages are displayed
# %
# %   Returns
# %     m_new     - new (improved) variable values
# %     num_evals - number of times objective was evaluated
# %     f_old     - objective value at m_old
# %     f_1       - objective value at m_old + step_len * d
# %     f_new     - objective value at m_new

    minDescentTol = 1.0e-7;
    smallStepTol = 1.0e-7;

    stepSize = step_len;

    # % Evaluate the current objective value.
    f_old = -1 * tt_loglikelihood_row(isSparse, x_row, m_old, Pi);
    num_evals = 1;
    count = 1;

    while (count <= max_steps){
        # % Compute a new step and project it onto the positive orthant.
        m_new = m_old + (stepSize * d);
        m_new = m_new * (m_new > 0);

        # % Check that it is a descent direction.
        gDotd = sum(grad * (m_new - m_old));
        if ((gDotd > 0) || (sum(m_new) < minDescentTol)){
            # % Don't evaluate the objective if not a descent direction
            # % or if all of the elements of m_new are close to zero.
            f_new = Inf;
            if (count == 1)
               f_1 = f_new;
            

            stepSize = stepSize * step_red;
            count = count + 1;
        }else{
            # % Evaluate objective function at new iterate.
            f_new = -1 * tt_loglikelihood_row(isSparse, x_row, m_new, Pi);
            num_evals = num_evals + 1;
            if (count == 1)
               f_1 = f_new;

            # % Check for sufficient decrease.
            if (f_new <= f_old + suff_decr * gDotd){
                break;
            }else{
                stepSize = stepSize * step_red;
                count = count + 1;
            }
        }
    }

    # % Check if the line search failed.
    if (f_1 == Inf){
        # % Unit step failed; return a value that yields ared = 0.
        f_1 = f_old;
    }
    if (   ((count >= max_steps) && (f_new > f_old)) || (sum(m_new) < smallStepTol) ){

        # % Fall back on a multiplicative update step (scaled steepest descent).
        # % Experiments indicate it works better than a unit step in the direction
        # % of steepest descent, which would be the following:
        # % m_new = m_old - (step_len * grad);     % steepest descent
        # % A simple update formula follows, but suffers from round-off error
        # % when phi_row is tiny:
        # % m_new = m_old - (m_old .* grad);
        # % Use this for best accuracy:
        m_new = m_old * phi_row;  #              % multiplicative update

        # % Project to the constraints and reevaluate the subproblem objective.
        m_new = m_new * (m_new > 0);
        f_new = -1 * tt_loglikelihood_row(isSparse, x_row, m_new, Pi);
        num_evals = num_evals + 1;

        # % Let the caller know the search direction made no progress.
        f_1 = f_old;

        if (disp_warn)
            print('WARNING: line search failed, using multiplicative update step');
        
    }

	return(list(m_new=m_new, f_old=f_old, f_1=f_1, f_new=f_new, num_evals=num_evals))
}

# %----------------------------------------------------------------------

# function f = 
tt_loglikelihood_row<-function(isSparse, x, m, Pi){
# %TT_LOGLIKELIHOOD_ROW Compute log-likelihood of one row subproblem.
# %
# %    The row subproblem for a given mode includes one row of matricized tensor
# %    data (x) and one row of the model (m) in the same matricized mode.
# %    Then
# %       (dense case)
# %          m:  R-length vector 
# %          x:  J-length vector
# %          Pi: R x J matrix
# %       (sparse case)
# %          m:  R-length vector
# %          x:  p-length vector, where p = nnz in row of matricized data tensor
# %          Pi: R x p matrix
# %       F = - (sum_r m_r - sum_j x_j * log (m * Pi_j)
# %           where Pi_j denotes the j^th column of Pi
# %           NOTE: Rows of Pi' must sum to one
# %
# %   isSparse - TRUE if x is sparse, FALSE if dense
# %   x        - vector of data values
# %   m        - vector of model values
# %   Pi       - matrix
# %
# %   Returns the log-likelihood probability f.
# %
# %   Intended for use by CP_PDN and CP_PQN.
# %   Similar to tt_loglikelihood() in CP_APR, which computes log likelihood
# %   for the entire tensor instead of a single row subproblem.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.

    term1 = -sum(m);

    if (isSparse){
        term2 = sum(Matrix::t(x) * log(m %*% Matrix::t(Pi)));
    }else{
        b_pi = m * Matrix::t(Pi);
        term2 = 0;
        for (i in 1:length(x)){
            if (x(i) == 0){
                # % Define zero times log(anything) to be zero.
            }else{
                term2 = term2 + x[i] * log(b_pi[i]);
            }
        }
    }

    f = term1 + term2;
	return(f)
}

# %----------------------------------------------------------------------

# function f = 
tt_loglikelihood<-function(X,M){
# %TT_LOGLIKELIHOOD Compute log-likelihood of data X with model M.
# %
# %   F = TT_LOGLIKELIHOOD(X,M) computes the log-likelihood of model M given
# %   data X, where M is a ktensor and X is a tensor or sptensor.
# %   Specifically, F = - (sum_i m_i - x_i * log_i) where i is a multiindex
# %   across all tensor dimensions.
# %
# %   See also cp_apr, tensor, sptensor, ktensor.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.

N = ndims(X);

# if ~isa(M, 'ktensor')
    # error('M must be a ktensor');
# end

M = normalize(M,1,1);

# if isa(X, 'sptensor')
    xsubs = X$subs;
    A = M$u[[1]][xsubs[,1],];
    for (n in 2:N){
       A = A * M$u[[n]][xsubs[,n],];
    }
    f = sum(X$vals * log(rowSums(A))) - sum(M$u[[1]]);
# else
# %{
# % Old code is probably faster, but returns NaN if X and M are both zero
# % for some element.
    # f = sum(sum(double(tenmat(X,1)) .* log(double(tenmat(M,1))))) - sum(sum(M.U{1}));
# %}
    # % The check for x==0 is also in tt_loglikelihood_row.
    # dX = double(tenmat(X,1));
    # dM = double(tenmat(M,1));
    # f = 0;
    # for i = 1:size(dX,1)
      # for j = 1:size(dX,2)
        # if (dX(i,j) == 0.0)
          # % Define zero times log(anything) to be zero.
        # else
          # f = f + dX(i,j) * log(dM(i,j));
        # end
      # end
    # end
    # f = f - sum(sum(M.U{1}));


return(f)	
}


#################################################################
#################################################################
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %  Main algorithm PDNR
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 tt_cp_apr_pdnr<-function (X, R, Minit, opts){#[M, out] 
# %TT_CP_APR_PDNR Compute nonnegative CP with alternating Poisson regression.
# %
# %   tt_cp_apr_pdnr(X, R, ...) computes an estimate of the best rank-R
# %   CP model of a tensor X using an alternating Poisson regression.
# %   The algorithm solves "row subproblems" in each alternating subproblem,
# %   using a Hessian of size R^2.
# %   The function is typically called by cp_apr.
# %
# %   The model is solved by nonlinear optimization, and the code literally
# %   minimizes the negative of log-likelihood.  However, printouts to the
# %   console reverse the sign to show maximization of log-likelihood.
# %
# %   The function call can specify optional parameters and values.
# %   Valid parameters and their default values are:
# %      'stoptol'       - Tolerance on the overall KKT violation {1.0e-4}
# %      'stoptime'      - Maximum number of seconds to run {1e6}
# %      'maxiters'      - Maximum number of iterations {1000}
# %      'maxinneriters' - Maximum inner iterations per outer iteration {10}
# %      'epsDivZero'    - Safeguard against divide by zero {1.0e-10}
# %      'printitn'      - Print every n outer iterations; 0 for no printing {1}
# %      'printinneritn' - Print every n inner iterations {0}
# %      'epsActive'     - Bertsekas tolerance for active set {1.0e-8}
# %      'mu0'           - Initial damping parameter {1.0e-5}
# %      'precompinds'   - Precompute sparse tensor indices to run faster {TRUE}
# %      'inexact'       - Compute inexact Newton steps {TRUE}
# %
# %   Return values are:
# %      M                 - ktensor model with R components
# %      out.fnEvals       - number of row obj fn evaluations per outer iteration
# %      out.kktViolations - maximum KKT violation per iteration
# %      out.nInnerIters   - number of inner iterations per outer iteration
# %      out.nZeros        - number of factor elements equal to zero per iteration
# %      out.obj           - final log-likelihood objective
# %                          (minimization objective is actually -1 times this)
# %      out.ttlTime       - time algorithm took to converge or reach max
# %      out.times         - cumulative time through each outer iteration
# %
# %   REFERENCE: Samantha Hansen, Todd Plantenga, Tamara G. Kolda.
# %   Newton-Based Optimization for Nonnegative Tensor Factorizations,
# %   arXiv:1304.4964 [math.NA], April 2013,
# %   URL: http://arxiv.org/abs/1304.4964. Submitted for publication.
# %
# %   See also CP_APR, KTENSOR, TENSOR, SPTENSOR.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.


	# %% Set algorithm parameters from input or by using defaults.
	epsActSet= ifelse(is.null(opts[['epsActive']]),1e-8,opts[['epsActive']]);
	epsDivZero= ifelse(is.null(opts[['epsDivZero']]),1e-10,opts[['epsDivZero']]);
	maxInnerIters= ifelse(is.null(opts[['maxinneriters']]),10,opts[['maxinneriters']]);
	maxOuterIters= ifelse(is.null(opts[['maxiters']]),400,opts[['maxiters']]);
	mu0= ifelse(is.null(opts[['mu0']]),1e-5,opts[['mu0']]);
	precomputeSparseIndices= ifelse(is.null(opts[['precompinds']]),TRUE,opts[['precompinds']]);
	inexactNewton= ifelse(is.null(opts[['inexact']]),TRUE,opts[['inexact']]);
	printInnerItn= ifelse(is.null(opts[['printinneritn']]),0,opts[['printinneritn']]);
	printOuterItn= ifelse(is.null(opts[['printitn']]),1,opts[['printitn']]);
	stoptime= ifelse(is.null(opts[['stoptime']]),1e6,opts[['stoptime']]);
	stoptol= ifelse(is.null(opts[['stoptol']]),1e-4,opts[['stoptol']]);

	# % Extract the number of modes in tensor X.
	N = ndims(X);
print(printInnerItn)
	# % If the initial guess has any rows of all zero elements, then modify
	# % so the row subproblem is not taking log(0).  Values will be restored to
	# % zero later if the unfolded X for the row has no nonzeros.
	for (n in 1:N){
		rs = rowSums(Minit$u[[n]]);#%sum(A,2)
		Minit$u[[n]][rs==0,1] = 1.0e-8;
	}

	# % Start with the initial guess, normalized using the vector L1 norm.
	M = normalize(Minit,normtype='1');

	# % Sparse tensor flag affects how Pi and Phi are computed.
	# if isa(X,'sptensor')
		isSparse = TRUE;
	# else
		# isSparse = FALSE;
	# end

	# % Initialize output arrays.
	fnEvals = rep(0,maxOuterIters)#zeros(maxOuterIters,1);
	kktViolations = rep(-1,maxOuterIters)#-ones(maxOuterIters,1);
	nInnerIters = rep(0,maxOuterIters);
	nzeros = rep(0,maxOuterIters);
	times = rep(0,maxOuterIters);

	if (printOuterItn > 0)
		print('CP_PDNR (alternating Poisson regression using damped Newton)');

	dispLineWarn = (printInnerItn > 0);

	# % Start the wall clock timer.
	tic=proc.time();


		if (isSparse && precomputeSparseIndices){
			# % Precompute sparse index sets for all the row subproblems.
			# % Takes more memory but can cut execution time significantly in some cases.
			if (printOuterItn > 0)
				print('  Precomputing sparse index sets...');
			
			sparseIx = list(N);
			for (n in 1:N){
				num_rows = nrow(M$u[[n]]);
				sparseIx[[n]] = list(num_rows,1);
				for (jj in 1:num_rows){
					sparseIx[[n]][[jj]] = which(X$subs[,n] == jj);
				}
			}
			if (printOuterItn > 0)
				print('done');
			
		}


	e_vec = rep(1,R);

	rowsubprobStopTol = stoptol;

	# %% Main Loop: Iterate until convergence or a max threshold is reached.
	  lno=0;
	  for( iter in 1:maxOuterIters){
		isConverged = TRUE;  
		kktModeViolations = rep(0,N);
		countInnerIters = rep(0,N);

		# % Alternate thru each factor matrix, A_1, A_2, ... , A_N.
		for(n in 1:N){

			# % Shift the weight from lambda to mode n.
			M = redistribute(M,n);

			# % Calculate Khatri-Rhao product of all matrices but the n-th.
			# if (isSparse == FALSE)
				# % Data is not a sparse tensor.
				# Pi = tt_calcpi_prowsubprob(X, isSparse, M, R, n, N, []);
				# X_mat = double(tenmat(X,n));
			# end
			num_rows = nrow(M$u[[n]]);
			isRowNOTconverged = rep(0,num_rows);

			# % Loop over the row subproblems in mode n.
			  for (jj in 1:num_rows){
				# % Initialize the damped Hessian parameter for the row subproblem.
				mu = mu0;

				# % Get data values for row jj of matricized mode n.
				if (isSparse){
					# % Data is a sparse tensor.
					if (precomputeSparseIndices == FALSE){
						sparse_indices = which(X$subs[,n] == jj);
					}else{
						sparse_indices = sparseIx[[n]][[jj]];
					}
					if (length(sparse_indices)==0){
						# % The row jj of matricized tensor X in mode n is empty.
						M$u[[n]][jj,] = 0;
						next;
					}
					x_row = X$vals[sparse_indices];

					# % Calculate just the columns of Pi needed for this row.
					Pi = tt_calcpi_prowsubprob(X, isSparse, M,R, n, N, sparse_indices);
				}else{
                    stop("dense tensors not implemented.")
					# x_row = X_mat[jj,];
				}
				
				# % Get current values of the row subproblem variables.
				m_row = M$u[[n]][jj,];

				# % Iteratively solve the row subproblem with projected Newton steps.
				innerIterMaximum = maxInnerIters;
				if (inexactNewton && (iter == 1))
					innerIterMaximum = 2;
				
				for (i in 1:innerIterMaximum){
					# % Calculate the gradient.
					tmp= calc_partials(isSparse, Pi, epsDivZero, x_row, m_row);
					phi_row=tmp[[1]]; ups_row=tmp[[2]] 
					gradM = Matrix::t(e_vec - phi_row);

					# % Compute the row subproblem kkt_violation.
					# % Experiments in the original paper used this:
					# %kkt_violation = norm(abs(min(m_row,gradM')),2);
					# % Now we use | KKT |_inf:
					# kkt_violation = max(abs(min(m_row,gradM')));
					min_v=ifelse(m_row < Matrix::t(gradM),m_row,Matrix::t(gradM))
					kkt_violation = max(abs(min_v));
					
					# % Report largest row subproblem initial violation.
					 if ((i == 1) && (kkt_violation > kktModeViolations[n])){
						 kktModeViolations[n] = kkt_violation;
					}
					
					 if (printInnerItn!=0 && (i %% printInnerItn) == 0){
						print(sprintf('    Mode = %1d, Row = %d, InnerIt = %d , RowKKT = %.2e%s', n, jj, i,kkt_violation,
									  ifelse(i==1,'',sprintf(', RowObj = %.4e', -f_new))));
					}

					# % Check for row subproblem convergence.
					if (kkt_violation < rowsubprobStopTol){
						break;
					}else{
						# % Not converged, so m_row will be modified.
						isRowNOTconverged[jj] = 1;
					}

					# % Calculate the search direction.
					tmp = getSearchDirPdnr(Pi, ups_row, R, gradM, m_row, mu, epsActSet);
					search_dir=tmp[[1]]; predicted_red=tmp[[2]]
					
					# % Perform a projected linesearch and update variables.
					# % Start from a unit step length, decrease by 1/2, stop with
					# % sufficient decrease of 1.0e-4 or at most 10 steps.
					tmp   = tt_linesearch_prowsubprob(Matrix::t(search_dir), Matrix::t(gradM), m_row, 
													1, 1/2, 10, 1.0e-4, 
													isSparse, x_row, Pi, 
													phi_row, dispLineWarn);
					m_rowNEW =tmp[[1]];f_old=tmp[[2]]; f_unit=tmp[[3]]; f_new=tmp[[4]]; num_evals=tmp[[5]] ;
					fnEvals[iter] = fnEvals[iter] + num_evals;

					m_row = m_rowNEW;

					# % Update damping parameter mu based on the unit step length,
					# % which is returned in f_unit.
					
					actual_red = f_old - f_unit;
					rho = actual_red / (-predicted_red);
					# if(lno==162900) browser()
					
					if (predicted_red == 0){
						mu = 10 * mu;
					}else{
						if (!is.nan(rho) && rho < 1/4){
							mu = (7/2) * mu;
						}else{
						 if (!is.nan(rho) && rho > 3/4){
							mu = (2/7) * mu;
						}}}
					lno=lno+1
					if(printInnerItn>2)
						print(sprintf('-lno=%d-i=%d--jj:%d----m_row:%f,---rho:%f,mu:%f, kkt_vn:%f, f_new:%f',lno,i,jj,m_row[1],rho,mu,
									kkt_violation,f_new));
					# if((lno%%20)==0) browser();
				}
				M$u[[n]][jj,] = m_row;
				countInnerIters[n] = countInnerIters[n] + i;
			}
			
			# % Test if all row subproblems have converged, which means that
			# % no variables in this mode were changed.
			 if (any(isRowNOTconverged != 0) )
				isConverged = FALSE;

			# % Shift weight from mode n back to lambda.
			M = normalize(M,normtype=1,pmode=n);

			# % Total number of inner iterations for a given outer iteration,
			# % totalled across all modes and all row subproblems in each mode.
			nInnerIters[iter] = nInnerIters[iter] + countInnerIters[n];
		}

		# % Save output items for the outer iteration.
		num_zero = 0;
		for (n in 1:N){
			num_zero = num_zero + sum(M$u[[n]] == 0.0);
		}
		nzeros[iter] = num_zero;
		kktViolations[iter] = max(kktModeViolations); 
		if (inexactNewton)
			rowsubprobStopTol = max(stoptol, kktViolations[iter] / 100.0);
		
		times[iter] = (proc.time()-tic)[3];
		
		# % Print outer iteration status.
		if ((iter%%printOuterItn) == 0){
			print(sprintf('%4d. Ttl Inner Its: %d, KKT viol = %.2e, obj = %.8e, nz: %d, time:%.2f', 
			iter, nInnerIters[iter], kktViolations[iter], tt_loglikelihood(X,M), 
			num_zero, times[iter]));
		}


		# % Check for convergence
		if (isConverged && (inexactNewton == FALSE))
			break;


		if (isConverged && (inexactNewton == TRUE) && (rowsubprobStopTol <= stoptol))
			break;
		
		if (times[iter] > stoptime){
			print('Exiting because time limit exceeded');
			break;
		}

	}

	t_stop = proc.time();

	# %% Clean up final result and set output items.
	M = normalize(M,'sort',1);
	loglike = tt_loglikelihood(X,M);

	if (printOuterItn > 0){
		# % For legacy reasons, compute "fit", the fraction explained by the model.
		# % Fit is in the range [0,1], with 1 being the best fit.
		normX = sqrt(sum(X$vals^2));   
		normresidual = sqrt( normX^2 + ktensor_norm(M)^2 - 2 * tt_innerprod(M,X) );
		fit = 1 - (normresidual / normX);

		print(sprintf('==========================================='));
		print(sprintf(' Final log-likelihood = %e ', loglike));
		print(sprintf(' Final least squares fit = %e ', fit));
		print(sprintf(' Final KKT violation = %7.7e', kktViolations[iter]));
		print(sprintf(' Total inner iterations = %d', sum(nInnerIters)));
		print(sprintf(' Total execution time = %.2f secs', (t_stop-tic)[3]));
	}


	out = list();
	# out.params = params.Results;
	out[['obj']] = loglike;
	out[['kktViolations']] = kktViolations[1:iter];
	out[['fnEvals']] = fnEvals[1:iter];
	out[['nInnerIters']] = nInnerIters[1:iter];
	out[['nZeros']] = nzeros[1:iter];
	out[['times']] = times[1:iter];
	out[['ttlTime']] = (t_stop-tic)[3];#t_stop;

	return(list(M=M,out=out))
}


# %----------------------------------------------------------------------

  calc_partials<-function(isSparse, Pi, eps_div_zero, x_row, m_row){#[phi_row, ups_row] ...
# % Compute derivative quantities for a PDNR row subproblem.
# %
# %   isSparse     - TRUE if x_row is sparse, FALSE if dense
# %   Pi           - matrix
# %   eps_div_zero - safeguard value to prevent division by zero
# %   x_row        - row vector of data values for the row subproblem
# %   m_row        - vector of variables for the row subproblem
# %
# %   Returns two vectors for a row subproblem:
# %     phi_row - gradient of row subproblem, except for a constant
# %     ups_row - intermediate quantity (upsilon) used for second derivatives

    if (isSparse){
        v = m_row %*% Matrix::t(Pi);
		max_v=ifelse(v > eps_div_zero,v,eps_div_zero)
        w = Matrix::t(x_row) / max_v;
        phi_row = w %*% Pi;
        u = v ^ 2;
		max_v=ifelse(u > eps_div_zero,u,eps_div_zero)
        ups_row = Matrix::t(x_row) / max_v;
 
    }# else
        # v = m_row * Pi';
        # w = x_row ./ max(v, eps_div_zero);
        # phi_row = w * Pi;
        # u = v .^ 2;
        # ups_row = x_row ./ max(u, eps_div_zero);
    # end

	return(list(phi_row=phi_row, ups_row=ups_row))
}

# %----------------------------------------------------------------------

 getHessian<-function(upsilon, Pi, free_indices){
# % Return the Hessian for one PDNR row subproblem of M{n}, for just the rows and
# % columns corresponding to the free variables.
    
    num_free = length(free_indices);
    H = matrix(0,num_free,num_free);
    for (i in 1:num_free){
        for (j in i:num_free){
            cc = free_indices[i];
            d = free_indices[j];
            # val = sum(upsilon' .* Pi(:,c) .* Pi(:,d));
            val = colSums(Matrix::t(upsilon) * Pi[,cc] * Pi[,d]);
            H[i,j] = val;
            H[j,i] = val;
        }
    }

	return(H)
}

# %----------------------------------------------------------------------

 getSearchDirPdnr<-function (Pi, ups_row, R, gradM, m_row, mu, epsActSet){
# % Compute the search direction for PDNR using a two-metric projection
# % with damped Hessian.
# %
# %   Pi        - matrix
# %   ups_row   - intermediate quantity (upsilon) used for second derivatives
# %   R         - number of variables for the row subproblem
# %   gradM     - gradient vector for the row subproblem
# %   m_row     - vector of variables for the row subproblem
# %   mu        - damping parameter
# %   epsActSet - Bertsekas tolerance for active set determination
# %
# %   Returns:
# %     search_dir - search direction vector
# %     pred_red   - predicted reduction in quadratic model

    search_dir = rep(0,R);
    # projGradStep = (m_row - gradM') .* (m_row - gradM' > 0);
	projGradStep = (m_row - Matrix::t(gradM)) * (m_row - Matrix::t(gradM) > 0);
    wk = norm(m_row - projGradStep);

    # % Determine active and free variables.
    num_free = 0;
    free_indices_tmp = rep(0,R);
    for (r in 1:R){
        if ((m_row[r] <= min(epsActSet,wk)) && (gradM[r] > 0) ){
            # % Variable is not free (belongs to set A or G).
            if (m_row[r] != 0){
                # % Variable moves according to the gradient (set G).
                search_dir[r] = -gradM[r];
            }
        }else{
            # % Variable is free (set F).
            num_free = num_free + 1;
            free_indices_tmp[num_free] = r;
        }
    }
    free_indices = free_indices_tmp[1:num_free];

    # % Compute the Hessian for free variables.
    Hessian_free = getHessian(ups_row, Pi, free_indices);
    grad_free = -gradM[free_indices];

    # % Compute the damped Newton search direction over free variables.
    # search_dir(free_indices)   = linsolve(Hessian_free + (mu * eye(num_free)), grad_free); 
    search_dir[free_indices]   = solve(Hessian_free + (mu * diag(num_free)), grad_free); 
		
    # % If the Hessian is too ill-conditioned, use gradient descent.
    # msgid = names(last.warning)#lastwarn('MATLAB:noWarning'); 
    # if (msgid == 'MATLAB:nearlySingularMatrix'){## trace!!!
        # print('WARNING: damped Hessian is nearly singular');
        # search_dir = -gradM;
    # }

    # % Calculate expected reduction in the quadratic model of the objective.
    q1 = Matrix::t(search_dir[free_indices]) %*% (Hessian_free + (mu * diag(num_free))) %*% search_dir[free_indices];
    pred_red = (Matrix::t(search_dir[free_indices]) %*% gradM[free_indices]) + (0.5 * q1);
    if (pred_red > 0){
        print('ERROR: expected decrease is positive');
        search_dir = -gradM;
    }

	return(list(search_dir=search_dir, pred_red=pred_red))
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %  Main algorithm MU
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tt_cp_apr_mu<-function(X, R, Minit, opts){ #[M, output] 
# %TT_CP_APR_MU Compute nonnegative CP with alternating Poisson regression.
# %
# %   tt_cp_apr_mu(X, R, ...) computes an estimate of the best rank-R
# %   CP model of a tensor X using an alternating Poisson regression.
# %   The algorithm solves each alternating subproblem using multiplicative
# %   updates with adjustments for values near zero.
# %   The function is typically called by cp_apr.
# %
# %   The function call can specify optional parameters and values.
# %   Valid parameters and their default values are:
# %      'stoptol'       - Tolerance on the overall KKT violation {1.0e-4}
# %      'stoptime'      - Maximum number of seconds to run {1e6}
# %      'maxiters'      - Maximum number of iterations {1000}
# %      'maxinneriters' - Maximum inner iterations per outer iteration {10}
# %      'epsDivZero'    - Safeguard against divide by zero {1.0e-10}
# %      'printitn'      - Print every n outer iterations; 0 for no printing {1}
# %      'printinneritn' - Print every n inner iterations {0}
# %      'kappatol'      - Tolerance on complementary slackness {1.0e-10}
# %      'kappa1'         - Offset to fix complementary slackness {100}
# %
# %   Return values are:
# %      M                 - ktensor model with R components
# %      out.kktViolations - maximum KKT violation per iteration
# %      out.nInnerIters   - number of inner iterations per outer iteration
# %      out.nViolations   - number of factor matrices needing complementary
# %                          slackness adjustment per iteration
# %      out.obj           - final log-likelihood objective
# %      out.ttlTime       - time algorithm took to converge or reach max
# %      out.times         - cumulative time through each outer iteration
# %
# %   REFERENCE: E. C. Chi and T. G. Kolda. On Tensors, Sparsity, and
# %   Nonnegative Factorizations, arXiv:1112.2414 [math.NA], December 2011,
# %   URL: http://arxiv.org/abs/1112.2414. Submitted for publication.
# %
# %   See also CP_APR, KTENSOR, TENSOR, SPTENSOR.
# %
# %MATLAB Tensor Toolbox.
# %Copyright 2015, Sandia Corporation.


# %% Set algorithm parameters from input or by using defaults.

	epsilon= ifelse(is.null(opts[['epsDivZero']]),1e-10,opts[['epsDivZero']]);
	kappa1= ifelse(is.null(opts[['kappa']]),1e-2,opts[['kappa']]);
	kappaTol= ifelse(is.null(opts[['kappatol']]),1e-10,opts[['kappatol']]);
	maxInnerIters= ifelse(is.null(opts[['maxinneriters']]),10,opts[['maxinneriters']]);
	maxOuterIters= ifelse(is.null(opts[['maxiters']]),400,opts[['maxiters']]);
	printInnerItn= ifelse(is.null(opts[['printinneritn']]),0,opts[['printinneritn']]);
	printOuterItn= ifelse(is.null(opts[['printitn']]),1,opts[['printitn']]);
	stoptime= ifelse(is.null(opts[['stoptime']]),1e6,opts[['stoptime']]);
	tol= ifelse(is.null(opts[['stoptol']]),1e-4,opts[['stoptol']]);
	kktViolations=rep(-1,maxOuterIters)
	times = rep(0,maxOuterIters);
	nInnerIters = rep(0,maxOuterIters);
	

# %% Extract dimensions of X and number of dimensions of X.
N = ndims(X);

# % Set up and error checking on initial guess for U.
# if isa(Minit,'ktensor')
    # if ndims(Minit) ~= N
        # error('Initial guess does not have the right number of dimensions');
    # end
    
    # if ncomponents(Minit) ~= R
        # error('Initial guess does not have the right number of components');
    # end
    
    # for n = 1:N
        # if size(Minit,n) ~= size(X,n)
            # error('Dimension %d of the initial guess is the wrong size',n);
        # end
    # end
# elseif strcmp(Minit,'random')
    # F = cell(N,1);
    # for n = 1:N
        # F{n} = rand(size(X,n),R);
    # end
    # Minit = ktensor(F);
# else
    # error('The selected initialization method is not supported');
# end


# %% Set up for iterations - initializing M and Phi.
	# M = normalize(Minit,[],1);
    
	M = normalize(Minit,normtype='1');
	Phi = list(0,0,0);
	kktModeViolations = rep(0,N)#zeros(N,1);

	if (printOuterItn > 0)
	  print('CP_APR: mu');

	nViolations = rep(0,maxOuterIters)#zeros(maxOuterIters,1);

# % Start the wall clock timer.
    tic=proc.time();

# % PDN-R and PQN-R benefit from precomputing sparse indices of X for each
# % mode subproblem.  However, MU execution time barely changes, so the
# % precompute option is not offered.

lno=0;
# %% Main Loop: Iterate until convergence.
for (iter in 1:maxOuterIters){
    
    isConverged = TRUE;   
    for(n in 1:N){

        # % Make adjustments to entries of M{n} that are violating
        # % complementary slackness conditions.
        if (iter > 1){
            V = (Phi[[n]] > 1) & (M$u[[n]] < kappaTol);
            if (any(V)){
                nViolations[iter] = nViolations[iter] + 1;
                M$u[[n]][V>0] = M$u[[n]][V>0] + kappa1;
            }
        }       

        # % Shift the weight from lambda to mode n
		t0=proc.time()
		M = redistribute(X=M,pmode=n);
        t1=proc.time()
        # % Calculate product of all matrices but the n-th
        # % (Sparse case only calculates entries corresponding to nonzeros in X)
		Pi = calculatePi(X, M, R, n, N);
        t2=proc.time()
        # % Do the multiplicative updates
        for( i in 1:maxInnerIters){

            # % Count the inner iterations
            nInnerIters[iter] = nInnerIters[iter] + 1;
                                  
            # % Calculate matrix for multiplicative update
        t3=proc.time()
            Phi[[n]] = calculatePhi(X, M, R, n, Pi, epsilon);
        t4=proc.time()
            
            # % Check for convergence
			min_v=as.vector(M$u[[n]])
			min_v=ifelse(min_v< (1-Phi[[n]]),min_v,(1-Phi[[n]]))
			kktModeViolations[n] = max(abs(min_v));
			if (kktModeViolations[n] < tol){
                break;
            }else{
                isConverged = FALSE;
            }                   
            
            # % Do the multiplicative update
            M$u[[n]] = M$u[[n]] * Phi[[n]];
            # print(sprintf(' Pi t:%.1f, Phi t:%.1f',(t2-t1)[3],(t4-t3)[3]))
            # % Print status
			lno=lno+1
            if (printInnerItn!=0 && (i %% printInnerItn) == 0){
                 print(sprintf('  lno=%d  Mode = %1d, Inner Iter = %2d, KKT violation = %.6e',lno, n, i, kktModeViolations[n]));
            }
			# if(lno==40)browser()
        }
        
        # % Shift weight from mode n back to lambda
        # M = normalize(M,[],1,n);
        t5=proc.time()
		M = normalize(M,normtype=1,pmode=n);
        t6=proc.time()
    }

    kktViolations[iter] = max(kktModeViolations);    

    if ((iter %% printOuterItn)==0){
        print(sprintf(' Iter %4d: Inner Its = %2d KKT violation = %.6e, nViolations = %2d', 
        iter, nInnerIters[iter], kktViolations[iter], nViolations[iter]));            
    }
    times[iter] = (proc.time()-tic)[3];
    
    # % Check for convergence
    if (isConverged){
        if (printOuterItn>0)
            print('Exiting because all subproblems reached KKT tol.');
        
        break;
    }    
    if (times[iter] > stoptime){
        if (printOuterItn>0)
            print('Exiting because time limit exceeded.');
        
        break;
    }
}
# t_stop = toc;
t_stop = proc.time();

# %% Clean up final result
M = normalize(M,'sort',1);

obj = tt_loglikelihood(X,M);
if( printOuterItn>0){
    normX = sqrt(sum(X$vals^2));   
    normresidual = sqrt( normX^2 + ktensor_norm(M)^2 - 2 * tt_innerprod(M,X) );
    fit = 1 - (normresidual / normX);#%fraction explained by model

    print(sprintf('==========================================='));
    print(sprintf(' Final log-likelihood = %e ', obj));
    print(sprintf(' Final least squares fit = %e ', fit));
    print(sprintf(' Final KKT violation = %7.7e', kktViolations[iter]));
    print(sprintf(' Total inner iterations = %d', sum(nInnerIters)));
    print(sprintf(' Total execution time = %.2f secs', (t_stop-tic)[3]));
}

	out = list();
	out[['obj']] = obj;
	out[['kktViolations']] = kktViolations[1:iter];
	out[['nInnerIters']] = nInnerIters[1:iter];
	out[['nViolations']] = nViolations[1:iter]
	out[['nTotalIters']] = sum(nInnerIters);
	out[['times']] = times[1:iter];
	out[['ttlTime']] = (t_stop-tic)[3];#t_stop;

	return(list(M=M,out=out))

}

calculatePi<-function(X, M, R, n, N){

	# if (isa(X,'sptensor'))
		# Pi = matrix(1,X$nnz, R);
		Pi = matrix(1,nrow(X$subs), R);
		# for nn = [1:n-1,n+1:N]
		for(nn in (1:N)[-n]){
			# Pi = M{nn}(X.subs(:,nn),:) .* Pi;
			Pi = M$u[[nn]][X$subs[,nn],] * Pi
		}
	# else
		# U = M.U;
		# Pi = khatrirao(U{[1:n-1,n+1:N]},'r');
	# end
return(Pi)
}

calculatePhi<-function(X, M, R, n, Pi, epsilon){

# if (isa(X,'sptensor'))
t0=proc.time()
    Phi = matrix(-1,X$size[n],R);
    xsubs = X$subs[,n];
    v = rowSums(M$u[[n]][xsubs,] *Pi);
	max_v=ifelse(v>epsilon,v,epsilon)
    wvals = X$vals / max_v#max(v, epsilon);
    t1=proc.time()
	for (r in 1:R){
        # Yr = accumarray(xsubs, wvals * Pi[,r], [size(X,n) 1]);
        # Yr = pracma::accumarray(xsubs,wvals * Pi[,r],sz=X$size[n])#accumarray(, , [size(X,n) 1]);
		Yr=aggregate(wvals * Pi[,r],list(xsubs),sum)# 10 times faster
        Phi[,r] = Yr[,2];
    }
	t2=proc.time()
# else
    # Xn = double(tenmat(X,n));
    # V = M.U{n}*Pi';
    # W = Xn ./ max(V, epsilon);
    # Y = W * Pi;
    # Phi = Y;
# end
# print(sprintf(' Phi tb4:%.1f,  tloop:%.1f',(t1-t0)[3],(t2-t1)[3]))
  return(Phi)
}

# %----------------------------------------------------------------------

# function y = vectorizeForMu(x)
# y = x(:);
# end
