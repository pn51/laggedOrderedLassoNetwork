#' @title .timeLagLassoNetworkLaggedDataCBG
#' @description .timeLagLassoNetworkLaggedDataCBG
#' @keywords internal
#' @param xData lagged expression matrix (n x (p * maxLag))
#' @param yData output expression or change-in-expression matrix (n x p)
#' @param maxLag maximum predictor lag
#' @param lambda a scalar or p x p matrix of penalization parameters. If a scalar, all coefficients are subject to the same penalization. Otherwise,\code{lambda[j,i]} is the penalization on variable \code{j} in the model for gene \code{i} 
#' @param beta_pos list of p optional vectors containing positive parts of the model coefficients (\code{NULL} or vector of length p * maxLag)
#' @param beta_neg list of p optional vectors containing negative parts of the model coefficients (\code{NULL} or vector of length p * maxLag)
#' @param strongly.ordered  if \code{TRUE} (\code{FALSE}), use the strongly (weakly) ordered lasso
#' @param method underlying ordered lasso optimization method ('Solve.QP' or 'GG')
#' @param self if \code{TRUE}, include loops in time-lagged regression models
#' @param maxiter maximum number of  time lagged ordered lasso iterations
#' @param inneriter maximum number ordered lasso iterations
#' @param iter.gg maximum number of generalized gradient iterations
#' @param cores number of parallel cores
#' @import parallel
#' @return lists, one for each gene, containing a list of ordered lasso coefficients corresponding to the positive part of the weakly ordered solution (\code{bp}), the negative part of the weakly ordered solution (\code{bn}), the weakly ordered solution (\code{beta}), the weakly ordered intercept (\code{b0}), the strongly ordered intercept (\code{b0.ordered}), and the strongly ordered solution (\code{beta.ordered}). 

.timeLagLassoNetworkLaggedDataCBG <- function(x, y, maxLag=2, lambda=1,
	beta_pos=NULL, beta_neg=NULL,
	strongly.ordered=FALSE,
	method='Solve.QP', self=TRUE,
	maxiter=500, inneriter=100, iter.gg=100,
	cores=1){

	#same as timeLagLassoNetworkLaggedData, but returns a list of the output of timeLagLassoLaggedData
	#useful for looping through a bunch of close lambda values and computing an AUC

	# assume that x and y have been obtained using functions from dataLagging.R

	#if one value of lambda is given, use it for all genes
	if(length(lambdas)==1){
		lambdas <- rep(lambdas, ncol(y))
	}

	if(!self){
		# no loops in the regression model
		coeffsByGene <- mclapply(seq(ncol(y)), mc.cores=cores, FUN=function(ii){
			# iterate through the genes, learn a model for each
			xInd <- ((ii-1) * maxLag + 1) : (ii * maxLag)
			.timeLagLassoLaggedData(x[,-xInd], y[,ii], maxLag, lambda[-ii, ii],
				beta_pos=beta_pos[[ii]], beta_neg=beta_neg[[ii]],
				method=method, strongly.ordered=strongly.ordered,
				maxiter=maxiter, inneriter=inneriter, iter.gg=iter.gg)
		})
	} else{
		#include loops in the regression models
		coeffsByGene <- mclapply(seq(ncol(y)), mc.cores=cores, FUN=function(ii){
			.timeLagLassoLaggedData(x, y[,ii], maxLag, lambda[, ii], 
				beta_pos=beta_pos[[ii]], beta_neg=beta_neg[[ii]],
				method=method, strongly.ordered=strongly.ordered,
				maxiter=maxiter, inneriter=inneriter, iter.gg=iter.gg)
		})
	}

	coeffsByGene # useful since it contains b0, beta_pos, beta_neg; makes looping through lambda and computing an AUC a little faster
}
