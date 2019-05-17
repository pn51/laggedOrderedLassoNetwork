
##largely adapted/modified from Suo/Tibshirani code
#https://cran.r-project.org/web/packages/orderedLasso/index.html
#version 1.7

# some additions/modifications:
# comments for major changes
# simplified/reformated code
# fixed/shortened some code
# removed redundant code/variables/parameters
# removed unused variables
# removed functions not used for network reconstruction
# some changes for speed, where possible
# modifications for multiple penalization parameters

#' @title .ObjMinFun
#' @description .ObjMinFun
#' @keywords internal
#' @param x input matrix (n x p)
#' @param y output vector (length n)
#' @param bp vector containing positive parts of the model coefficients
#' @param bn vector containing negative parts of the model coefficients
#' @param lambda vector of penalization parameters (length 1 or p)
#' @return lasso objective value
.ObjMinFun <- function(x, y, bp, bn, lambda){
 	# objective function

 	yhat <- x %*% (bp-bn)
 	# modified for separate lambda for each input variable:
  	if(length(lambda) == 1){ #assume length p otherwise
 		rep(lambda, length(bp))
 	}
 	sum((y - yhat)^2) / 2 + drop(crossprod(lambda, bp + bn))
}

#' @title .ObjMinSignFun
#' @description .ObjMinSignFun
#' @keywords internal
#' @param x input matrix (n x p)
#' @param y output vector (length n)
#' @param b0.ordered strongly ordered intercept
#' @param beta.ordered strongly ordered coefficients
#' @param lambda vector of penalization parameters (length 1 or p)
#' @param sign sign vector (usually for beta.ordered)
#' @return lasso objective value
.ObjMinSignFun <- function(x, y, b0.ordered, beta.ordered, lambda, sign){

 	yhat <- x %*% beta.ordered + b0.ordered
 	# modified for separate lambda for each input variable:
 	if(length(lambda) == 1){ #assume length p otherwise
 		rep(lambda, length(beta.ordered))
 	}
 	sum((y - yhat)^2) / 2 + drop(crossprod(lambda, sign * beta.ordered))
}

#' @title .prox
#' @description .prox
#' @keywords internal
#' @param y input vector (length p)
#' @param lambda vector of penalization parameters (length p)
#' @importFrom Iso pava
#' @return pool adjacent violators algorithm fitted values
.prox <- function(y, lambda){
	aa <- pava(y-lambda, decreasing=TRUE)
	aa * (aa>0)
}

#' @title .fastgg
#' @description .fastgg
#' @keywords internal
#' @param x input matrix (n x p)
#' @param y output vector (length n)
#' @param beta vector of current model coefficients
#' @param lambda vector of penalization parameters (length 1 or p)
#' @param t starting step
#' @param gam step length 
#' @param epsilon convergence error tolerance
#' @param iter.gg maximum number of generalized gradient iterations
#' @return generalized gradient-based ordered lasso update
.fastgg <- function(x, y, beta, lambda, t=1, gam=0.8, epsilon=1e-6, inneriter.gg=100){
	# cleaned up original code a bit
	# removed sig and replaced beta_pos/beta_neg with just beta
	b_old_old <- b_old <- b <- beta
	ii <- 0
	while(ii < inneriter.gg){
		ii <- ii + 1
		alpha <- b_old + (ii - 2) / (ii + 1) * (b_old - b_old_old)
		leastsquareconst <- drop(y - x %*% alpha)
		const <- -drop(crossprod(x, leastsquareconst))
		constobjalpha <- drop(crossprod(leastsquareconst) / 2)
		const3 <- drop(crossprod(const, alpha))

		g0 <- 1
		g1 <- 0

		while(g0 > g1 && t > epsilon){
		  b <- .prox(alpha - t * const, t * lambda)
		  residualb <- drop(y - x %*% b)
		  g0 <- drop(crossprod(residualb)) / 2
		  g1 <- drop(constobjalpha + crossprod(const, b) - const3  + crossprod(b - alpha) / (2 * t))
		  t <- gam * t
		}
		if(sum(abs(b - b_old)) < epsilon){
			break
		}
		b_old_old <- b_old
		b_old <- b
	}

	# if(ii == inneriter.gg){
	# 	cat('generalized gradient failed to converge\n')
	# }

	b
}

#' @title .ordLasSignPos
#' @description .ordLasSignPos
#' @keywords internal
#' @param x input matrix (n x p)
#' @param y output vector (length n)
#' @param lambda vector of penalization parameters (length 1 or p)
#' @param signvec sign vector (length p)
#' @importFrom quadprog solve.QP
#' @return strongly ordered lasso solution
.ordLasSignPos <- function(x, y, lambda, signvec){
	# cleaned up original code
	p <- ncol(x)

	A <- diag(signvec) #note: possible problems with signvec==0; see comment in .timeLagLassoEstOrdered
	A[col(A) == (row(A) + 1)] <- -signvec[-1]
	A <- A[-p, , drop=F]
	Amat <- t(rbind(A, diag(signvec)))

	emin <- -min(eigen(crossprod(x))$val)
	emin <- max(emin, 1e-4) #max(.001,emin)
	dmat <- crossprod(x) + diag(rep(emin, p))
	dvec <- drop(crossprod(x, y)) - lambda * signvec
	bvec <- rep(0, ncol(Amat))

	solve.QP(dmat, dvec, Amat, bvec)$sol
}

#' @title .ordLas2
#' @description .ordLas2
#' @keywords internal
#' @param x input matrix (n x p)
#' @param y output variable vector
#' @param lambda vector of penalization parameters (length 1 or p)
#' @importFrom quadprog solve.QP
#' @return ordered lasso solution
.ordLas2 <- function(x, y, lambda){
	# cleaned up original code
	# also modified for separate lambda for each input variable
 	if(length(lambda) == 1){
 		rep(lambda, ncol(x))
 	}

	p <- ncol(x)
	xx <- cbind(x,-x)

	A <- diag(p)
	A[col(A) == (row(A) + 1)] <- -1
	A <- A[-p, , drop = F]
	AA <- matrix(0, nrow=2 * (p - 1), ncol=2 * p)
	AA[1:(p - 1), 1:p] <- A
	AA[p:(2 * (p - 1)), (p +1 ):(2 * p)] <- A
	AA <- rbind(AA, diag(2 * p))
	Amat <- t(AA)

	emin <- -min(eigen(crossprod(xx))$val)
	emin <- max(emin, 1e-4) #max(.001,emin)
	dmat <- crossprod(xx) + diag(rep(emin, 2 * p))
	dvec <- drop(crossprod(xx, y) - rep(lambda, 2))
	bvec <- rep(0, length=ncol(Amat))

	a <- solve.QP(dmat, dvec, Amat, bvec)$sol
	bp <- a[1:p]
	bn <- a[-(1:p)]

	list(b=bp-bn, bp=bp, bn=bn)
}

#' @title .orderedLasso
#' @description .orderedLasso
#' @keywords internal
#' @param x input matrix (n x p)
#' @param y output vector (length n)
#' @param lambda vector of penalization parameters, containing one or p elements
#' @param intercept if \code{TRUE}, include a model intercept
#' @param beta_pos optional vector containing positive parts of the model coefficients
#' @param beta_neg optional vector containing negative parts of the model coefficients
#' @param method underlying ordered lasso optimization method ('Solve.QP' or 'GG')
#' @param strongly.ordered if \code{TRUE} (\code{FALSE}), use the strongly (weakly) ordered lasso
#' @param standardize if \code{TRUE}, standardize x
#' @param niter number of iterations
#' @param iter.gg maximum number of generalized gradient iterations
#' @param epsilon convergence error tolerance
#' @importFrom stats sd
#' @return ordered lasso solution; ordered lasso coefficients corresponding to the positive part of the weakly ordered solution (\code{bp}), the negative part of the weakly ordered solution (\code{bn}), the weakly ordered solution (\code{beta}), the weakly ordered intercept (\code{b0}), the strongly ordered intercept (\code{b0.ordered}), and the strongly ordered solution (\code{beta.ordered}). 
.orderedLasso <- function(x, y, lambda, intercept=TRUE, beta_pos=NULL, beta_neg=NULL,
	method='Solve.QP', strongly.ordered=FALSE,
	standardize=TRUE, niter=500, iter.gg=100, epsilon=1e-6){

	#cleaned up original code
	#modified for a vector of lambda
	#fixed some problems on beta tolerance
	#removed unused input, fixed default method parameter
	#added sd(y) > 0 vs. sd(y)==0 cases

	stopifnot(nrow(x) == length(y), lambda >= 0)
	stopifnot(class(lambda) == 'numeric' | class(lambda) == 'integer')
	stopifnot(is.finite(x), is.finite(y), is.finite(lambda))
	stopifnot(length(lambda)==1 | length(lambda)==ncol(x))

	if(is.null(beta_pos)){
		beta_pos <- rep(0, ncol(x))
	}
	if(is.null(beta_neg)){
		beta_neg <- rep(0, ncol(x))
	}

	if(standardize){
		stdeviation_inverse <- 1 / apply(x, 2, sd)
	}
	if(intercept){
		mean_y <- mean(y)
		y <- y - mean_y
		mean_x <- apply(x, 2, mean)
	}
	x <- scale(x, intercept, standardize)

	if(sd(y) > 0){
		if(method == 'GG'){
			ii <- 0
			go <- TRUE
			while(go & ii<niter){
				ii <- ii + 1
				beta_pos_old <- beta_pos
				beta_neg_old <- beta_neg

				r <- y + x %*% beta_neg
				beta_pos <- .fastgg(x, r, beta_pos, lambda=lambda, inneriter.gg=iter.gg, epsilon=epsilon)
				r <- y - x %*% beta_pos
				beta_neg <- .fastgg(-x, r, beta_neg, lambda=lambda, inneriter.gg=iter.gg, epsilon=epsilon)

				go <- sum(abs(beta_pos_old-beta_pos)) > epsilon | sum(abs(beta_neg_old - beta_neg)) > epsilon
			}
		} else{ #Solve.QP
			update <- .ordLas2(x, y, lambda)
			beta_pos <- update$bp
			beta_neg <- update$bn
		}

	} else{
		beta_pos <- rep(0, ncol(x))
		beta_neg <- rep(0, ncol(x))
	}

	beta <- beta_pos - beta_neg
	beta.diff <- -diff(abs(beta)) #mistake in original code: abs missing
	if(ncol(x) == 1){
		ordered.test <- FALSE
	} else if(any(beta.diff < -1e-5)){ #allow |b[k]| to be at less than |b[k+1]| by at most 1e-5
		ordered.test <- TRUE
	} else{
		ordered.test <- FALSE
	}

	if(strongly.ordered & ordered.test){
		signvec <- sign(beta)
		signvec[signvec==0] <- sample(c(-1,1), sum(signvec==0), TRUE) #not in original code; see note in .timeLagEstOrdered
		beta.ordered <- .ordLasSignPos(x=x, y=y, lambda=lambda, signvec=signvec)
	} else if(strongly.ordered & !ordered.test){
		beta.ordered <- beta
	} else{ #!strongly.ordered
		beta.ordered <- b0.ordered <-  NULL
	}

	if(standardize){
		beta_pos <- stdeviation_inverse * beta_pos
		beta_neg <- stdeviation_inverse * beta_neg
		beta <- beta_pos - beta_neg
		if(strongly.ordered){
			beta.ordered <- stdeviation_inverse * beta.ordered
		}
	}

	if(intercept){
		b0 <- as.numeric(mean_y - mean_x %*% beta)
		if(strongly.ordered){
			b0.ordered <- as.numeric(mean_y - mean_x %*% beta.ordered)
		}
	} else{
		b0 <- NULL
		if(strongly.ordered){
			b0.ordered <- NULL
		}
	}

	list(bp=beta_pos, bn=beta_neg, beta=beta, b0=b0, b0.ordered=b0.ordered, beta.ordered=beta.ordered)
}

###################
###################

#' @title .timeLagLassoEstOrdered
#' @description .timeLagLassoEstOrdered
#' @keywords internal
#' @param xData lagged matrix (n x (p * maxLag))
#' @param yData output vector (length n)
#' @param lambda vector of penalization parameters, containing one or p elements
#' @param maxLag maximum predictor lag
#' @param intercept if \code{TRUE}, include a model intercept
#' @param beta_pos optional vector containing positive parts of the model coefficients
#' @param beta_neg optional vector containing negative parts of the model coefficients
#' @param stdeviation_inverse_scaled vector of inverse standard deviations of the prescaled xData
#' @param standardize if \code{TRUE}, unscale the coefficients
#' @param method underlying ordered lasso optimization method ('Solve.QP' or 'GG')
#' @param strongly.ordered if \code{TRUE} (\code{FALSE}), use the strongly (weakly) ordered lasso
#' @param maxiter maximum number of  time lagged ordered lasso iterations
#' @param inneriter maximum number ordered lasso iterations
#' @param iter.gg maximum number of generalized gradient iterations
#' @param epsilon convergence error tolerance
#' @importFrom stats sd
#' @return time-lagged ordered lasso solution; time-lagged ordered lasso coefficients corresponding to the positive part of the weakly ordered solution (\code{bp}), the negative part of the weakly ordered solution (\code{bn}), the weakly ordered solution (\code{beta}), the weakly ordered intercept (\code{b0}), the strongly ordered intercept (\code{b0.ordered}), and the strongly ordered solution (\code{beta.ordered}). 
.timeLagLassoEstOrdered <- function(x, y, lambda, maxlag, intercept, beta_pos, beta_neg, stdeviation_inverse_scaled, standardize,
	method, strongly.ordered, maxiter=500, inneriter=100, iter.gg=100, epsilon=1e-6){

	#cleaned up original code
	#modified for a vector of lambda
	#fixed some problems on beta tolerance
	#removed unused input, fixed default method parameter
	#added sd(y) > 0 vs. sd(y)==0 cases

	p <- ncol(x) / maxlag
	beta_pos_new <- beta_pos
	beta_neg_new <- beta_neg

	if(intercept){
		mean_x <- apply(x, 2, mean)
		x <- scale(x, TRUE, FALSE)
		mean_y <- mean(y)
		y <- y - mean_y
	}

	val_new <- val <- .ObjMinFun(x=x, y=y, bp=beta_pos_new, bn=beta_neg_new, lambda=lambda)

	if(sd(y) > 0){
		go <- TRUE
		ii <- 0
		while(go & ii < maxiter){
			ii <- ii + 1
			for (j in seq(0, p-1)){
				subsetindex <- (j * maxlag + 1) : (j * maxlag + maxlag)

				r_j <- y  - (x[, -subsetindex, drop=FALSE] %*% (beta_pos_new[-subsetindex, drop=FALSE] - beta_neg_new[-subsetindex, drop=FALSE]))
				beta_cal <- .orderedLasso(x[, subsetindex], r_j, lambda[subsetindex], intercept=FALSE,
					beta_pos=beta_pos_new[subsetindex], beta_neg=beta_neg_new[subsetindex],
					method=method, strongly.ordered=FALSE, standardize=FALSE,
					niter=inneriter, iter.gg=iter.gg,  epsilon=epsilon)
				beta_pos_new[subsetindex] <- beta_cal$bp
				beta_neg_new[subsetindex] <- beta_cal$bn
			}

			val_new <- .ObjMinFun(x=x, y=y, bp=beta_pos_new, bn=beta_neg_new, lambda=lambda)
			go <- abs(val_new - val) > abs(val * epsilon)
			val <- val_new
		}
		beta_pos <- beta_pos_new
		beta_neg <- beta_neg_new
	} else{
		beta_pos <- rep(0, length(beta_pos))
		beta_neg <- rep(0, length(beta_neg))
	}

	beta <- beta_pos - beta_neg
	if (intercept){
		b0 <- as.numeric(mean_y - mean_x %*% beta)
	} else{
		b0 <- NULL
	}

  ####strongly.ordered part#########################################################################################################
	b0.ordered <- b0
	for(t in seq(0, p-1)){
		index <- (t * maxlag+1):(t*maxlag+maxlag)
		beta.diff <- -diff(abs(beta[index])) #abs missing in original code
		if(maxlag == 1){
			ordered.test <- FALSE
		} else if(any(beta.diff < -1e-5)){
			ordered.test <- TRUE
			break
		} else{
			ordered.test <- FALSE
		}
	}

	if(strongly.ordered){
		beta_ordered <- beta
		if(ordered.test){
			signvec <- sign(beta_ordered)
			val_ordered <- .ObjMinSignFun(x=x, y=y, b0.ordered=0, beta.ordered=beta_ordered, lambda=lambda, sign=signvec)

			jj <- 1
			go <- TRUE

			while(go & jj < maxiter){
				signvec <- sign(beta_ordered) #not in the original code; include because the signs may have changed in the previous iteration!
				signvec[signvec==0] <- sample(c(-1,1), sum(signvec==0), TRUE) #not in original code; but,
				#putting it in allows lags to possibly change from 0 quickly and not get stuck in local minima
				#otherwise, the algorithm has to wait for earlier lags to change from 0 first
				#if sign(b_i)=0, then the old code would force all b_j, j > i to 0 in the next iteration.
				#we want sparsity through the lasso penalization, not through some quirk with the optimization algorithm!
				jj <- jj + 1
				go <- FALSE
				for (j in seq(0, p-1)){
					subsetindex <- (j * maxlag + 1) : (j * maxlag + maxlag)
					r_j <- y - (x[, -subsetindex] %*% (beta_ordered[-subsetindex]))
					beta_ordered[subsetindex] <- .ordLasSignPos(x=x[, subsetindex], y=r_j, lambda=lambda[subsetindex], signvec=signvec[subsetindex])
				}

				val_ordered_new <- .ObjMinSignFun(x=x, y=y, b0.ordered=0, beta.ordered=beta_ordered, lambda=lambda, sign=signvec)
				go <- abs(val_ordered_new - val_ordered) > abs(val * epsilon)
				val_ordered <- val_ordered_new
			}
			if(intercept){
				b0.ordered <- as.numeric(mean_y - mean_x %*% beta_ordered)
			}

		}
	} else{
		beta_ordered <- NULL
	}

	if(standardize){
		beta <- stdeviation_inverse_scaled  * beta
		beta_pos <- stdeviation_inverse_scaled * beta_pos
		beta_neg <- stdeviation_inverse_scaled * beta_neg
		if(strongly.ordered){
			beta_ordered <- stdeviation_inverse_scaled * beta_ordered
		}
	}

	list(bp=beta_pos, bn=beta_neg, beta=beta, b0=b0, b0.ordered=b0.ordered, beta.ordered=beta_ordered)
}
