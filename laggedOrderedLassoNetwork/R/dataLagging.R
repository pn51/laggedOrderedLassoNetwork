
#functions for lagging data and constructing expression and change-in-expression output variables

#' @title .transformMultiLag
#' @description .transformMultiLag
#' @keywords internal
#' @param exprData gene expression matrix (timepoints x genes)
#' @param maxLag maximum predictor lag
#' @return a list containing two matrices (\code{xData}), (\code{yData}) containing the lagged expression data and the output expression data at time t+1, respectively
.transformMultiLag <- function(exprData, maxLag=1){
	# expression at time t+1
	yData <- exprData[-(1:maxLag),]
	k <- nrow(yData)

	if(k > 0){
		genes <- colnames(exprData)

		# lag xData
		if(maxLag > 1){
			xData <- lapply(maxLag:1, function(ii){
				laggedMatrix <- exprData[(ii):(ii+k-1),]
				colnames(laggedMatrix) <- paste(genes,maxLag-ii+1,sep='-')
				laggedMatrix
			})
			xData <- do.call(cbind, xData)
		} else{
			xData <- exprData[1:k,]
			colnames(xData) <- paste(genes, 1, sep='-')
		}
		#reorder so the each gene's lags are blocked together and ordered by increasing lag
		reordered <- paste(rep(genes, each=maxLag), rep(seq(maxLag), ncol(exprData)), sep='-')
		xData <- xData[,reordered]

		lags <- rep(1:maxLag, length(genes))
		genes <- rep(genes, maxLag)

		rownames(yData) <- NULL
		rownames(xData) <- NULL

		list(xData=xData, yData=yData)
	} else{
		list(xData=NULL, yData=NULL)
	}

}

#' @title .transformListMultiLag
#' @description .transformListMultiLag
#' @keywords internal
#' @param exprDataList list of gene expression matrices (timepoints x genes)
#' @param maxLag maximum predictor lag
#' @return a list containing two matrices (\code{xData}), (\code{yData}) containing the lagged expression data and the output expression data at time t+1, respectively
.transformListMultiLag <- function(exprDataList, maxLag=1){
	#apply transformMultiLag to a list of expr. datasets with the same set of genes in the same order
	genes <- colnames(exprDataList[[1]])
	transformed <- lapply(exprDataList, function(curExpr){
		.transformMultiLag(curExpr[, genes], maxLag)
	})

	list(xData=do.call(rbind, lapply(transformed, function(ii){ii$xData})),
		yData=do.call(rbind, lapply(transformed, function(ii){ii$yData})))
}

#' @title .transformMultiLagChange
#' @description .transformMultiLagChange
#' @keywords internal
#' @param exprData gene expression matrix (timepoints x genes)
#' @param maxLag maximum predictor lag
#' @return a list containing two matrices (\code{xData}), (\code{yData}) containing the lagged expression data and the output change-in-expression data at time t+1, respectively
.transformMultiLagChange <- function(exprData, maxLag=1){
	# change in expression at time t+1
	yData <- exprData[(maxLag+1):nrow(exprData),] - exprData[maxLag:(nrow(exprData)-1),]
	k <- nrow(yData)

	if(k > 0){
		genes <- colnames(exprData)

		# lag xData
		if(maxLag > 1){
			xData <- lapply(maxLag:1, function(ii){
				laggedMatrix <- exprData[(ii):(ii+k-1),]
				colnames(laggedMatrix) <- paste(genes,maxLag-ii+1,sep='-')
				laggedMatrix
			})
			xData <- do.call(cbind, xData)
		} else{
			xData <- exprData[1:k,]
			colnames(xData) <- paste(genes, 1, sep='-')
		}
		#reorder so the each gene's lags are blocked together and ordered by increasing lag
		reordered <- paste(rep(genes, each=maxLag), rep(seq(maxLag), ncol(exprData)), sep='-')
		xData <- xData[,reordered]

		lags <- rep(1:maxLag, length(genes))
		genes <- rep(genes, maxLag)

		rownames(yData) <- NULL
		rownames(xData) <- NULL

		list(xData=xData, yData=yData)
	} else{
		list(xData=NULL, yData=NULL)
	}

}

#' @title .transformListMultiLagChange
#' @description .transformListMultiLagChange
#' @keywords internal
#' @param exprDataList list of gene expression matrices (timepoints x genes)
#' @param maxLag maximum predictor lag
#' @return a list containing two matrices (\code{xData}), (\code{yData}) containing the lagged expression data and the output change-in-expression data at time t+1, respectively
.transformListMultiLagChange <- function(exprDataList, maxLag=1){
	#apply transformMultiLagChange to a list of expr. datasets with the same set of genes in the same order
	genes <- colnames(exprDataList[[1]])
	transformed <- lapply(exprDataList, function(curExpr){
		.transformMultiLagChange(curExpr[, genes], maxLag)
	})

	list(xData=do.call(rbind, lapply(transformed, function(ii){ii$xData})),
		yData=do.call(rbind, lapply(transformed, function(ii){ii$yData})))
}
