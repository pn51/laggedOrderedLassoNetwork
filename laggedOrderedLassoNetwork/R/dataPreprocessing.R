
#data preprocessing functions

#' @title .rescaleData
#' @description .rescaleData
#' @keywords internal
#' @param datasets list of datasets (samples x variables)
#' @return list of datasets with variables scaled by the overall variable means and standard deviations across all input datasets
.rescaleData <- function(datasets){
	# given a list of datasets, combine
	# and rescale each row to mean 0, sd 1
	# and separate
	genes <- colnames(datasets[[1]])
	if(is.list(datasets)){
		datasets <- lapply(datasets, function(ii){ii[,genes]})
		allData <- do.call(rbind, datasets) #assuming genes x time
		samples <- sapply(datasets, nrow)
		samples <- rep(seq_along(samples), samples)
		allData <- scale(allData)
		lapply(unique(samples), function(ii){
			allData[samples==ii, ]
		})
	} else{
		scale(datasets)
	}
}

#' @title .rescaleDataSeparate
#' @description .rescaleDataSeparate
#' @keywords internal
#' @param datasets list of datasets (samples x variables)
#' @return list of datasets, each separately scaled
.rescaleDataSeparate <- function(datasets){
	# given a list of datasets, treat the datasets separately
	# and rescale each row to mean 0, sd 1
	lapply(datasets, function(ii){
		scale(ii)
	})
}
