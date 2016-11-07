#' List Covariance Matrices for Every Subject
#'
#' This function is generic, and methods exists for character, list, function,
#' and various model fit classes.
#'
#' @param obj an R object of class character, function, or a model fit
#' @param t list of vectors of observation times, one for each subject
#' @return a list containing covariance matrices of appropriate dimensions
list_covariances <- function(obj, t) {
	UseMethod("list_covariances", obj)
}

#' @rdname list_covariances
list_covariances.character <- function(obj, t) {
	switch(obj,
		identity = lapply(t, function(ti) {
			diag(length(ti))
		}),
		randomwalk = lapply(t, function(ti) {
			mi <- length(ti)
			outer(1:mi, 1:mi, pmin)
		}),
		brownian = lapply(t, function(ti) {
			outer(ti, ti, pmin)
		}),
		pascal = lapply(t, function(ti) {
			mi <- length(ti)
			outer(1:mi, 1:mi, function(i, j) {
				choose(i + j - 1, i - 1)
			})
		}))
}

#' @rdname list_covariances
list_covariances.list <- function(obj, t) {
	obj
}

#' @rdname list_covariances
list_covariances.function <- function(obj, t) {
	lapply(t, function(ti) {
		outer(ti, ti, obj)
	})
}

#' @rdname list_covariances
list_covariances.jmcmMod <- function(obj, t) {
	lapply(1:length(t), function(i) jmcm::getJMCM(obj, "Sigma", i))
}

#' @rdname list_covariances
list_covariances.lmerMod <- function(obj, t) {
	i <- rep(names(t), sapply(t, length))
	i <- factor(i, levels = unique(i)) # to preserve the order in the data
	mmList <- lme4::getME(obj, "mmList")
	zList <- lapply(mmList, function(mm) split.data.frame(mm, i))
	vList <- unclass(lme4::VarCorr(obj))
	zvzList <- Map(function(z, v) {
		lapply(z, function(zi) {
			zi %*% v %*% t(zi)
		})}, zList, vList)
	zvzSumList <- Reduce(function(x, y) mapply("+", x, y), zvzList) 
	resVar <- attr(vList, "sc")^2 
	lapply(zvzSumList, function(zvzSum) {
		zvzSum + diag(resVar, nrow = nrow(zvzSum))
	})
}
