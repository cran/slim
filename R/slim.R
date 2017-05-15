#' Singular linear models for longitudinal data.
#'
#' The slim package fits singular linear models to longitudinal data. Singular
#' linear models are useful when the number, or timing, of longitudinal
#' observations may be informative about the observations themselves. They are
#' described in Farewell (2010) <doi:10.1093/biomet/asp068>, and are extensions
#' of the linear increments model of Diggle et al. (2007)
#' <doi:10.1111/j.1467-9876.2007.00590.x> to general longitudinal data.
#'
#' The most important function is slim, whose formula interface is similar to
#' that of lm.
#' @docType package
#' @name slim-package
#' @seealso \code{\link{slim}}
#' @import data.table
NULL

#' Laurent Expansion of Inverse of Linear Matrix Function 
#'
#' This function computes the first two terms of the Laurent expansion of the
#' inverse of a linear matrix function.
#'
#' @param V for some integer m >= 1, an array of dimension (m, m, 2), where V[,
#' , 1] is the intercept and V[, , 2] is the slope of the linear matrix
#' function.
#' @param zapsmall logical: should zapsmall be called on the result? Default TRUE.
#' @return array of dimension (m, m, 2), where W[, , 1] corresponds to the
#' exponent -1, and W[, , 2] corresponds to the exponent 0.
compute_laurent <- function(V, zapsmall = TRUE) {
	m <- dim(V)[1]
	zero_matrix <- matrix(0, nrow = m, ncol = m)
	augmented_matrix <- rbind(cbind(V[, , 1], zero_matrix),
		cbind(V[, , 2], V[, , 1]))
	G <- MASS::ginv(augmented_matrix)
	G00 <- G[1:m, 1:m]
	W0 <- G01 <- G[1:m, (m + 1):(2 * m)]
	W1 <- G00 %*% (diag(m) - V[, , 2] %*% G01)
	W <- c(W0, W1)
	if(zapsmall) W <- zapsmall(W)
	array(W, dim = c(m, m, 2))
}

#' Fitter Function for Singular Linear Models
#'
#' This function computes the limiting solution to the estimating equation
#' sum(x' V^{-1} (y - x beta)) = 0 as the covariance V tends from V[, , 1] +
#' V[, , 2] to V[, , 1]. 
#'
#' @param x list of design matrices, one for each subject, all having the same
#' number of columns.
#' @param V list of covariance arrays, one for each subject, matching the
#' dimensions of y.
#' @param y list of response vectors, one for each subject.
#' @return a list with components coefficients (the limiting solution),
#' residuals, fitted_values, vcov_empirical and vcov_modelled. 
fit_slim <- function(x, V, y) {
	p <- dim(x[[1]])[2]
	W <- lapply(V, compute_laurent)
	Q <- Map(function(x, W) {
		Q0 <- t(x) %*% W[, , 1] %*% x
		Q1 <- t(x) %*% W[, , 2] %*% x
		array(c(Q0, Q1), dim = c(p, p, 2))
	}, x, W)
	R <- Map(function(x, W, y) {
		R0 <- t(x) %*% W[, , 1] %*% y
		R1 <- t(x) %*% W[, , 2] %*% y
		array(c(R0, R1), dim = c(p, 1, 2))
	}, x, W, y)
	Qsum <- Reduce('+', Q)
	Rsum <- Reduce('+', R)
	Qaug <- rbind(cbind(Qsum[, , 1], matrix(0, nrow = p, ncol = p)),
			cbind(Qsum[, , 2], Qsum[, , 1]))
	Raug <- c(Rsum[, , 1], Rsum[, , 2])
	# solve(qr(Qaug, LAPACK = TRUE), Raug)[1:p]
	G <- MASS::ginv(Qaug)
	beta <- (G %*% Raug)[1:p]
	pred <- lapply(x, function(x) as.vector(x %*% beta))
	resd <- Map(function(y, pred) y - pred, y, pred) 
	covy <- list(empirical = lapply(resd, function(r) outer(r, r)),
		modelled = lapply(V, function(V) V[, , 1] + V[, , 2]))
	covMat <- lapply(covy, function(covy) {
		VarR <- Map(function(x, W, covy) {
			mult <- cbind(W[, , 1] %*% x, W[, , 2] %*% x)
			t(mult) %*% covy %*% mult
		}, x, W, covy)
		VarRSum <- Reduce('+', VarR)
		(G %*% VarRSum %*% t(G))[1:p, 1:p]
	})
	list(coefficients = beta,
		residuals = unlist(resd),
		fitted_values = unlist(pred),
		vcov_empirical = covMat$empirical,
		vcov_modelled = covMat$modelled)
}

#' Fit Singular Linear Models
#'
#' Fit a singular linear model to longitudinal data.
#'
#' @param formula a model formula for the fixed effects
#' @param data a 'data.table' with two keys, respectively identifying subjects
#' and observation times
#' @param covariance an R object for which a 'list_covariances' method exists.
#' Options include a character string such as "identity", "randomwalk" (the
#' default), "brownian" or "pascal"; a list of covariance matrices; a function
#' to be used in 'outer' and applied to the observation times; or a 'jmcmMod'
#' or 'lmerMod' model fit.
#' @param limit a one-sided model formula for the (thin) Cholesky factor of the
#' limiting covariance matrix (default ~ 1, so the limiting covariance matrix
#' is the matrix of ones)
#' @param contrasts an optional list. See the 'contrasts.arg' argument of
#' 'model.matrix.default'.
#' @return an object of class 'slim'
#' @examples
#' slim_fit <- slim(renalfn ~ group + month, dialysis)
#' summary(slim_fit)
#'
#' if(require("lme4")) {
#'   lmer_fit <- lmer(renalfn ~ group + month + (1 + month | id), dialysis)
#'   slim_fit <- slim(renalfn ~ 1 + group + month, dialysis, covariance = lmer_fit)
#'   summary(slim_fit)
#'   summary(slim_fit, empirical = FALSE)
#' }
#'
#' if(require("jmcm")) {
#'   jmcm_fit <- jmcm(renalfn | id | month ~ group | 1, dialysis,
#'     triple = rep(2L, 3), cov.method = "mcd")
#'   slim_fit <- slim(renalfn ~ group + month, dialysis, covariance = jmcm_fit)
#'   summary(slim_fit)
#'   summary(slim_fit, empirical = FALSE)
#' }
#' @export
slim <- function(formula, data, covariance = "randomwalk", limit = ~ 1, 
	contrasts = NULL) {
	cl <- match.call()
  if(!haskey(data)) stop("\'data\' must be a *keyed* \'data.table\'.")
	i <- data[[key(data)[1]]]
	i <- factor(i, levels = unique(i)) # to preserve the order in the data
	t <- split(data[[key(data)[2]]], i)
	m <- lapply(t, length)
	mf <- lapply(c(x = formula, z = limit), function(f) {
		stats::lm(f, data, method = "model.frame")
	})
	tm <- lapply(mf, function(fr) {
		stats::terms(fr)
	})
	mm <- Map(function(a, b) {
		stats::model.matrix(a, b, contrasts)
	}, tm, mf)
	x <- split.data.frame(mm$x, i)
	y <- split(stats::model.response(mf$x), i) 
	z <- split.data.frame(mm$z, i)
	V0 <- lapply(z, function(zi) zi %*% t(zi)) 
	Vstar <- list_covariances(covariance, t)
	V <- Map(function(m, V0, Vstar) {
		array(c(V0, Vstar - V0), dim = c(m, m, 2))
	}, m, V0, Vstar)
	out <- fit_slim(x, V, y)
	class(out) <- "slim"
	out$contrasts <- contrasts 
	out$call <- cl
	out$terms <- tm$x
	names(out$coefficients) <- colnames(mm$x)
	rownames(out$vcov_empirical) <- colnames(out$vcov_empirical) <-
		rownames(out$vcov_modelled) <- colnames(out$vcov_modelled) <-
		colnames(mm$x)
	out
}
