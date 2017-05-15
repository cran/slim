#' Methods for Singular Linear Model Fits
#'
#' @name slim.methods
#' @param object an object of class 'slim', usually, a result of a call to
#' 'slim'. 
#' @param empirical logical indicating if empirical variances of y should be
#' used in estimating standard errors (the default). Empirical standard errors
#' should be used unless covariances have been well modelled.
#' @param ... arguments passed to or from other methods.
NULL

#' Extract Model Coefficients from Singular Linear Model
#'
#' @inheritParams slim.methods
#' @return a vector of model coefficients.
#' @export
coef.slim <- function(object, ...) {
	object$coefficients
}

#' Confidence Intervals for Model Parameters from Singular Linear Model
#' 
#' @inheritParams slim.methods
#' @inheritParams stats::confint
#' @return A matrix (or vector) with columns giving lower and upper confidence
#' limits for each parameter.
#' @export
confint.slim <- function(object, parm, level = 0.95, empirical = TRUE, ...) {
	cf <- coef.slim(object)
	pnames <- names(cf)
	if(missing(parm))
		parm <- pnames
	else if(is.numeric(parm))
		parm <- pnames[parm]
	a <- (1 - level) / 2
	a <- c(a, 1 - a)
	fac <- stats::qnorm(a)
	pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE,
		digits = 3), "%")	
	ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm,
			pct))
	ses <- sqrt(diag(vcov.slim(object, empirical)))[parm]
	ci[] <- cf[parm] + ses %o% fac
	ci
}

#' Extract Model Fitted Values from Singular Linear Model
#'
#' @inheritParams slim.methods
#' @return a vector of fitted values from the model fit.
#' @export
fitted.slim <- function(object, ...) {
	object$fitted_values
}

#' Model Predictions from Singular Linear Model
#'
#' @inheritParams slim.methods
#' @inheritParams stats::predict.lm
#' @return a vector of model predictions.
#' @export
predict.slim <- function(object, newdata, ...) {
	f <- stats::formula(object)[-2]
	if (missing(newdata) || is.null(newdata)) {
		mf <- stats::model.frame(object)
	} else {
		mf <- stats::model.frame(f, newdata)
	}
	mm <- stats::model.matrix(f, mf, contrasts.arg = object$contrasts)
	as.vector(mm %*% coef.slim(object))
}

#' Print 'slim' Objects
#'
#' 'print' methods for class 'slim' and 'slim_summary'. 'print.slim_summary'
#' differs only in its default value of 'empirical'.
#'
#' @param x an object of class 'slim' or 'slim_summary', as appropriate.
#' @inheritParams slim.methods
#' @inheritParams base::print
#' @inheritParams stats::summary.lm
#' @return x, invisibly.  
#' @export
print.slim <- function(x, empirical = TRUE,
	digits = max(3, getOption("digits") - 3),
	signif.stars = getOption("show.signif.stars"), ...) {
	cat("Singular Linear Model Fit\n")
	cat("Call:\n", paste(deparse(x$call), sep = "\n",
		collapse = "\n"), "\n\n", sep = "")
	stats::printCoefmat(summary(x, empirical)$coefficient_matrix,
		digits = digits, signif.stars = signif.stars,
		na.print = "NA", ...)
	invisible(x)
}

#' @rdname print.slim
#' @export
print.slim_summary <- function(x, empirical = x$empirical, ...) {
	print.slim(x, empirical, ...)
}

#' Extract Model Residuals from Singular Linear Model 
#' 
#' @inheritParams slim.methods
#' @return a vector of model residuals.
#' @export
residuals.slim <- function(object, ...) {
	object$residuals
}

#' Summarizing Singular Linear Model Fits
#'
#' 'summary' method for class 'slim'.
#'
#' @inheritParams slim.methods
#' @return an object with class c("slim_summary", "slim") and, in addition to
#' the usual 'slim' components, coefficient_matrix (the matrix of estimated
#' coefficients, standard errors, z- and p-values) and empirical (logical
#' indicating if empirical standard errors have been used)
#' @export
summary.slim <- function(object, empirical = TRUE, ...) {
	est <- coef.slim(object)
	se <- sqrt(diag(vcov.slim(object, empirical)))
	zval <- est / se
	coefMat <- cbind(est, se, zval, 2 * stats::pnorm(abs(zval),
			lower.tail = FALSE))
	dimnames(coefMat) <- list(names(est),
		c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
	class(object) <- c("slim_summary", "slim")
	object$coefficient_matrix <- coefMat
	object$empirical <- empirical
	object
}

#' Extract Variance-Covariance Matrix from a 'slim' Object
#'
#' 'vcov' method for class 'slim'.
#'
#' @inheritParams slim.methods
#' @return a matrix of the estimated covariances between the parameter
#' estimates.
#' @export
vcov.slim <- function(object, empirical = TRUE, ...) {
	if(empirical) object$vcov_empirical
	else object$vcov_modelled
}
