#' Imputation by partial least square regression
#'
#' Imputes univariate missing data by using partial least square regression.
#'
#' @aliases mice.impute.pls pls
#' @inheritParams mice.impute.norm.boot
#' @param nlvs The number of latent variables to consider in the PLS regression.
#' @return Vector with imputed data, same type as \code{y}, and of length
#' \code{sum(wy)}
#' @details
#' Imputation of \code{y} by PLS. The procedure is as follows:
#' \enumerate{
#' \item A bootstrap sample of the observed part of the data is taken
#' \item PLSR is estimated on this sample for a given number of latent variables (nlvs).
#' \item Based on this PLSR model, predictions on the new data x[!ry, ] are returned
#' \item Imputations are then obtained by adding noise to these predictions.
#' }
#' 
#' Noise is added by sampling from a standard normal distribution rescaled based 
#' on the residuals standard error of the PLSR.
#' The residual degrees of freedom (df) are computed based on Krylov representation 
#' described by Kramer and Sugiyama (2011).
#' 
#' @author Edoardo Costantini, 2022
#' @references
#' 
#' Krämer, N., & Sugiyama, M. (2011). The degrees of freedom of partial least 
#' squares regression. Journal of the American Statistical Association, 106(494), 697-705.
#' 
#' Wehrens, R., & Mevik, B. H. (2007). The pls package: principal component and 
#' partial least squares regression in R.
#'
#' @family univariate imputation functions
#' @keywords datagen
#' @export
mice.impute.pls <- function(y, ry, x, wy = NULL, nlvs = 1L, ...) {

    # Set up    
    install.on.demand("plsdof", ...)
    if (is.null(wy)) wy <- !ry

    # Take bootstrap sample for model uncertainty
    n1 <- sum(ry)
    s <- sample(n1, n1, replace = TRUE)
    dotxobs <- x[ry, , drop = FALSE][s, ]
    dotyobs <- y[ry][s]

    # Train PLS and get predictions
    pls_fit <- plsdof::pls.model(
        X = dotxobs,
        y = dotyobs,
        m = ncol(dotxobs),
        compute.DoF = TRUE,
        Xtest = x[!ry, , drop = FALSE], 
        ytest = NULL
    )

    # Extract residual standard error (sd of residuals with df as denominator) for PLS
    sigma <- pls_fit$sigmahat[nlvs + 1]

    # Add noise for imputation uncertainty
    imputes <- pls_fit$prediction[, nlvs + 1] + rnorm(sum(wy)) * sigma

    # Return
    return(imputes)
}
