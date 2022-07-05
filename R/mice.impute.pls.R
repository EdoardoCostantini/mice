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
#' Imputation of \code{y} by partial least square regression (Wold, 1975).
#' The method consists of the following steps:
#' \enumerate{
#' \item For a given \code{y} variable under imputation, draw a bootstrap version y*
#' with replacement from the observed cases \code{y[ry]}, and stores in x* the
#' corresponding values from \code{x[ry, ]}.
#' \item PLSR is estimated on this sample for a given number of latent variables (nlvs).
#' \item Based on this PLSR model, predictions on the new data x[wy, ] are returned
#' \item Imputations are then obtained by adding noise to these predictions.
#' }
#' 
#' Noise is added by sampling from a standard normal distribution rescaled based 
#' on the residuals standard error of the PLSR model.
#' The residual degrees of freedom are computed based on the Krylov representation 
#' described by Kramer and Sugiyama (2011) as implemented in the R package \code{plsdof} 
#' (Krämer and Sugiyama, 2011).
#' 
#' @author Edoardo Costantini, 2022
#' @references
#' 
#' Krämer, N., & Sugiyama, M. (2011). The degrees of freedom of partial least 
#' squares regression. Journal of the American Statistical Association, 106(494), 697-705.
#' 
#' Wold, H. (1975). Path models with latent variables: The NIPALS approach. 
#' In Quantitative sociology (pp. 307-357). Academic Press.
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
