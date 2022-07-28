#' Imputation by principal component regression
#'
#' Imputes univariate missing data by extracting PCs from all potential auxiliary variables.
#'
#' @aliases mice.impute.pcr pcr
#' @inheritParams mice.impute.norm.boot
#' @param npcs The number of principal components to extract in the PC regression.
#' @return Vector with imputed data, same type as \code{y}, and of length
#' \code{sum(wy)}
#' @details
#' Imputation of \code{y} by principal component regression (Jolliffe, 1986, p. 168).
#' The method consists of the following steps:
#' \enumerate{
#' \item For a given \code{y} variable under imputation, draw a bootstrap version y*
#' with replacement from the observed cases \code{y[ry]}, and stores in x* the
#' corresponding values from \code{x[ry, ]}.
#' \item Fit a PC regression with y* as the outcome, and x* as predictors using 
#' \code{npcs} components.
#' \item Calculate the estimated residual standard error \code{sigma} based on the residuals 
#' obtained from the PC regression and \code{npcs - 1} degrees of freedom.
#' \item Obtain predicted values for \code{y} based on the fitted PC regression 
#' and the new data \code{x[wy, ]} 
#' \item Obtain imputations by adding noise scaled by \code{sigma} to these
#' predictions.
#' }
#' 
#' @author Edoardo Costantini, 2022
#' @references
#'
#' Jolliffe, I. (1986) Principal Component Analysis. Springer.
#'
#' @family univariate imputation functions
#' @keywords datagen
#' @export
mice.impute.pcr <- function(y, ry, x, wy = NULL, npcs = 1L, ...) {

    # Set up
    install.on.demand("pls", ...)
    if (is.null(wy)) wy <- !ry

    # Take bootstrap sample for model uncertainty
    n1 <- sum(ry)
    s <- sample(n1, n1, replace = TRUE)
    dotxobs <- x[ry, , drop = FALSE][s, ]
    dotyobs <- y[ry][s]

    # Train PCR on dotxobs sample
    pcr_out <- pls::pcr(
        dotyobs ~ dotxobs,
        ncomp = npcs,
        scale = TRUE,
        center = TRUE,
        validation = "none"
    )

    # Compute the residual sum of squares
    res_ss <- sum(stats::resid(pcr_out)[, , npcs]^2)

    # Compute degrees of freedom
    res_df <- n1 - npcs

    # Compute sigma
    sigma <- sqrt(res_ss / res_df)

    # Get prediction on missing part
    yhat <- predict(pcr_out, newdata = x[wy, ], ncomp = npcs, type = "response")

    # Add noise for imputation uncertainty
    imputes <- yhat + rnorm(sum(wy)) * sigma

    # Return
    return(imputes)
}
