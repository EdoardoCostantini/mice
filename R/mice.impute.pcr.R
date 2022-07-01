#' Imputation by principal component regression
#'
#' Imputes univariate missing data by extracting PCs from all potential auxiliary variables.
#'
#' @aliases mice.impute.pcr pcr
#' @inheritParams mice.impute.norm.boot
#' @param npcs The number of principal components to extract for PC regression.
#' @return Vector with imputed data, same type as \code{y}, and of length
#' \code{sum(wy)}
#' @details
#' Extracts components from x and uses them as predictors in a regular call of
#' norm.boot().
#' @author Edoardo Costantini, 2022
#' @references
#'
#' Jolliffe, I. (2002) Principal Component Analysis (2nd ed). Springer.
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

    # Define sigma
    RSS <- sqrt(sum(pcr_out$residuals^2))
    sigma <- RSS / (n1 - npcs - 1)

    # Get prediction on missing part
    yhat <- predict(pcr_out, newdata = x[wy, ], ncomp = npcs, type = "response")

    # Add noise for imputation uncertainty
    imputes <- yhat + rnorm(sum(wy)) * sigma

    # Return
    return(imputes)
}
