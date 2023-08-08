#' Imputation by generalized supervised principal component logistic regression
#'
#' Imputes univariate missing values on a binary variable using generalized supervised principal component regression.
#'
#' @aliases mice.impute.gspcr.logreg gspcr.logreg
#' @inheritParams mice.impute.norm.boot
#' @param thrs character vector of length 1 storing the type of threshold to be used (see below for available options)
#' @param nthrs numeric vector of length 1 storing the number of threshold values to be used
#' @param fit_measure character vector of length 1 indicating the type of fit measure to be used in the cross-validation procedure
#' @param npcs_range numeric vector defining the numbers of principal components to be used
#' @param K numeric vector of length 1 storing the number of folds for the K-fold cross-validation procedure
#' @return Vector with imputed data, the same type as \code{y}, and of length
#' \code{sum(wy)}
#' @details
#' See [mice.impute.gspcr.logreg].
#'
#' @author Edoardo Costantini, 2022
#' @references
#'
#' Bair, E., Hastie, T., Paul, D., & Tibshirani, R. (2006). Prediction by
#' supervised principal components. Journal of the American Statistical
#' Association, 101(473), 119-137.
#'
#' @family univariate imputation functions
#' @keywords imputation
#' @export
mice.impute.gspcr.logreg <- function(y, ry, x, wy = NULL,
                                     thrs = "PR2",
                                     fit_measure = "BIC",
                                     nthrs = 10,
                                     npcs_range = 1:5,
                                     K = 1,
                                     ...) {
    # Set up -------------------------------------------------------------------

    # Create wy if not there
    if (is.null(wy)) wy <- !ry

    # Bootstrap sample for model uncertainty -----------------------------------

    # Sample size of responses
    n1 <- sum(ry)

    # Define bootstrap sample
    s <- sample(n1, n1, replace = TRUE)

    # Create bootstrap sample observed predictors
    dotxobs <- x[ry, , drop = FALSE][s, ]

    # Create bootstrap sample of observed values of variable under imputation
    dotyobs <- y[ry][s]

    # escape with same input if the dependent does not vary
    cat.has.all.obs <- table(dotyobs) == sum(ry)
    if (any(cat.has.all.obs)) {
        return(y[ry])
    }

    # GSPCR --------------------------------------------------------------------

    # Train model to tune parameters
    gscpr_fit <- gspcr::cv_gspcr(
        dv = dotyobs,
        ivs = dotxobs,
        fam = "binomial",
        thrs = thrs,
        nthrs = nthrs,
        npcs_range = npcs_range,
        K = K,
        fit_measure = fit_measure,
        min_features = 1
    )

    # Estimate model with tuned parameters
    gspcr_est <- gspcr::est_gspcr(gscpr_fit)

    # Obtain imputations -------------------------------------------------------

    # Obtain predictions
    p_hat <- predict(
        object = gspcr_est,
        newdata = x[!ry, ]
    )

    # Draw imputations
    imputes <- runif(sum(wy)) <= p_hat

    # Assign 0 or 1 as values
    imputes[imputes] <- 1

    # Return as factor if original was factor
    if (is.factor(y)) {
        imputes <- factor(imputes, c(0, 1), levels(y))
    }

    # Return imputations
    return(imputes)
}