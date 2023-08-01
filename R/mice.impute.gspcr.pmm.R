#' Imputation by predictive mean matching based on generalized supervised principal component regression
#'
#' Imputes univariate missing data using predictive mean matching and generalized supervised principal component regression.
#'
#' @aliases mice.impute.gspcr.pmm gspcr.pmm
#' @inheritParams mice.impute.pmm
#' @param thrs character vector of length 1 storing the type of threshold to be used (see below for available options)
#' @param nthrs numeric vector of length 1 storing the number of threshold values to be used
#' @param fit_measure character vector of length 1 indicating the type of fit measure to be used in the cross-validation procedure
#' @param npcs_range numeric vector defining the numbers of principal components to be used
#' @param K numeric vector of length 1 storing the number of folds for the K-fold cross-validation procedure
#' @param donors the size of the donor pool among which a draw is made (default is \code{5L}).
#' @param exclude Value or vector of values to exclude from the imputation donor pool in \code{y}
#' @return Vector with imputed data, the same type as \code{y}, and of length
#' \code{sum(wy)}
#' @details
#' For details on the gspcr framework, see [mice::mice.impute.gspcr.norm()].
#' For details on how the pmm framework see [mice::mice.impute.pmm()].
#' Note that \code{gspcr.pmm} uses match type 2, based on the bootstrap pmm presented in Heitjan and Little (1991).
#'
#' @author Edoardo Costantini, 2022
#' @references
#'
#' Bair, E., Hastie, T., Paul, D., & Tibshirani, R. (2006). Prediction by
#' supervised principal components. Journal of the American Statistical
#' Association, 101(473), 119-137.
#' 
#' Heitjan, D. F., & Little, R. J. (1991). Multiple imputation for the fatal accident reporting system. Journal of the Royal Statistical Society Series C: Applied Statistics, 40(1), 13-29.
#'
#' @family univariate imputation functions
#' @keywords imputation
#' @export
mice.impute.gspcr.pmm <- function(y, ry, x, wy = NULL,
                                  thrs = "PR2",
                                  fit_measure = "BIC",
                                  nthrs = 10,
                                  npcs_range = 1:3,
                                  K = 1,
                                  donors = 5L,
                                  exclude = -99999999,
                                  ...) {

    # Prepare data in pmm fashion ----------------------------------------------

    # Create an id vector for exclusion
    id.ex <- !ry | !y %in% exclude

    # leave out the excluded vector y's
    y <- y[id.ex] 

    # allow for one-dimensional x-space
    if (!is.null(dim(x))) {
        x <- x[id.ex, ]
    } else {
        x <- x[id.ex]
    }

    # leave out the exclude vector x's
    ry <- ry[id.ex]

    # If applicable adjust wy to match and exclude
    if (is.null(wy)) {
        wy <- !ry
    } else {
        wy <- wy[id.ex]
    }

    # Transform factor y to integer
    ynum <- y
    if (is.factor(y)) {
        ynum <- as.integer(y)
    }

    # Bootstrap sample for model uncertainty -----------------------------------

    # Sample size of responses
    n1 <- sum(ry)

    # Define bootstrap sample
    s <- sample(n1, n1, replace = TRUE)

    # Create bootstrap sample observed predictors
    dotxobs <- x[ry, , drop = FALSE][s, ]

    # Create bootstrap sample of observed values of variable under imputation
    dotyobs <- ynum[ry][s]

    # Train GSPCR --------------------------------------------------------------
    gscpr_fit <- gspcr::cv_gspcr(
        dv = dotyobs,
        ivs = dotxobs,
        fam = "gaussian",
        thrs = thrs,
        nthrs = nthrs,
        npcs_range = npcs_range,
        K = K,
        fit_measure = fit_measure,
        min_features = 1
    )

    # Estimate GSPCR -----------------------------------------------------------
    gspcr_est <- gspcr::est_gspcr(gscpr_fit)

    # Obtain imputations -------------------------------------------------------

    # Obtain predictions for the observed values of y
    y_hat_obs <- predict(
        object = gspcr_est,
        newdata = x[ry, ]
    )

    # Obtain predictions for the unobserved values of y
    if(sum(wy) > 0){
        y_hat_mis <- predict(
            object = gspcr_est,
            newdata = x[wy, ]
        )
    } else {
        y_hat_mis <- matrix(nrow = 0, ncol = 1)
    }

    # Use matcher
    idx <- matchindex(y_hat_obs, y_hat_mis, donors)

    # Define imputations
    imputes <- y[ry][idx]

    # Return imputations
    return(imputes)
}
