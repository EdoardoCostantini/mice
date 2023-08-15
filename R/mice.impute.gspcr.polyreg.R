#' Imputation by generalized supervised principal component polytomous regression.
#'
#' Imputes univariate missing values on an unorded categorical (nominal) variable using generalized supervised principal component proportional odds regression.
#'
#' @aliases mice.impute.gspcr.polyreg gspcr.polyreg
#' @inheritParams mice.impute.norm.boot
#' @param thrs character vector of length 1 storing the type of threshold to be used (see below for available options)
#' @param nthrs numeric vector of length 1 storing the number of threshold values to be used
#' @param fit_measure character vector of length 1 indicating the type of fit measure to be used in the cross-validation procedure
#' @param npcs_range numeric vector defining the numbers of principal components to be used
#' @param K numeric vector of length 1 storing the number of folds for the K-fold cross-validation procedure
#' @return Vector with imputed data, the same type as \code{y}, and of length
#' \code{sum(wy)}
#' @details
#' See [mice.impute.gspcr.norm].
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
mice.impute.gspcr.polyreg <- function(y, ry, x, wy = NULL,
                                      thrs = "PR2",
                                      fit_measure = "BIC",
                                      nthrs = 10,
                                      npcs_range = 1:3,
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
    dotyobs <- droplevels(as.factor(y[ry][s]))

    # escape with same impute if the dependent does not vary
    cat.has.all.obs <- table(dotyobs) == sum(ry)
    if (any(cat.has.all.obs)) {
        return(y[ry])
    }

    # GSPCR --------------------------------------------------------------------
    tryCatch(
        expr = {
            # Train model to tune parameters
            gscpr_fit <- gspcr::cv_gspcr(
                dv = dotyobs,
                ivs = dotxobs,
                fam = "baseline",
                thrs = thrs,
                nthrs = nthrs,
                npcs_range = npcs_range,
                K = K,
                fit_measure = fit_measure,
                min_features = 1,
                save_call = FALSE
            )
        },
        error = function(e) {
            saveRDS(
                list(
                    dv = dotyobs,
                    ivs = dotxobs,
                    fam = "baseline",
                    thrs = thrs,
                    nthrs = nthrs,
                    npcs_range = npcs_range,
                    K = K,
                    fit_measure = fit_measure,
                    min_features = 1,
                    e = e$message
                ),
                file = paste0(
                    "./",
                    format(Sys.time(), "%Y%m%d-%H%M%S"),
                    "-mice-call-gspcr-polyreg-error-cv.rds"
                )
            )
        }
    )

    # Estimate model with tuned parameters
    tryCatch(
        expr = {
            # Train model to tune parameters
            gspcr_est <- gspcr::est_gspcr(
                dv = dotyobs,
                ivs = dotxobs,
                fam = "baseline",
                active_set = gscpr_fit$solution$standard$active_set,
                ndim = gscpr_fit$solution$standard$Q
            )
        },
        error = function(e) {
            saveRDS(
                list(
                    dv = dotyobs,
                    ivs = dotxobs,
                    fam = "baseline",
                    thrs = thrs,
                    nthrs = nthrs,
                    npcs_range = npcs_range,
                    K = K,
                    fit_measure = fit_measure,
                    min_features = 1,
                    gscpr_fit = gscpr_fit,
                    e = e$message
                ),
                file = paste0(
                    "./",
                    format(Sys.time(), "%Y%m%d-%H%M%S"),
                    "-mice-call-gspcr-polyreg-error-est.rds"
                )
            )
        }
    )

    # Obtain imputations -------------------------------------------------------

    # Obtain predictions
    tryCatch(
        expr = {
            post <- predict(
                object = gspcr_est,
                newdata = x[!ry, ]
            )
        },
        error = function(e) {
            saveRDS(
                list(
                    dv = dotyobs,
                    ivs = dotxobs,
                    fam = "binomial",
                    thrs = thrs,
                    nthrs = nthrs,
                    npcs_range = npcs_range,
                    K = K,
                    fit_measure = fit_measure,
                    min_features = 1,
                    gscpr_fit = gscpr_fit,
                    x = x,
                    ry = ry,
                    gspcr_est = gspcr_est,
                    e = e$message
                ),
                file = paste0(
                    "./",
                    format(Sys.time(), "%Y%m%d-%H%M%S"),
                    "-mice-call-gspcr-polyref-error-pred.rds"
                )
            )
        }
    )

    # Reshape as matrix if there was a single value
    if (sum(wy) == 1) {
        post <- matrix(post, nrow = 1, ncol = length(post))
    }

    # Count the number of unique values in y
    nc <- nlevels(dotyobs)

    # Sample from a uniform distribution as many values as unique values in y
    un <- rep(runif(sum(wy)), each = nc)

    # If post happens to be a vector
    if (is.vector(post)) {
        post <- matrix(c(1 - post, post), ncol = 2)
    }

    # Define the draws based on the predicted and sampled probabilities
    draws <- un > apply(post, 1, cumsum)

    # Reshape draws as manageable information
    idx <- 1 + apply(draws, 2, sum)

    # Revert to characters
    imputes <- levels(dotyobs)[idx]

    # Return
    imputes
}