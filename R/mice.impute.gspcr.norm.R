#' Imputation by generalized supervised principal component regression
#'
#' Imputes univariate missing data using generalized supervised principal component regression.
#'
#' @aliases mice.impute.gspcr gspcr
#' @inheritParams mice.impute.norm.boot
#' @param thrs character vector of length 1 storing the type of threshold to be used (see below for available options)
#' @param nthrs numeric vector of length 1 storing the number of threshold values to be used
#' @param fit_measure character vector of length 1 indicating the type of fit measure to be used in the cross-validation procedure
#' @param npcs_range numeric vector defining the numbers of principal components to be used
#' @param K numeric vector of length 1 storing the number of folds for the K-fold cross-validation procedure
#' @return Vector with imputed data, the same type as \code{y}, and of length
#' \code{sum(wy)}
#' @details
#' Imputation of \code{y} by supervised principal component regression (Bair, 2006).
#' The method consists of the following steps:
#' \enumerate{
#' \item For a given \code{y} variable under imputation, draw a bootstrap version y*
#' with replacement from the observed cases \code{y[ry]}, and stores in x* the
#' corresponding values from \code{x[ry, ]}.
#' \item Compute bivariate association measure \eqn{\theta} between the observed part of the variable under and every potential predictor in the data.
#' \item Collect all the predictors with \eqn{\theta} higher than some threshold \eqn{\theta_t} and define them as the active set of predictors.
#' \item Regress \code{y*} on the Principal components computed on the active set.
#' \item Calculate the estimated residual standard error \code{sigma} based on the residuals obtained from the PC regression and \code{npcs - 1} degrees of freedom.
#' \item Obtain predicted values for \code{y} based on the fitted PC regression
#' and the new data \code{x[wy, ]}
#' \item Obtain imputations by adding noise scaled by \code{sigma} to these
#' predictions.
#' }
#'
#' The user specifies a number of association values to be checked as threshold values and a range of number of components to be checked.
#' Then, the range between the minimum and maximum value the bivairate association measure between y and the x is divided into equally spaced steps to produce the desired number of association measures.
#' Every association measure value is used to define a different active set.
#' For every active set, all requested PCs are computed and used to predict the dependent variable. If the sets npcs to 1, 3, and 5, then the dependent variable is regressed on 1 component, then on 1 to 3 components and finally on 1 to 5 components.
#' The combination of active set and npcs that returns the optimal value of a given fit measure is then selected.
#' This selection of the threshold value \eqn{\theta_t} and the number of components Q can be done through K-fold cross-validation by setting the value of the folds to anything larger than 1.
#' 
#' For details on the admissible values for \code{thrs}, \code{npcs_range}, and \code{fit_measure} consult the help file for the underlying gspcr engine [gspcr::cv_gspcr()].
#'
#' The user can specify a \code{predictorMatrix} in the \code{mice} call
#' to define which predictors are provided to this univariate imputation method.
#' Therefore, users may force the exclusion of a predictor from a given
#' imputation model by specifying a \code{0} entry.
#' However, a non-zero entry does not guarantee the variable will be used,
#' as this decision is ultimately made based on the k-fold cross-validation
#' procedure.
#' 
#' When using any gspcr-based method, it is recommended to set the following arguments in the mice call
#' - \code{eps = 0}, to bypasses \code{remove.lindep()}
#' - \code{threshold = 1L}, to bypass collinearity checks
#' These aspects are already handled by gspcr.
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
mice.impute.gspcr.norm <- function(y, ry, x, wy = NULL,
                                   thrs = "PR2",
                                   fit_measure = "BIC",
                                   nthrs = 10,
                                   npcs_range = 1:5,
                                   K = 1,
                                   ...) {
  # Set up ---------------------------------------------------------------------

  if (is.null(wy)) wy <- !ry

  # Bootstrap sample for model uncertainty -------------------------------------

  # Sample size of responses
  n1 <- sum(ry)

  # Define bootstrap sample
  s <- sample(n1, n1, replace = TRUE)

  # Create bootstrap sample observed predictors
  dotxobs <- x[ry, , drop = FALSE][s, ]

  # Create bootstrap sample of observed values of variable under imputation
  dotyobs <- y[ry][s]

  # Train GSPCR ----------------------------------------------------------------

  tryCatch(
    expr = {
      # Train model to tune parameters
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
    },
    error = function(e) {
      saveRDS(
        list(
          dv = dotyobs,
          ivs = dotxobs,
          fam = "gaussian",
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
          "-mice-call-gspcr-norm-error-cv.rds"
        )
      )
    }
  )

  # Estimate GSPCR -------------------------------------------------------------

  tryCatch(
    expr = {
      # Train model to tune parameters
      gspcr_est <- gspcr::est_gspcr(gscpr_fit)
    },
    error = function(e) {
      saveRDS(
        list(
          dv = dotyobs,
          ivs = dotxobs,
          fam = "gaussian",
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
          "-mice-call-gspcr-norm-error-est.rds"
        )
      )
    }
  )

  # Obtain imputations ---------------------------------------------------------

  # Obtain predictions
  tryCatch(
    expr = {
      y_hat <- predict(object = gspcr_est)
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
          gspcr_est = gspcr_est,
          e = e$message
        ),
        file = paste0(
          "./",
          format(Sys.time(), "%Y%m%d-%H%M%S"),
          "-mice-call-gspcr-norm-error-pred.rds"
        )
      )
    }
  )

  # Compute residuals
  dotyobs_res <- dotyobs - y_hat

  # Compute error variance
  sigma <- sqrt(sum(dotyobs_res^2) / (n1 - length(coef(gspcr_est$glm_fit))))

  # Obtain predictions
  y_hat <- predict(
    object = gspcr_est,
    newdata = x[!ry, ]
  )

  # Add noise for imputation uncertainty
  imputes <- y_hat + rnorm(sum(wy)) * sigma

  # Return imputations
  return(imputes)
}
