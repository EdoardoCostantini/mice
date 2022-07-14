#' Imputation by supervised principal component regression
#'
#' Imputes univariate missing data using supervised principal component regression.
#'
#' @aliases mice.impute.spcr spcr
#' @inheritParams mice.impute.norm.boot
#' @param npcs The number of principal components to keep or the automatic method
#' to use for decision.
#' @param nfolds The number of folds for the cross-validation of the association-threshold.
#' @param thresholds Vector of possible values for the association-threshold on
#' which cross-validation will be performed.
#' @return Vector with imputed data, the same type as \code{y}, and of length
#' \code{sum(wy)}
#' @details
#' Imputation of \code{y} by supervised principal component regression (Bair, 2006).
#' The method consists of the following steps:
#' \enumerate{
#' \item All of the potential predictors returning an R-squared larger than a 
#' threshold \eqn{\theta} are selected as an active set of predictors.
#' \item For a given \code{y} variable under imputation, draw a bootstrap version y*
#' with replacement from the observed cases \code{y[ry]}, and stores in x* the
#' corresponding values from \code{x[ry, ]}.
#' \item Regress \code{y*} on every potential predictor (simple linear 
#' regression), and store their R-squared.
#' \item Fit a PC regression with \code{y*} as the outcome, the predictors in 
#' \code{x*} whose simple regression R-square exceeds a defined threshold,
#' and \code{npcs} components.
#' \item Calculate the estimated residual standard error \code{sigma} based on the residuals
#' obtained from the PC regression and \code{npcs - 1} degrees of freedom.
#' \item Obtain predicted values for \code{y} based on the fitted PC regression
#' and the new data \code{x[wy, ]}
#' \item Obtain imputations by adding noise scaled by \code{sigma} to these
#' predictions.
#' }
#'
#' K-fold cross-validation is used to select a threshold value among a user-defined
#' grid of values.
#' For every provided threshold in the range [0, 1], all predictors with an R-square larger
#' than the threshold form an active set of predictors.
#' Then, \code{npcs} PCs are extracted from each active set and used to predict 
#' the dependent variable.
#' The active set giving best K-fold validation MSE is kept.
#'
#' This function allows the specification of a custom value for `npcs`.
#' In a dataset where few predictors are associated with the variables under imputation,
#' the requested number of PCs may be larger than the number of predictors
#' selected for a given cross-validated threshold.
#' If this is the case, the maximum number of PCs supported is used.
#' This means that the predictors selected based on the cross-validation threshold
#' are projected on a new space where they are independent, but no dimensionality
#' reduction is performed.
#' Similarly, it may happen that for a given threshold value, fewer predictors 
#' are kept than the number of PCs requested by the user.
#' In this case, the maximum number of npcs supported is used in the cross-validation procedure.
#'
#' The user can specify a \code{predictorMatrix} in the \code{mice} call
#' to define which predictors are provided to this univariate imputation method.
#' Therefore, users may force the exclusion of a predictor from a given
#' imputation model by speficing a \code{0} entry.
#' However, a non-zero entry does not guarantee the variable will be used,
#' as this decision is ultimately made based on the k-fold cross-validation
#' procedure.
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
mice.impute.spcr <- function(y, ry, x, wy = NULL,
                             thresholds = seq(0.1, .9, by = .1),
                             npcs = 1, 
                             nfolds = 10,
                             ...) {
  # Set up ---------------------------------------------------------------------

  if (is.null(wy)) wy <- !ry

  # Take bootstrap sample for model uncertainty
  n1 <- sum(ry)
  s <- sample(n1, n1, replace = TRUE)
  dotxobs <- x[ry, , drop = FALSE][s, ]
  dotyobs <- y[ry][s]

  # Cross-validate threshold ---------------------------------------------------

  # Obtain R-squared for all simple linear regression models
  r2_vec <- apply(x, 2, function(j) {
    sqrt(summary(lm(y ~ j))$r.squared)
  })

  # DEfine predictor groups (pred groups) based on different thresholds
  pred_groups <- lapply(thresholds, function(m) {
    preds <- colnames(x)[r2_vec >= m]
    if (length(preds) >= 1) {
      preds
    } else {
      NULL
    }
  })

  # If thresholds used lead only to empty pred groups, say so
  if(all(sapply(pred_groups, is.null)) == TRUE){
    stop(
      paste0(
        "The threshold values used are too high. Try using a lower range."
      )
    )
  }
  
  # Drop empty pred_groups slots
  pred_groups <- pred_groups[!sapply(pred_groups, is.null)]

  # Drop possible duplicated pred_groups slots
  pred_groups <- unique(pred_groups)

  # Drop preds groups that are smaller than required npcs
  pred_groups <- pred_groups[sapply(pred_groups, length) >= npcs]

  # If there is no pred group with enough predictors for the required npcs, say so
  if(length(pred_groups) == 0){
    stop(
      paste0(
        "There is no threshold value that can select enough predictors to extract ",
        npcs, " PCs. Try using a smaller npcs or lower thresholds."
      )
    )
  }

  # Create a partition vector
  part <- sample(rep(1:nfolds, ceiling(nrow(dotxobs) / nfolds)))[1:nrow(dotxobs)]

  # Obtain Cross-validation error
  cve_obj <- lapply(pred_groups, function(set) {
    .spcrCVE(
      dv = dotyobs,
      pred = dotxobs[, set, drop = FALSE],
      K = nfolds,
      part = part,
      npcs = npcs
    )
  })

  # Extract CVEs
  cve <- sapply(cve_obj, "[[", 1)
  preds_active <- pred_groups[[which.min(cve)]]

  # Train PCR on dotxobs sample
  pcr_out <- pls::pcr(
    dotyobs ~ dotxobs[, preds_active, drop = FALSE],
    ncomp = npcs,
    scale = TRUE,
    center = TRUE,
    validation = "none"
  )

  # Define sigma
  RSS <- sqrt(sum(pcr_out$residuals^2))
  sigma <- RSS / (n1 - npcs - 1)

  # Get prediction on (active) missing part
  yhat <- predict(
    object = pcr_out,
    newdata = x[wy, preds_active, drop = FALSE],
    ncomp = npcs,
    type = "response"
  )

  # Add noise for imputation uncertainty
  imputes <- yhat + rnorm(sum(wy)) * sigma

  # Return
  return(imputes)
}

#' Computes the cross-validation error
#'
#' @aliases .spcrCVE
#' @param dv Vector of observed values on the dv.
#' @param pred Matrix of complete covariates for a version of a model.
#' @param K Unit vector (double) indicating the number of folds for cross-validation.
#' @param part Vector indicating how to partition data in K-folds.
#' @param npcs Unit vector (double) indicating the number of components to retain.
#' @export
.spcrCVE <- function(dv, pred, part, K = 10, npcs = 1) {
  # Input examples
  # dv   = mtcars[, 1]
  # pred = mtcars[, -1]
  # K    = 10
  # npcs = 5
  # part = sample(rep(1 : K, ceiling(nrow(mtcars) / K)))[1 : nrow(mtcars)]

  # Install packages on demand for this function
  install.on.demand("MLmetrics")

  # Define a safe number of pcs
  q <- min(npcs, ncol(pred))

  # Create an empty storing object
  mse <- rep(NA, K)

  # Loop over K folds
  for (k in 1:K) {

    # Partition data:
    Xtr <- pred[part != k, , drop = FALSE]
    Xva <- pred[part == k, , drop = FALSE]
    ytr <- dv[part != k]
    yva <- dv[part == k]

    # Calibrate PCR on training datest
    pcr_out <- pls::pcr(
      ytr ~ Xtr,
      ncomp = q,
      scale = TRUE,
      center = TRUE,
      validation = "none"
    )

    # Get prediction on validation data set
    yva_hat <- predict(pcr_out, newdata = Xva, ncomp = q, type = "response")

    # Save MSE
    mse[k] <- MLmetrics::MSE(
      y_pred = yva_hat,
      y_true = yva
    )

  }

  # Return the CVE:
  cve <- sum(mse * (table(part) / length(part)))

  # Return
  return(list(
    cve = cve,
    npcs = q
  ))

}