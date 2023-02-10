#' Imputation by generalized supervised principal component regression
#'
#' Imputes univariate missing data using generalized supervised principal component regression.
#'
#' @aliases mice.impute.gspcr gspcr
#' @inheritParams mice.impute.norm.boot
#' @param maxnpcs The maximum number of PCs allowed.
#' @param K The number of folds for the cross-validation of the association-threshold.
#' @param nthrs Number of threshold values to be used for cross-validation.
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
mice.impute.gspcr <- function(y, ry, x, wy = NULL,
                              # theta = seq(0.05, .95, by = .05),
                              nthrs = 10,
                              maxnpcs = 3,
                              K = 5,
                              ...) {
  # Set up ---------------------------------------------------------------------

  if (is.null(wy)) wy <- !ry

  # Take bootstrap sample for model uncertainty
  n1 <- sum(ry)
  s <- sample(n1, n1, replace = TRUE)
  dotxobs <- x[ry, , drop = FALSE][s, ]
  dotyobs <- y[ry][s]

  # Cross-validate threshold ---------------------------------------------------

  # Family model
  fam <- "gaussian"

  # Fit null model
  glm0 <- glm(dotyobs ~ 1, family = fam)

  # Fit univariate models
  glm.fits <- lapply(1:ncol(dotxobs), function(j) {
      glm(dotyobs ~ dotxobs[, j], family = fam)
  })

  # Extract Log-likelihood values
  ll0 <- as.numeric(logLik(glm0))
  lls <- sapply(glm.fits, function(m) as.numeric(logLik(m)))

  # Define log-likelihood threshold values
  ll_thrs <- seq(from = min(lls), to = max(lls), length.out = nthrs)

  # Create a map of active predictors for threshold values
  map_active <- sapply(1:nthrs, function(a) lls > ll_thrs[a])

  # Create an object to store k-fold cross-validation log-likelihoods
  map_kfll <- array(dim = c(maxnpcs, nthrs, K))

  # Create a fold partitioning object
  part <- sample(rep(1:K, ceiling(nrow(dotxobs) / K)))[1:nrow(dotxobs)]

  # Compute cross-validation measures
  # Loop over K folds
  for (k in 1:K) {
    # Create fold data:
    Xtr <- dotxobs[part != k, , drop = FALSE]
    Xva <- dotxobs[part == k, , drop = FALSE]
    ytr <- dotyobs[part != k]
    yva <- dotyobs[part == k]

    # Loop over threshold values
    for (thr in 1:nthrs) {
      # thr <- 1
      # Define the active set of predictors based on the current threshold value
      aset <- map_active[, thr]

      # If there is more than 1 active variable
      if (sum(aset) > 1) {
        # Scale Xs
        Xtr_thr <- scale(Xtr[, aset], center = TRUE, scale = TRUE)
        Xva_thr <- scale(Xva[, aset],
          center = attributes(Xtr_thr)$`scaled:center`,
          scale = attributes(Xtr_thr)$`scaled:scale`
        )

        # Perform PCA
        svd_Xtr <- svd(Xtr_thr)

        # Project validation data on the PCs
        PCsva <- Xva_thr %*% svd_Xtr$v # TODO: check project is done simply with V

        # Check how many components are available (effective number)
        q.eff <- min(sum(aset), maxnpcs)

        # if (sum(aset) >= maxnpcs) {
          # q.eff <- maxnpcs

          # Select the PC scores that are available
          PCsva.eff <- PCsva[, 1:q.eff, drop = FALSE]

          # Compute the F-statistic for the possible additive PCRs
          for (Q in 1:q.eff) {
            # For the k-th PCs
            # Q <- 1

            # Estimate GLM
            glm.fit <- glm(yva ~ PCsva.eff[, 1:Q], family = fam)

            # Store the F statistic (used as a scaled value of the chi-square stat)
            map_kfll[Q, thr, k] <- as.numeric(logLik(glm.fit))
          }
        # }
      }
    }
  }

  # K-fold cross-validation log-likelihood scores
  map_kfll_avg <- apply(map_kfll, c(1, 2), mean, na.rm = FALSE)

  # Put on likelihood scale
  map_kfll_avg_l <- exp(map_kfll_avg)

  # K-fold Cross-Validation Choice
  coord_KCVC <- which(
    map_kfll_avg_l == max(map_kfll_avg_l, na.rm = TRUE), #TODO: min or max?
    arr.ind = TRUE
  )

  # Which threshold has been selected?
  thr.cv <- coord_KCVC[, "col"]

  # How many npcs have been selected
  Q.cv <- coord_KCVC[, "row"]
 
  # Train PCR on dotxobs sample
  pcr_out <- pls::pcr(
    dotyobs ~ dotxobs[, map_active[, thr.cv], drop = FALSE],
    ncomp = Q.cv,
    scale = TRUE,
    center = TRUE,
    validation = "none"
  )

  # Compute the residual sum of squares
  res_ss <- sum(stats::resid(pcr_out)[, , Q.cv]^2)

  # Compute degrees of freedom
  res_df <- n1 - Q.cv

  # Compute sigma
  sigma <- sqrt(res_ss / res_df)

  # Get prediction on (active) missing part
  yhat <- predict(
    object = pcr_out,
    newdata = x[wy, map_active[, thr.cv], drop = FALSE],
    ncomp = Q.cv,
    type = "response"
  )

  # Add noise for imputation uncertainty
  imputes <- yhat + rnorm(sum(wy)) * sigma

  # Return
  return(imputes)
}
