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
#' The method consists of the following steps:
#' \enumerate{
#' \item The observed part of the variable under imputation is regressed on every potential predictor (simple linear regression), and the square root of the R-squared is stored.
#' \item All of the potential predictors returning an R-squared larger than a threshold \theta are selected as an active set of predictors.
#' \item `npcs` principal components are extracted from this active set.
#' \item These principal components are used as input for a `norm.boot` univariate imputation algorithm
#' }
#' 
#' K-fold cross-validation is used to select a threshold value among a user-defined
#' vector of possible values.
#' For every given value in the range [0, 1], all predictors with an R-square larger 
#' than the threshold form an active set of predictors.
#' Then, `npcs` PCs are extracted from each active set and used to predict the dependent variable.
#' The active set giving best validation MSE for `npcs` PCs is kept.
#' 
#' This function allows the specification of a custom value for `npcs`.
#' In a dataset where few predictors are associated with the variables under imputation,
#' the requested number of PCs may be larger than the number of predictors 
#' selected for a given cross-validated threshold.
#' If this is the case, the maximum number of PCs supported is used.
#' This means that the predictors selected based on the cross-validation threshold
#' are projected on a new space where they are independent but no dimensionality 
#' reduction is performed.
#' 
#' The user can specify a \code{predictorMatrix} in the \code{mice} call
#' to define which predictors are provided to this univariate imputation method.
#' Therefore, users may force the exclusion of a predictor from a given
#' imputation model by speficing a \code{0} entry.
#' However, a non-zero entry does not guarantee the variable will be used,
#' as this decision is ultimately made based on the k-fold cross-validation
#' procedure.
#'
#' The method is based on the supervised principal component prediction approach proposed
#' by Bair et. al. (2006).
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
                             thresholds = seq(.1, .9, by = .1),
                             npcs = 1, nfolds = 10,
                             ...) {
  # Set up ---------------------------------------------------------------------

  # Body
  if (is.null(wy)) wy <- !ry

  # Cross-validate threshold ---------------------------------------------------
  # Fit univariate regression models (efficient correlation coefficient)
  r2_vecs <- apply(x, 2, function (j){
    sqrt(summary(lm(y ~ j))$r.squared)
  })

  # Cross-validate threshold ---------------------------------------------------

  # Predictors based on different thrasholds
  mods  <- lapply(thresholds, function (r2) {
    mod <- x[, r2_vecs >= r2, drop = FALSE]
    if(ncol(mod) >= 1){
      mod
    } else {
      NULL
    }
  })
  # TODO: what if thresholds are too high and nothing is selected?

  # Drop empty mods
  mods <- mods[!sapply(mods, is.null)]

  # Create a partition vector:
  part <- sample(rep(1 : nfolds, ceiling(nrow(x) / nfolds)))[1 : nrow(x)]

  # Obtain Cross-validation error
  cve <- sapply(mods, .spcrCVE, dv = y, K = nfolds, part = part)

  # Select predictors giving model with smallest error
  x_sub <- mods[[which.min(cve)]]

  # Extract PCs from the predictors involved in this model
  pcr_out <- stats::prcomp(x_sub,
                           center = TRUE,
                           scale = TRUE)

  # Compute Explained Variance by each principal component
  pc_var_exp <- prop.table(pcr_out$sdev^2)

  # Keep PCs based on npcs object
  x_pcs <- pcr_out$x[, 1:npcs, drop = FALSE]
  pca_exp <- sum(pc_var_exp[1:npcs])

  # Impute ---------------------------------------------------------------------

  # Use traditional norm.boot machinery to obtain prediction
  imputes <- mice.impute.norm.boot(y = y, ry = ry, x = x_pcs, wy = NULL)

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
.spcrCVE <- function(dv, pred, K, part, npcs = 1) {
    # Input examples
    # dv   = mtcars[, 1]
    # pred = mtcars[, -1]
    # K    = 10
    # part = sample(rep(1 : K, ceiling(nrow(mtcars) / K)))[1 : nrow(mtcars)]

    # Install packages on demand for this function
    install.on.demand("MLmetrics")

    # Extract PCs from the predictors involved in this model
    pcr_out <- stats::prcomp(pred,
                             center = TRUE,
                             scale = TRUE)

    # Extract PCs from this set of predictors
    x_pcs <- pcr_out$x[, 1:npcs, drop = FALSE]

    # Put them together
    data <- data.frame(y = dv, x_pcs)

    # Create model
    model <- paste0("y ~ ", paste0(colnames(x_pcs), collapse = " + "))

    # Create empty storing object
    mse <- rep(NA, K)

    # Loop over K repititions:
    for(k in 1 : K) {

        # Partition data:
        train <- data[part != k, ]
        valid <- data[part == k, ]

        # Fit model, generate predictions, and save the MSE:
        fit    <- lm(model, data = train)

        # Generate predictions
        dv_hat   <- predict(fit, newdata = valid)

        # Save MSE
        mse[k] <- MLmetrics::MSE(y_pred = dv_hat,
                                 y_true = dv[part == k])

    }

    # Return the CVE:
    sum((table(part) / length(part)) * mse)
}