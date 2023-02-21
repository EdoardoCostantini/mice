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

  # Cross-validate the npcs and threshold
  spcr.cv <- .spcr.cv(
    dv = dotyobs,
    ivs = dotxobs,
    fam = "gaussian",
    nthrs = 10
  )

  # Train PCR on dotxobs sample
  pcr_out <- pls::pcr(
    dotyobs ~ dotxobs[, spcr.cv$pred.active, drop = FALSE],
    ncomp = spcr.cv$Q.cv,
    scale = TRUE,
    center = TRUE,
    validation = "none"
  )
  # TODO: make this your function

  # Compute the residual sum of squares
  res_ss <- sum(stats::resid(pcr_out)[, , Q.cv]^2)

  # Compute degrees of freedom
  res_df <- n1 - Q.cv

  # Compute sigma
  sigma <- sqrt(res_ss / res_df)

  # Get prediction on (active) missing part
  yhat <- predict(
    object = pcr_out,
    newdata = x[wy, pred.map[, thr.cv], drop = FALSE],
    ncomp = Q.cv,
    type = "response"
  )

  # Add noise for imputation uncertainty
  imputes <- yhat + rnorm(sum(wy)) * sigma

  # Return
  return(imputes)
}

# Cross-validate function

.spcr.cv <- function(
  dv, 
  ivs, 
  fam = "gaussian",
  thrs = c("LLS", "pseudoR2", "normalized")[1],
  nthrs = 10,
  maxnpcs = 3,
  K = 5,
  test = c("LRT", "F", "MSE")[2],
  max.features = ncol(ivs),
  min.features = 5
  ) {

  # Example inputs
  # dv <- mtcars[, 1]
  # ivs <- mtcars[, -1]
  # thrs = c("LLS", "pseudoR2")[1]
  # nthrs = 10
  # fam <- c("gaussian", "binomial", "poisson")[1]
  # maxnpcs <- 5
  # K = 5
  # test = c("LRT", "F", "MSE")[3]
  # max.features = ncol(ivs)
  # min.features = 1

  # Sample size
  n <- nrow(ivs)

  # Fit null model
  glm0 <- glm(dv ~ 1, family = fam)

  # Fit univariate models
  glm.fits <- lapply(1:ncol(ivs), function(j) {
    glm(dv ~ ivs[, j], family = fam)
  })

  # Extract Log-likelihood values
  ll0 <- as.numeric(logLik(glm0))
  lls <- sapply(glm.fits, function(m) as.numeric(logLik(m)))

  # Create active sets based on threshold type

  if(thrs == "LLS"){

    # Use the logLikelihoods as bivariate association scores
    ascores <- lls

    # Give it good names
    names(ascores) <- colnames(ivs)

    # Define the upper and lower bounds of the association
    lower <- min(ascores)
    upper <- max(ascores)

  }

  if(thrs == "pseudoR2"){

    # Compute pseudo R-squared
    CNR2 <- 1 - exp(-2 / n * (lls - ll0))

    # Give it good names
    names(CNR2) <- colnames(ivs)

    # Make them correlation coefficients
    ascores <- sqrt(CNR2)

    # Define upper and lower bounds of the association
    lower <- quantile(ascores, 1 - (max.features / ncol(ivs)))
    upper <- quantile(ascores, 1 - (min.features / ncol(ivs)))

  }

  if (thrs == "normalized") {
    
    # Set objects to the required dimension
    x <- t(as.matrix(ivs))
    y <- dv
    featurenames <- colnames(ivs)

    # Empty
    s0.perc <- NULL

    # Sample size
    n <- length(y)

    # Compute vector of feature means
    xbar <- x %*% rep(1 / n, n)

    # Same as computing the row means
    cbind(xbar, rowMeans(x))

    # Compute the diagonal of the cross-product matrix between variables
    sxx <- ((x - as.vector(xbar))^2) %*% rep(1, n)

    # Which is the mid step for variance
    cbind(sxx, apply(x - as.vector(xbar), 1, var) * (n - 1))

    # Compute the cross-product matrix between X and Y
    sxy <- (x - as.vector(xbar)) %*% (y - mean(y))

    # Which is the mid step for covariance between the two
    cbind(sxx, apply(x - as.vector(xbar), 1, var) * (n - 1))

    # Total sum of squares
    syy <- sum((y - mean(y))^2)

    # Ratio of the two
    numer <- sxy / sxx

    # Compute sd?
    sd <- sqrt((syy / sxx - numer^2) / (n - 2))

    # add "fudge"(?) to the denominator
    if (is.null(s0.perc)) {
      fudge <- median(sd)
    }
    if (!is.null(s0.perc)) {
      if (s0.perc >= 0) {
        fudge <- quantile(sd, s0.perc)
      }
      if (s0.perc < 0) {
        fudge <- 0
      }
    }

    # Ratio between numerator and sd
    tt <- numer / (sd + fudge)

    # Store the normalized correlation scores
    ascores <- abs(tt)[, 1]

    # Define upper and lower bounds of the normalized correlation
    lower <- quantile(abs(ascores), 1 - (max.features / nrow(x)))
    upper <- quantile(abs(ascores), 1 - (min.features / nrow(x)))

  }

  # Define threshold values
  thrs_values <- seq(from = lower, to = upper, length.out = nthrs)

  # Create a map of active predictors based on threshold values
  pred.map <- sapply(1:nthrs, function(a) ascores > thrs_values[a])

  # Use thresholds as names
  colnames(pred.map) <- round(thrs_values, 3)

  # If two thresholds are giving the same result reduce the burden
  pred.map <- pred.map[, !duplicated(t(pred.map))]

  # Get rid of thresholds that are keeping too few predictors
  pred.map <- pred.map[, colSums(pred.map) >= min.features]

  # Get rid of thresholds that are keeping too many predictors
  pred.map <- pred.map[, colSums(pred.map) <= max.features]

  # And update the effective number of the thresholds considered
  nthrs.eff <- ncol(pred.map)
  
  # Create an object to store k-fold cross-validation log-likelihoods
  map_kfcv <- array(
    dim = c(maxnpcs, nthrs.eff, K),
    dimnames = list(NULL, colnames(pred.map), NULL)
  )

  # Create a fold partitioning object
  part <- sample(rep(1:K, ceiling(nrow(ivs) / K)))[1:nrow(ivs)]

  # Loop over K folds
  for (k in 1:K) {
    # k <- 5

    # Create fold data:
    Xtr <- ivs[part != k, , drop = FALSE]
    Xva <- ivs[part == k, , drop = FALSE]
    ytr <- dv[part != k]
    yva <- dv[part == k]

    # Null model
    glm.fit0 <- glm(yva ~ 1, family = fam)

    # Loop over threshold values
    for (thr in 1:nthrs.eff) {
      # thr <- 1
      # Define the active set of predictors based on the current threshold value
      aset <- pred.map[, thr]

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
        PCsva <- Xva_thr %*% svd_Xtr$v

        # Check how many components are available (effective number)
        q.eff <- min(sum(aset), maxnpcs)

        # Select the PC scores that are available
        PCsva.eff <- PCsva[, 1:q.eff, drop = FALSE]

        # Compute the F-statistic for the possible additive PCRs
        for (Q in 1:q.eff) {
          # Q <- 1

          # Estimate GLM models
          glm.fit <- glm(yva ~ PCsva.eff[, 1:Q], family = fam)

          # Extract desired statistic
          if (test == "LRT") {
            map_kfcv[Q, thr, k] <- as.numeric(- 2 * (logLik(glm.fit0) - logLik(glm.fit)))
          }
          if (test == "F") {
            map_kfcv[Q, thr, k] <- anova(glm.fit0, glm.fit, test = "F")$F[2]
          }
          if (test == "PR2") {
            map_kfcv[Q, thr, k] <- as.numeric(1 - exp(-2 / n * (logLik(glm.fit) - logLik(glm.fit0))))
          }
          if (test == "MSE") {
            map_kfcv[Q, thr, k] <- MLmetrics::MSE(y_pred = predict(glm.fit), y_true = yva)
          }
          if (test == "BIC") {
            map_kfcv[Q, thr, k] <- as.numeric(log(n) * (Q + 1 + 1) - 2 * logLik(glm.fit))
          }
        }
      }

    }
  }

  # Average selected score across folds

  if(test == "F"){
    # average F scores on a more symmetrical scale
    lscor <- apply(log(map_kfcv), c(1, 2), mean, na.rm = FALSE)

    # revert to correct scale
    scor <- exp(lscor)

    # K-fold Cross-Validation solution
    kfcv_sol <- which(
      scor == max(scor, na.rm = TRUE), # TODO: min or max?
      arr.ind = TRUE
    )
  }

  if(test == "LRT" | test == "PR2"){
    # Mean of the likelihood ratio test statistics
    scor <- apply(map_kfcv, c(1, 2), mean, na.rm = FALSE)

    # K-fold Cross-Validation solution
    kfcv_sol <- which(
      scor == max(scor, na.rm = TRUE), # TODO: min or max?
      arr.ind = TRUE
    )

  }

  if (test == "MSE" | test == "BIC") {
    # Mean of the likelihood ratio test statistics
    scor <- apply(map_kfcv, c(1, 2), mean, na.rm = FALSE)

    # K-fold Cross-Validation solution
    kfcv_sol <- which(
      scor == min(scor, na.rm = TRUE),
      arr.ind = TRUE
    )
  }

  # Which threshold has been selected?
  thr.cv <- as.numeric(names(scor[kfcv_sol[1], kfcv_sol[2]]))

  # How many npcs have been selected?
  Q.cv <- as.numeric(kfcv_sol[, "row"])

  # Return
  return(
    list(
      thr.cv = thr.cv,
      thr = thrs_values,
      Q.cv = Q.cv,
      scor = scor,
      pred.map = pred.map,
      pred.active = rownames(pred.map)[pred.map[, kfcv_sol[, "col"]]]
    )
  )
}

.spcr.cv.pkg <- function(dv, ivs, fam = "gaussian", nthrs = 10, maxnpcs = 3, K = 5) {
  # Example inputs
  # dv <- mtcars[, 1]
  # ivs <- t(mtcars[, -1])
  # nthrs = 10
  # fam <- c("gaussian", "binomial", "poisson")[2]
  # maxnpcs <- 3
  # K = 5

  # Process data
  x = t(as.matrix(ivs))
  y = dv
  featurenames = colnames(ivs)

  #
  s0.perc <- NULL

  # Sample size
  n <- length(y)

  # Compute vector of feature means
  xbar <- x %*% rep(1 / n, n)

  # Same as computing the row means
  cbind(xbar, rowMeans(x))

  # Compute the diagonal of the cross-product matrix between variables
  sxx <- ((x - as.vector(xbar))^2) %*% rep(1, n)

  # Which is the mid step for variance
  cbind(sxx, apply(x - as.vector(xbar), 1, var) * (n - 1))

  # Compute the cross-product matrix between X and Y
  sxy <- (x - as.vector(xbar)) %*% (y - mean(y))

  # Which is the mid step for covariance between the two
  cbind(sxx, apply(x - as.vector(xbar), 1, var) * (n - 1))

  # Total sum of squares
  syy <- sum((y - mean(y))^2)

  # Ratio of the two
  numer <- sxy / sxx

  # Compute sd?
  sd <- sqrt((syy / sxx - numer^2) / (n - 2))

  # add "fudge"(?) to the denominator
  if (is.null(s0.perc)) {
    fudge <- median(sd)
  }
  if (!is.null(s0.perc)) {
    if (s0.perc >= 0) {
      fudge <- quantile(sd, s0.perc)
    }
    if (s0.perc < 0) {
      fudge <- 0
    }
  }

  # Ratio between numerator and sd
  tt <- numer / (sd + fudge)

  # Compute normalized correlation between y and every x
  feature.scores <- tt

  # Set up the same arguments
  fit <- train.obj
  data <- data.train
  n.threshold <- nthrs
  n.fold <- K
  folds <- NULL
  n.components <- maxnpcs
  min.features <- 5
  max.features <- nrow(data.train$x)
  compute.fullcv <- TRUE
  compute.preval <- TRUE
  xl.mode <- c(
    "regular",
    "firsttime",
    "onetime",
    "lasttime"
  )[1]
  xl.time <- NULL
  xl.prevfit <- NULL

  # Type of fit
  type <- fit$type

  # Number of components
  n.components <- min(5, n.components)

  # Sample size
  n <- ncol(data$x)

  # Store the normalized correlation scores
  cur.tt <- fit$feature.scores

  # Define upper and lower bounds of the normalized correlation
  lower <- quantile(abs(cur.tt), 1 - (max.features / nrow(data$x)))
  upper <- quantile(abs(cur.tt), 1 - (min.features / nrow(data$x)))

  # Number of folds
  folds <- vector("list", n.fold)
  breaks <- round(seq(from = 1, to = (n + 1), length = (n.fold + 1)))
  cv.order <- sample(1:n)
  for (j in 1:n.fold) {
    folds[[j]] <- cv.order[(breaks[j]):(breaks[j + 1] - 1)]
  }

  # Storing objects
  featurescores.folds <- matrix(nrow = nrow(data$x), ncol = n.fold)

  # Thresholds
  thresholds <- seq(from = lower, to = upper, length = n.threshold)

  # Check selected predictors
  map_active_pkg <- sapply(1:n.threshold, function(a) sort(cur.tt) > thresholds[a])

  # Other important objects
  nonzero <- rep(0, n.threshold)
  scor <- array(NA, c(n.components, n.threshold, n.fold))
  scor.preval <- matrix(NA, nrow = n.components, ncol = n.threshold)
  scor.lower <- NULL
  scor.upper <- NULL
  v.preval <- array(NA, c(n, n.components, n.threshold))

  # Define objects for the CV procedure
  first <- 1
  last <- n.fold

  # Then, let's go through the folds
  for (j in first:last) {
    # For the j-th fold
    # j <- 1

    # Create a temporary training data
    data.temp <- list(
      x = data$x[, -folds[[j]]],
      y = data$y[-folds[[j]]],
      censoring.status = data$censoring.status[-folds[[j]]]
    )

    # What are the normalized correlations for this fold
    cur.tt <- superpc.train(data.temp, type = type, s0.perc = fit$s0.perc)$feature.scores

    # Store that
    featurescores.folds[, j] <- cur.tt

    # For every threshold check
    for (i in 1:n.threshold) {
      # For the i-th threshold
      # i <- 1

      # Define active set (check which normalized correlations are smaller)
      cur.features <- (abs(cur.tt) > thresholds[i])

      # If there are more then 1 variables, we do this
      if (sum(cur.features) > 1) {
        # Store this number
        nonzero[i] <- nonzero[i] + sum(cur.features) / n.fold

        # Compute the SVD of the active set
        cur.svd <- mysvd(
          x = data$x[cur.features, -folds[[j]]],
          n.components = n.components
        )

        # Scale unseen data
        xtemp <- data$x[cur.features, folds[[j]], drop = FALSE]
        xtemp <- t(scale(t(xtemp),
          center = cur.svd$feature.means,
          scale = FALSE
        ))
        cur.v.all <- scale(t(xtemp) %*% cur.svd$u,
          center = FALSE,
          scale = cur.svd$d
        )

        # Check how many components are available (effective number)
        n.components.eff <- min(sum(cur.features), n.components)

        # Select the PC scores that are available
        cur.v <- cur.v.all[, 1:n.components.eff, drop = FALSE]

        # Store them for this threshold value
        v.preval[folds[[j]], 1:n.components.eff, i] <- cur.v

        # Compute the F-statistic for the possible additive PCRs
        for (k in 1:ncol(cur.v)) {
          # For the k-th PCs
          # k <- 2

          # Compute the linear model
          junk <- summary(lm(data$y[folds[[j]]] ~ cur.v[, 1:k]))

          # Store the F statistic (used as a scaled value of the chi-square stat)
          scor[k, i, j] <- junk$fstat[1]
        }
      }
    }
  }

  # Compute the log-likelihood scores?
  lscor <- apply(log(scor), c(1, 2), mean.na) # average on a more symmetrical scale
  se.lscor <- apply(log(scor), c(1, 2), se.na)
  scor.lower <- exp(lscor - se.lscor)
  scor.upper <- exp(lscor + se.lscor)
  scor <- exp(lscor)

  # Give names that make sense
  colnames(scor) <- round(thresholds, 3)

  # Compute the pre-validation
  if (compute.preval) {
    for (i in 1:n.threshold) {
      # i <- 1
      for (j in 1:n.components) {
        # j <- 1
        if (sum(is.na(v.preval[, 1:j, i])) == 0) {
          junk <- summary(lm(data$y ~ v.preval[, 1:j, i]))
          scor.preval[j, i] <- junk$fstat[1]
        }
      }
    }
  }

  # Store objects
  junk <- list(
    thresholds = thresholds,
    n.threshold = n.threshold,
    nonzero = nonzero,
    scor.preval = scor.preval,
    scor = scor,
    scor.lower = scor.lower,
    scor.upper = scor.upper,
    folds = folds,
    n.fold = n.fold,
    featurescores.folds = featurescores.folds,
    v.preval = v.preval,
    compute.fullcv = compute.fullcv,
    compute.preval = compute.preval,
    type = type,
    call = NULL
  )

  # K-fold Cross-Validation Choice
  coord_KCVC <- which(
    scor == max(scor, na.rm = TRUE), # TODO: min or max?
    arr.ind = TRUE
  )

  # Which threshold has been selected?
  thr.cv <- thresholds[coord_KCVC[, "col"]]

  # How many npcs have been selected?
  Q.cv <- coord_KCVC[, "row"]

  # Return
  return(
    list(
      thr.cv = thr.cv,
      # pred.active = pred.map[, coord_KCVC[, "col"]],
      thr = as.numeric(thresholds),
      Q.cv = Q.cv,
      scor = scor,
      # pred.map = pred.map,
      pred.active = rownames(cur.tt)[cur.tt > thr.cv]
    )
  )
}

# Load a special SVD function
mysvd <- function(x,
                  n.components = NULL) {
  # finds PCs of matrix x

  p <- nrow(x)
  n <- ncol(x)

  # center the observations (rows)
  feature.means <- rowMeans(x)
  x <- t(scale(t(x), center = feature.means, scale = FALSE))

  if (is.null(n.components)) {
    n.components <- min(n, p)
  }
  if (p > n) {
    a <- eigen(t(x) %*% x)
    v <- a$vec[, 1:n.components, drop = FALSE]
    d <- sqrt(a$val[1:n.components, drop = FALSE])
    u <- scale(x %*% v, center = FALSE, scale = d)

    return(list(
      u = u,
      d = d,
      v = v,
      feature.means = feature.means
    ))
  } else {
    junk <- svd(x)
    nc <- min(ncol(junk$u), n.components)

    return(list(
      u = junk$u[, 1:nc],
      d = junk$d[1:nc],
      v = junk$v[, 1:nc],
      feature.means = feature.means
    ))
  }
}

# Define two functions that compute the mean and the sd ignoring the NAs
mean.na <- function(x) {
  mean(x[!is.na(x)])
}

se.na <- function(x) {
  val <- NA
  if (sum(!is.na(x)) > 0) {
    val <- sqrt(var(x[!is.na(x)]) / sum(!is.na(x)))
  }
  return(val)
}

# Get test statistics given model of interest
getTestStat <- function(glm.fit, glm.fit0 = NULL, test = "LRT", y_true = NULL) {
  # glm.fit0 <- lm(mpg ~ 1, data = mtcars)
  # glm.fit <- lm(mpg ~ ., data = mtcars)
  if (test == "LRT") {
    test.stat <- as.numeric(logLik(glm.fit0)) - as.numeric(logLik(glm.fit))
  }
  if (test == "F") {
    test.stat <- anova(glm.fit, glm.fit0, test = "F")$F[2]
  }
  if (test == "MSE") {
    test.stat <- MLmetrics::MSE(y_pred = predict(glm.fit), y_true = y_true)
  }
  return(test.stat)
}
