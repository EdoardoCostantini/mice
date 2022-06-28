#' Imputation by principal covariates regression
#'
#' Imputes univariate missing data by using principal covariates regression (PCovR).
#'
#' @aliases mice.impute.pcovr pcovr
#' @inheritParams mice.impute.norm.boot
#' @param npcs The number of principal components to extract for PC regression.
#' @return Vector with imputed data, same type as \code{y}, and of length
#' \code{sum(wy)}
#' @details
#' Imputation of \code{y} by PCovR. The procedure is as follows:
#' \enumerate{
#' \item Draws a bootstrap sample from \code{x[ry,]} and \code{y[ry]}
#' \item Computes an optimal \code{alpha} value for the PCovR model on the basis 
#' of maximum likelihood principles as described in Vervloet et. al. (2015).
#' \item Estimates W and Py from the PCovR model fitted to the bootsrap samples of \code{x[ry,]} and \code{y[ry]}
#' \item Predicts y values based on \code{x[wy,]} and adds noise.
#' }
#' 
#' @author Edoardo Costantini, 2022
#' @references
#'
#' De Jong, S., & Kiers, H. A. (1992). Principal covariates regression: part I. 
#' Theory. Chemometrics and Intelligent Laboratory Systems, 14(1-3), 155-164.
#' 
#' Vervloet, M., Kiers, H. A., Van den Noortgate, W., & Ceulemans, E. (2015). 
#' PCovR: An R package for principal covariates regression. Journal of Statistical Software, 65, 1-14.
#'
#' @family univariate imputation functions
#' @keywords datagen
#' @export
mice.impute.pcovr <- function(y, ry, x, wy = NULL, npcs = 1L, ...) {

    # Set up
    install.on.demand("PCovR", ...)
    if (is.null(wy)) wy <- !ry

    # Scale X
    x_sc <- scale(x, center = TRUE, scale = TRUE)

    # Take bootstrap sample for model uncertainty
    n1 <- sum(ry)
    s <- sample(n1, n1, replace = TRUE)
    dotxobs <- x_sc[ry, , drop = FALSE][s, ]
    dotyobs <- y[ry][s] - mean(y[ry][s])
    xmis <- x_sc[wy, ]
    
    # Compute error ratio components
    lm_mod <- lm(dotyobs ~ -1 + dotxobs)
    ery <- 1 - summary(lm_mod)$r.squared

    # Compute Erx
    svd_erx <- svd(dotxobs)
    vec <- 1:npcs
    vec <- c(vec[1] - 1, vec, vec[length(vec)] + 1)
    VAF <- c(0, cumsum(svd_erx$d^2) / sum(svd_erx$d^2))
    VAF <- VAF[vec + 1]
    scr <- array(NA, c(1, length(vec)))
    for (u in 2:(length(vec) - 1)) {
        scr[, u] <- (VAF[u] - VAF[u - 1]) / (VAF[u + 1] - VAF[u])
    }
    erx <- 1 - VAF[which.max(scr)]

    # Find alpha ML
    alpha <- sum(dotxobs^2) / (sum(dotxobs^2) + sum(dotyobs^2) * erx / ery)

    # Estiamte PCovR on observed data
    Hx <- dotxobs %*% solve(t(dotxobs) %*% dotxobs) %*% t(dotxobs)
    G <- alpha * dotxobs %*% t(dotxobs) / sum(dotxobs^2) + (1 - alpha) * Hx %*% dotyobs %*% t(dotyobs) %*% Hx / sum(dotyobs^2)
    EG <- eigen(G) # eigen-decomposition of matrix
    Ts <- EG$vectors[, 1:npcs]
    W <- solve(t(dotxobs) %*% dotxobs) %*% t(dotxobs) %*% Ts
    Py <- t(W) %*% t(dotxobs) %*% dotyobs

    # Predict new observations
    y_hat <- mean(y[ry]) + xmis %*% W %*% Py

    # Compute residual standard error (sd of residuals with df as denominator) for PLS
    res_ss <- sum((y[ry] - mean(y[ry]) + dotxobs %*% W %*% Py)^2)
    res_df <- nrow(dotxobs) - npcs - 1
    sigma <- sqrt(res_ss / res_df)

    # Add noise to make imputations
    imputes <- y_hat + rnorm(sum(wy)) * sigma

    # Return
    return(imputes)
}
