#' Imputation by principal covariates regression
#'
#' Imputes univariate missing data by using principal covariates regression (PCovR).
#'
#' @aliases mice.impute.pcovr pcovr
#' @inheritParams mice.impute.norm.boot
#' @param npcs The number of principal components to extract for PC regression.
#' @param DoF The method to compute the degrees of freedom, either \code{"naive"}, or \code{"kramer"}.
#' @return Vector with imputed data, same type as \code{y}, and of length
#' \code{sum(wy)}
#' @details
#' Imputation of \code{y} by principal covariates regression (De Jong and Kiers, 1992).
#' The method consists of the following steps:
#' \enumerate{
#' \item For a given \code{y} variable under imputation, draw a bootstrap version y*
#' with replacement from the observed cases \code{y[ry]}, and stores in x* the
#' corresponding values from \code{x[ry, ]}.
#' \item Compute an optimal \code{alpha} value for the PCovR model on the basis 
#' of maximum likelihood principles described in Vervloet et. al. (2015).
#' \item Estimate W and Py from the PCovR model fitted to y* and x*
#' \item Predict y values based on \code{x[wy,]} and add noise.
#' }
#' 
#' Residual degrees of freedom are estimated by using the Krylov representation
#' described by Kramer and Sugiyama (2011) as implemented in the \code{plsdof} 
#' package, or with the naive approach (n - npcs - 1).
#' 
#' @author Edoardo Costantini, 2022
#' @references
#'
#' De Jong, S., & Kiers, H. A. (1992). Principal covariates regression: part I. 
#' Theory. Chemometrics and Intelligent Laboratory Systems, 14(1-3), 155-164.
#' 
#' Krämer, N., & Sugiyama, M. (2011). The degrees of freedom of partial least
#' squares regression. Journal of the American Statistical Association, 106(494), 697-705.
#' 
#' Vervloet, M., Kiers, H. A., Van den Noortgate, W., & Ceulemans, E. (2015). 
#' PCovR: An R package for principal covariates regression. Journal of Statistical Software, 65, 1-14.
#'
#' @family univariate imputation functions
#' @keywords datagen
#' @export
mice.impute.pcovr <- function(y, ry, x, wy = NULL, npcs = 1L, DoF = "kramer", ...) {

    # Set up
    install.on.demand("PCovR", ...)
    if (is.null(wy)) wy <- !ry
    
    # Take bootstrap sample for model uncertainty
    n1 <- sum(ry)
    s <- sample(n1, n1, replace = TRUE)
    dotxobs <- x[ry, , drop = FALSE][s, ]
    dotyobs <- y[ry][s] - mean(y[ry][s])
    xmis <- x[wy, ]

    # Scale Xs
    dotxobs <- scale(dotxobs, center = TRUE, scale = TRUE)
    xmis <- scale(x[wy, ],
        center = attributes(dotxobs)$`scaled:center`,
        scale = attributes(dotxobs)$`scaled:scale`
    )
    
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
    Ts <- EG$vectors [, 1:ncol(dotxobs), drop = FALSE]
    W <- solve(t(dotxobs) %*% dotxobs) %*% t(dotxobs) %*% Ts
    Py <- t(W) %*% t(dotxobs) %*% dotyobs
    
    # Predict new observations
    y_hat <- mean(y[ry]) + xmis %*% W[, 1:npcs, drop = FALSE] %*% Py[1:npcs, , drop = FALSE]

    # Compute residual sum of squares
    res_ss <- sum((y[ry] - mean(y[ry]) + dotxobs %*% W %*% Py)^2)

    # Compute degrees of freedom
    if (DoF == "naive") {
        res_df <- nrow(dotxobs) - npcs - 1
    }
    if (DoF == "kramer") {
        # Extract fitted values for all numbers of pcs
        Yhat <- mean(y[ry]) + dotxobs %*% W[, 1:npcs, drop = FALSE] %*% Py[1:npcs, , drop = FALSE]

        # Compute DoFs
        DoF_plsr <- .dofPLS(
            X = dotxobs,
            y = dotyobs,
            q = npcs,
            TT = Ts,
            Yhat = Yhat
        )
        res_df <- nrow(dotxobs) - DoF_plsr
    }

    # Compute sigma    
    sigma <- sqrt(res_ss / res_df)

    # Add noise to make imputations
    imputes <- y_hat + rnorm(sum(wy)) * sigma

    # Return
    return(imputes)
}