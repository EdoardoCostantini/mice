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
mice.impute.pcovr <- function(y, ry, x, wy = NULL, npcs = 1L, DoF = "naive", ...) {

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
        # Extract fitted values for all numbers of pcs up untill npcs
        Yhat <- sapply(1:ncol(dotxobs), function(j){
            mean(y[ry]) + dotxobs %*% W[, 1:j, drop = FALSE] %*% Py[1:j, , drop = FALSE]
        })
        # Compute DoFs
        DoFs <- .dofPLS(
            X = dotxobs,
            y = dotyobs,
            TT = Ts,
            Yhat = Yhat,
            m = ncol(dotxobs),
            DoF.max = ncol(dotxobs) + 1
        )
        res_df <- nrow(dotxobs) - DoFs[npcs + 1]
    }

    # Compute sigma    
    sigma <- sqrt(res_ss / res_df)

    # Add noise to make imputations
    imputes <- y_hat + rnorm(sum(wy)) * sigma

    # Return
    return(imputes)
}

# Degrees of freedom for supervised derived input models
.dofPLS <- function(X, y, TT, Yhat, m = ncol(X), DoF.max = ncol(X) + 1){
    # Example inputs
    # X = scale(mtcars[, -1])
    # y = mtcars[, 1]
    # m = ncol(X)
    # DoF.max = m + 1
    # TT <- linear.pls.fit(X, y, m, DoF.max = DoF.max)$TT # normalizezs PC scores
    # Yhat <- linear.pls.fit(X, y, m, DoF.max = DoF.max)$Yhat[, 2:(m + 1)]

    # Body
    n <- nrow(X)

    # Scale data
    mean.X <- apply(X, 2, mean)
    sd.X <- apply(X, 2, sd)
    sd.X[sd.X == 0] <- 1
    X <- X - rep(1, nrow(X)) %*% t(mean.X)
    X <- X / (rep(1, nrow(X)) %*% t(sd.X))
    K <- X %*% t(X)

    # pls.dof
    DoF.max <- DoF.max - 1
    TK <- matrix(, m, m)
    KY <- krylov(K, K %*% y, m)
    lambda <- eigen(K)$values
    tr.K <- vector(length = m)
    for (i in 1:m) {
        tr.K[i] <- sum(lambda^i)
    }
    BB <- t(TT) %*% KY
    BB[row(BB) > col(BB)] <- 0
    b <- t(TT) %*% y
    DoF <- vector(length = m)
    Binv <- backsolve(BB, diag(m))
    tkt <- rep(0, m)
    ykv <- rep(0, m)
    KjT <- array(dim = c(m, n, m))
    dummy <- TT
    for (i in 1:m) {
        dummy <- K %*% dummy
        KjT[i, , ] <- dummy
    }
    trace.term <- rep(0, m)

    for (i in 1:m) {
        Binvi <- Binv[1:i, 1:i, drop = FALSE]
        ci <- Binvi %*% b[1:i]
        Vi <- TT[, 1:i, drop = FALSE] %*% t(Binvi)
        trace.term[i] <- sum(ci * tr.K[1:i])
        ri <- y - Yhat[, i]
        for (j in 1:i) {
            KjTj <- KjT[j, , ]
            tkt[i] <- tkt[i] + ci[j] * 
                sum(diag(t(TT[, 1:i, drop = FALSE]) %*% KjTj[, 1:i, drop = FALSE]))
            ri <- K %*% ri
            ykv[i] <- ykv[i] + sum(ri * Vi[, j])
        }
    }

    DoF <- trace.term + 1:m - tkt + ykv

    DoF[DoF > DoF.max] <- DoF.max
    DoF <- c(0, DoF) + 1
    DoF

}

krylov <- function(A, b, m) {
    K <- matrix(, length(b), m)
    dummy <- b
    for (i in 1:m) {
        K[, i] <- dummy
        dummy <- A %*% dummy
    }
    return(K)
}