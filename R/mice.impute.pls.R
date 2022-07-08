#' Imputation by partial least square regression
#'
#' Imputes univariate missing data by using partial least square regression.
#'
#' @aliases mice.impute.pls pls
#' @inheritParams mice.impute.norm.boot
#' @param nlvs The number of latent variables to consider in the PLS regression.
#' @param DoF The method to compute the degrees of freedom, either \code{"naive"}, or \code{"kramer"}.
#' @return Vector with imputed data, same type as \code{y}, and of length
#' \code{sum(wy)}
#' @details
#' Imputation of \code{y} by partial least square regression (Wold, 1975).
#' The method consists of the following steps:
#' \enumerate{
#' \item For a given \code{y} variable under imputation, draw a bootstrap version y*
#' with replacement from the observed cases \code{y[ry]}, and stores in x* the
#' corresponding values from \code{x[ry, ]}.
#' \item PLSR is estimated on this sample for a given number of latent variables (nlvs).
#' \item Based on this PLSR model, predictions on the new data x[wy, ] are returned
#' \item Imputations are then obtained by adding noise to these predictions.
#' }
#' 
#' Noise is added by sampling from a standard normal distribution rescaled based 
#' on the residuals standard error of the PLSR model.
#' The residual degrees of freedom are computed based on the Krylov representation 
#' described by Kramer and Sugiyama (2011) as implemented in the R package \code{plsdof} 
#' (Krämer and Sugiyama, 2011), or with the naive approach (n - nlvs - 1).
#' 
#' @author Edoardo Costantini, 2022
#' @references
#' 
#' Krämer, N., & Sugiyama, M. (2011). The degrees of freedom of partial least 
#' squares regression. Journal of the American Statistical Association, 106(494), 697-705.
#' 
#' Wold, H. (1975). Path models with latent variables: The NIPALS approach. 
#' In Quantitative sociology (pp. 307-357). Academic Press.
#'
#' @family univariate imputation functions
#' @keywords datagen
#' @export
mice.impute.pls <- function(y, ry, x, wy = NULL, nlvs = 1L, DoF = "naive", ...) {

    # Set up    
    install.on.demand("plsdof", "pls", ...)
    if (is.null(wy)) wy <- !ry

    # Take bootstrap sample for model uncertainty
    n1 <- sum(ry)
    s <- sample(n1, n1, replace = TRUE)
    dotxobs <- x[ry, , drop = FALSE][s, ]
    dotyobs <- y[ry][s]

    # Train PLS on observed (bootsrap) data
    pls_out <- pls::plsr(
        y ~ .,
        data = data.frame(y = dotyobs, dotxobs),
        ncomp = ncol(dotxobs),
        validation = "none",
    )

    # Predict new data
    pls_pred <- as.matrix(
        predict(pls_out, newdata = x[!ry, , drop = FALSE], ncomp = nlvs)
    )

    # Compute residual sum of squares
    res_ss <- sum((resid(pls_out)[, , nlvs])^2)

    # Compute degrees of freedom
    if(DoF == "naive"){
        res_df <- nrow(dotxobs) - nlvs - 1
    }
    if(DoF == "kramer"){
        # Extract fitted values for all numbers of pcs
        Yhat <- sapply(1:ncol(dotxobs), function(j) {
            predict(pls_out, ncomp = j)
        })
        # Compute DoFs
        DoF_plsr <- .dofPLS(
            X = dotxobs,
            y = dotyobs,
            TT = apply(pls::scores(pls_out), 2, function(j) j / sqrt(sum(j^2))),
            Yhat = Yhat,
            m = ncol(dotxobs),
            DoF.max = ncol(dotxobs) + 1
        )
        res_df <- nrow(dotxobs) - DoF_plsr[nlvs + 1]
    }

    # Compute sigma
    sigma <- sqrt(res_ss / res_df)

    # Add noise for imputation uncertainty
    imputes <- pls_pred + rnorm(sum(wy)) * sigma

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