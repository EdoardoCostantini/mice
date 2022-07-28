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
    install.on.demand("pls", ...)
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
      ncomp = nlvs,
      validation = "none",
      center = TRUE
    )

    # Compute degrees of freedom
    if(DoF == "naive"){
        res_df <- nrow(dotxobs) - nlvs
    }
    if(DoF == "kramer"){
        # Extract fitted values for all numbers of pcs
        Yhat <- predict(pls_out, ncomp = nlvs)

        # Compute DoFs
        DoF_plsr <- .dofPLS(
            X = dotxobs,
            y = dotyobs,
            q = nlvs,
            TT = apply(pls::scores(pls_out), 2, function(j) j / sqrt(sum(j^2))),
            Yhat = Yhat
        )

        res_df <- nrow(dotxobs) - DoF_plsr
        # If the computation fails, say so
        if(is.na(res_df) | is.nan(res_df)){
            stop(
              paste0(
                "Could not compute DoF. Try using a smaller npcs or DoF naive computation."
              )
            )
        }
    }

    # Compute residual sum of squares
    res_ss <- sum((stats::resid(pls_out)[, , nlvs])^2)

    # Compute sigma
    sigma <- sqrt(res_ss / res_df)

    # Predict new data
    pls_pred <- as.matrix(
        predict(pls_out, newdata = x[!ry, , drop = FALSE], ncomp = nlvs)
    )

    # Add noise for imputation uncertainty
    imputes <- pls_pred + rnorm(sum(wy)) * sigma

    # Return
    return(imputes)
}

# Degrees of freedom for supervised derived input models
.dofPLS <- function(X, y, q = 1, TT, Yhat){
    # Example inputs
    # X = scale(mtcars[, -1])
    # y = mtcars[, 1]
    # q <- 3 # desired component / latent variable
    # TT <- linear.pls.fit(X, y, m, DoF.max = DoF.max)$TT # normalizezs PC scores
    # Yhat <- linear.pls.fit(X, y, m, DoF.max = DoF.max)$Yhat[, (q + 1)]

    # Body
    TT <- TT[, 1:q, drop = FALSE]
    m <- q
    DoF.max <- ncol(X) + 1
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
    KY <- .krylov(K, K %*% y, m)
    lambda <- eigen(K)$values
    tr.K <- vector(length = m)
    for (i in 1:m) {
        tr.K[i] <- sum(lambda^i)
    }
    BB <- t(TT) %*% KY
    BB[row(BB) > col(BB)] <- 0
    b <- t(TT) %*% y
    Binv <- backsolve(BB, diag(m))
    tkt <- 0
    ykv <- 0
    KjT <- array(dim = c(q, n, m))

    dummy <- TT
    for (i in 1:q) {
        dummy <- K %*% dummy
        KjT[i, , ] <- dummy
    }

    Binvi <- Binv[1:q, 1:q, drop = FALSE]
    ci <- Binvi %*% b[1:q]
    Vi <- TT[, 1:q, drop = FALSE] %*% t(Binvi)
    trace.term <- sum(ci * tr.K[1:q])
    ri <- y - Yhat

    for (j in 1:q) {
        KjTj <- as.matrix(KjT[j, , ])
        tkt <- tkt + ci[j] * sum(diag((t(TT[, 1:q, drop = FALSE]) %*%
            KjTj[, 1:q, drop = FALSE])))
        ri <- K %*% ri
        ykv <- ykv + sum(ri * Vi[, j])
    }

    DoF <- trace.term + q - tkt + ykv
    DoF <- ifelse(DoF > DoF.max, DoF.max, DoF)
    DoF <- DoF + 1
    DoF
}

.krylov <- function(A, b, m) {
    K <- matrix(, length(b), m)
    dummy <- b
    for (i in 1:m) {
        K[, i] <- dummy
        dummy <- A %*% dummy
    }
    return(K)
}