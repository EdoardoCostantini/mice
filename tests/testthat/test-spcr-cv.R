context("SPCR cross validation")

# Functions --------------------------------------------------------------------

generateXTP <- function(I, J, VAFr = c(.5, .4, .2), VAFsum = 100, CPVE = 0.9) {
    # Internals -------------------------------------------------------------

    # I    = 100 # sample size
    # J    = 9 # number of variables
    # VAFr = c(.5, .3, .2) # relative variance of each components
    # VAFsum = 100 # total variance of the components
    # CPVE = 0.9 # proportion of explained variance by the R components

    # Body ------------------------------------------------------------------
    # Number of components
    R <- length(VAFr)

    # Random sample U
    U <- matrix(
        data = rnorm(I * R),
        nrow = I,
        ncol = R
    )
    U <- scale(U, center = TRUE, scale = FALSE)
    U <- orthmat(U, verbose = FALSE)
    U <- normmat(U)

    # Random sample P
    V <- matrix(
        data = runif(J * R),
        nrow = J,
        ncol = R
    )
    V <- orthmat(V, verbose = FALSE)
    P <- normmat(V)

    # Define D
    D <- diag(c(VAFsum * VAFr))

    # Create X
    Xtrue <- U %*% D %*% t(P)

    # sample from normal distribution (Ex = Error of X)
    Ex <- MASS::mvrnorm(n = I, mu = rep(0, J), Sigma = diag(J))

    # centering and scaling the Ex matrix
    Ex <- scale(Ex, center = TRUE, scale = FALSE)

    # sum of squares
    ssqXtrue <- sum(Xtrue^2)
    ssqEx <- sum(Ex^2)

    # Compute noise scaling constant
    Escale <- sqrt(ssqXtrue * (1 - CPVE) / (CPVE * ssqEx))

    # Add scaled noise
    X <- Xtrue + Escale * Ex

    # Scale data for estimation
    X <- scale(X)

    # Define outputs
    return(list(
        X = data.frame(X),
        T = U %*% D,
        P = P,
        U = U,
        D = D
    ))
}



generateDV <- function(X = matrix(), R2 = 0.90, beta = 1) {
    # Internals -------------------------------------------------------------

    # X     = matrix(rnorm(1e3 * 4), nrow = 1e3, ncol = 4)
    # R2    = .9
    # beta  = 1

    # Body ------------------------------------------------------------------
    # Generate a dependent variable on true line
    y_true <- as.vector(X %*% rep(beta, ncol(X)))

    # Generate random error
    error <- rnorm(nrow(X))

    # Make the error orthogonal to the Xs
    e_ortho <- orthmat(cbind(X, error))[, "error"]

    # sum of squares
    ssqy <- sum(X^2)
    ssqe <- sum(e_ortho^2)

    # Rescale noise to desired level
    e_scale <- sqrt(ssqy * (1 - R2) / (R2 * ssqe))

    # Generate samples for y
    y <- y_true + e_scale * e_ortho

    # What to return
    return(y)
}

normmat <- function(X) {
    # Internals -------------------------------------------------------------

    # X    = matrix(rnorm(1e3 * 4), nrow = 1e3, ncol = 4)

    # Body ------------------------------------------------------------------
    X <- apply(X, 2, function(j) j / sqrt(sum(j^2)))
    return(X)
}

orthmat <- function(X, verbose = FALSE) {
    # Internals -------------------------------------------------------------

    # X    = matrix(rnorm(1e3 * 4), nrow = 1e3, ncol = 4)

    # Body ------------------------------------------------------------------
    for (i in 2:ncol(X)) {
        for (j in 1:(i - 1)) {
            if (verbose == TRUE) {
                print(paste0("Adjusting piar ", i, "-", j))
            }
            A <- X[, j]
            b <- X[, i]

            # Find projection of b on A
            B <- as.vector(b - (t(A) %*% b / t(A) %*% A) %*% A)

            # Replace in original X the orthogonalized columns
            X[, j] <- A
            X[, i] <- B
        }
    }
    return(X)
}

# Generate data ----------------------------------------------------------------

library(superpc)
library(reshape2)
library(ggplot2)
library(PCovR)
library(nFactors)
library(patchwork)
library(dplyr)

set.seed(2026)

# Define the desired number of components
Q <- 5

# Desired sample size
N <- 1e3

# Desired number of variables
P <- 50

# Generate data
XTP <- generateXTP(
    I <- N, # sample size
    J <- P, # number of variables
    VAFr <- diff(seq(0, 1, len = Q + 1)), # relative variance of each components
    VAFsum <- 100, # total variance of the components
    CPVE <- 0.9 # proportion of explained variance by the R components
)
X <- XTP$X

# Attach junk predictors?
X <- cbind(X, junk = matrix(rnorm(N * P), nrow = N, ncol = P))

# Compute PCA with prcomp
PCX <- prcomp(X)

# Extract eigenvalues
eigenvalues <- PCX$sdev^2

# Cumulative proportion of explained variance
cumsum(prop.table(eigenvalues))

# Screeplot
plotuScree(x = eigenvalues)

# Non-graphical solutions
nScree(x = eigenvalues)

# generate DV based on the component scores
y <- generateDV(
    X = XTP$T[, c(1, 2), drop = FALSE],
    R2 = 0.99,
    beta = 1
)

# Check the model is as expected
summary(lm(y ~ XTP$T))

# Estimate PCovR model ---------------------------------------------------------

PCovR_out <- pcovr(
    X = X,
    Y = as.data.frame(y),
    rot = "none",
    modsel = "seq" # fastest option
)

# Npcs selected by PCovR
PCovR_out$R

# Value of alpha selected
PCovR_out$alpha

# Estimate SPCR ----------------------------------------------------------------

# Define a train data
data.train <- list(x = t(as.matrix(X)), y = y, featurenames = colnames(X))

# Train the model (computes the scores for each feature)
train.obj <- superpc.train(
    data = data.train,
    type = "regression"
)

# Cross-validate the model
cv.obj <- superpc.cv(
    fit = train.obj,
    data = data.train,
    n.fold = 10,
    n.threshold = 20,
    n.components = 5
)

# Objects 
cv.obj$thresholds
cv.obj$scor

# Plot the cross-validation curves
superpc.plotcv(cv.obj)

# Estimate SPCR with my implementation of the superpc method -------------------

# Obtain CV estimates
out_F <- .spcr.cv(
    dv = y,
    ivs = X,
    fam = "gaussian",
    nthrs = 20,
    maxnpcs = 10,
    K = 10,
    test = c("LRT", "F", "MSE")[2],
    thrs = c("LLS", "pseudoR2", "normalized")[3],
    min.features = 1,
    max.features = ncol(X)
)

# Obtain CV estimates
out_LRT <- .spcr.cv(
    dv = y,
    ivs = X,
    fam = "gaussian",
    nthrs = 20,
    maxnpcs = 10,
    K = 5,
    test = c("LRT", "F", "AIC", "BIC", "PR2")[1],
    thrs = c("LLS", "pseudoR2", "normalized")[3],
    min.features = 1,
    max.features = ncol(X)
)

# Obtain CV estimates
out_BIC <- .spcr.cv(
    dv = y,
    ivs = X,
    fam = "gaussian",
    nthrs = 20,
    maxnpcs = 10,
    K = 5,
    test = c("LRT", "F", "AIC", "BIC")[4],
    thrs = c("LLS", "pseudoR2", "normalized")[3],
    min.features = 1,
    max.features = ncol(X)
)

# Check scores
t(out_F$scor)
t(out_LRT$scor)

# Check thresholds
data.frame(
    F = c(
        thr = which(out_F$thr.cv == round(out_F$thr, 3)),
        Q = out_F$Q.cv
    ),
    LRT = c(
        thr = which(out_LRT$thr.cv == round(out_LRT$thr, 3)),
        Q = out_LRT$Q.cv
    ),
    BIC = c(
        thr = which(out_BIC$thr.cv == round(out_BIC$thr, 3)),
        Q = out_BIC$Q.cv
    )
)

# Plot trends in a similar way to the package
df <- reshape2::melt(out_F$scor) # the function melt reshapes it from wide to long
plot_norm <- ggplot(df, aes(Var2, value, group = factor(Var1), label = factor(Var1))) +
    geom_line() +
    geom_point() +
    geom_label() +
    theme_bw()

# Plot trends in a similar way to the package
df <- reshape2::melt(out_LRT$scor) # the function melt reshapes it from wide to long
plot_LRT <- ggplot(df, aes(Var2, value, group = factor(Var1), label = factor(Var1))) +
    geom_line() +
    geom_point() +
    geom_label() +
    theme_bw()

# Plot trends in a similar way to the package
df <- reshape2::melt(out_BIC$scor) # the function melt reshapes it from wide to long
plot_BIC <- ggplot(df, aes(Var2, value, group = factor(Var1), label = factor(Var1))) +
    geom_line() + 
    geom_point() + 
    geom_label() + 
    # coord_cartesian(ylim = c(min(out_BIC$scor, na.rm = TRUE), 0)) + 
    theme_bw()

plot_norm / plot_LRT + plot_BIC

# Function based on PCKG

spcr.out.pkg <- .spcr.cv.pkg(
    dv = y,
    ivs = X,
    fam = "gaussian",
    nthrs = 20,
    maxnpcs = 5,
    K = 10
)

# Check thresholds are the same
cbind(
    pack = cv.obj$thresholds,
    aspack = spcr.out.pkg$thr,
    mine = out$thr
)

# Check scores are the same
t(cv.obj$scor)
t(spcr.out.pkg$scor)
t(out$scor)

spcr.out.pkg$thr

t(out_LRT$scor)
t(out_F$scor)

df <- reshape2::melt(spcr.out.pkg$scor) # the function melt reshapes it from wide to long
ggplot(df, aes(Var2, value, group = factor(Var1), label = factor(Var1))) +
    geom_line() + 
    geom_point() + 
    geom_label() + 
    theme_bw()

# Estimate SPCR with my full CV method -----------------------------------------

# Use the functions with a given method
out1 <- .spcr.cv.full(
    dv = y,
    ivs = X,
    fam = "gaussian",
    nthrs = 20,
    maxnpcs = 10,
    K = 10,
    test = "F",
    thrs = "normalized",
    min.features = 1,
    max.features = ncol(X)
)

# Vector of desired methods
vmeth <- c("LRT", "F", "AIC", "BIC", "PR2", "MSE")#[1:2]

# Estimate all of the methods desired
out <- lapply(
    vmeth, 
    function(meth) {
    .spcr.cv.full(
        dv = y,
        ivs = X,
        fam = "gaussian",
        nthrs = 20,
        maxnpcs = 10,
        K = 10,
        test = meth,
        thrs = "normalized",
        min.features = 1,
        max.features = ncol(X),
        oneSE = TRUE
    )
})

# Give names for simplicity
names(out) <- vmeth

# Create plots for all of the desired methods
plots <- lapply(
    1:length(vmeth),
    function(meth) {

        # Plot trends in a similar way to the package
        df <- reshape2::melt(out[[meth]]$scor) # the function melt reshapes it from wide to long

        # Add error bars
        df$low <- reshape2::melt(out[[meth]]$scor.lwr)[, "value"]
        df$high <- reshape2::melt(out[[meth]]$scor.upr)[, "value"]

    # Make plot
    store_plot <- df %>%
        filter(Var1 %in% unique(Var1)) %>%
        ggplot(
            aes(Var2, value, group = factor(Var1), label = factor(Var1))
        ) +
        geom_line() +
        # geom_point() +
        geom_errorbar(aes(ymin = low, ymax = high),
            width = .2
        ) +
        geom_label() +
        ggtitle(paste0("Measure used for cv: ", vmeth[meth])) +
        theme_bw()

        # And return the plot
        return(store_plot)
    }
)

# Plots
(plots[[1]] + plots[[2]] + plots[[3]]) / (plots[[4]] + plots[[5]] + plots[[6]])

# default Solutions
res <- lapply(
    1:length(vmeth),
    function(meth) {
        c(
            thr_value = out[[meth]]$thr.cv,
            thr_number = which(out[[meth]]$thr.cv == round(out[[meth]]$thr, 3)),
            Q = out[[meth]]$Q.cv
        )
    }
)

# 1se solutions
res.1se <- lapply(
    1:length(vmeth),
    function(meth) {
        c(
            thr_value = out[[meth]]$thr.cv.1se,
            thr_number = which(out[[meth]]$thr.cv.1se == round(out[[meth]]$thr, 3)),
            Q = out[[meth]]$Q.cv.1se
        )
    }
)

# Present them neatly
res <- do.call(rbind, res)
rownames(res) <- vmeth
res.1se <- do.call(rbind, res.1se)
rownames(res.1se) <- vmeth

res
res.1se

# SPCR using CHull -------------------------------------------------------------

out_chull <- .bic.selection(
    dv = y,
    ivs = X,
    fam = "gaussian",
    nthrs = 20,
    maxnpcs = 10,
    test = "BIC",
    thrs = c("LLS", "pseudoR2", "normalized")[3],
    min.features = 1,
    max.features = ncol(X)
)

out_chull$thr
out_chull$thr.cv
out_chull$Q.cv

# Plot scores
df <- reshape2::melt(out_chull$scor) # the function melt reshapes it from wide to long
ggplot(df, aes(Var2, value, group = factor(Var1), label = factor(Var1))) +
    geom_line() +
    geom_point() +
    geom_label() +
    theme_bw()

# Decision based on CHull
out_chull$thr.cv
out_chull$Q.cv
length(out_chull$pred.active)
plot(out_chull$solution)

# Check out pcovr data ---------------------------------------------------------

# Load data
data(alexithymia)

# PCA
# Compute PCA with prcomp
PCX <- prcomp(alexithymia$X)

# Extract eigenvalues
eigenvalues <- PCX$sdev^2

# Cumulative proportion of explained variance
cumsum(prop.table(eigenvalues))

# Screeplot
plotuScree(x = eigenvalues)

# Non-graphical solutions
nScree(x = eigenvalues)

# PCovR
results <- pcovr(alexithymia$X, alexithymia$Y[, 2, drop = FALSE])
summary(results)
plot(results)

# SPCR
out_alex <- .spcr.cv(
    dv = alexithymia$Y[, 2],
    ivs = alexithymia$X,
    fam = "gaussian",
    nthrs = 20,
    maxnpcs = 10,
    K = 5,
    test = c("LRT", "F", "MSE", "BIC")[4],
    thrs = c("LLS", "pseudoR2", "normalized")[2],
    min.features = 1,
    max.features = ncol(alexithymia$X)
)

out_alex$Q.cv
out_alex$thr.cv
out_alex$pred.map
out_alex$pred.active

df <- reshape2::melt(out_alex$scor) # the function melt reshapes it from wide to long
ggplot(df, aes(Var2, value, group = factor(Var1), label = factor(Var1))) +
    geom_line() +
    geom_point() +
    geom_label() +
    theme_bw()

library("survival")
# Create the simplest test data set
test1 <- list(
    time = c(4, 3, 1, 1, 2, 2, 3),
    status = c(1, 1, 1, 0, 1, 1, 0),
    x = c(0, 2, 1, 1, 1, 0, 0),
    sex = c(0, 0, 0, 0, 1, 1, 1)
)

# Fit a stratified model
survival::coxph(Surv(time, status) ~ x + strata(sex), test1)
# Create a simple data set for a time-dependent model
test2 <- list(
    start = c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
    stop = c(2, 3, 6, 7, 8, 9, 9, 9, 14, 17),
    event = c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
    x = c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0)
)
summary(coxph(Surv(start, stop, event) ~ x, test2))

out <- coxph(Surv(start, stop, event) ~ x, test2)


# Survival data plot?

# generate some synthetic survival data.
#
# there are 1000 features and 60 samples
#  the outcome is highly correlated with the first principal component of the
#  first 80 features


set.seed(464)

p <- 100
n <- 1e3

x <- matrix(rnorm(p * n), ncol = n)
v1 <- svd(x[1:80, ])$v[, 1]
X <- t(x)
y <- 2 + 5 * v1 + .05 * rnorm(n)

data <- list(x = x, y = y, featurenames = featurenames)

# train  the model. This step just computes the  scores for each feature

train.obj <- superpc.train(data, type = "regression")

# note for regression (non-survival) data, we leave the component "censoring.status"
# out of the data object, and call superpc.train with type="regression".
# otherwise the superpc commands are all the same

# cross-validate the model

cv.obj <- superpc.cv(train.obj, data, n.components = 5)

# plot the cross-validation curves. From this plot we see that the 1st
# principal component is significant and the best threshold  is around 0.7

superpc.plotcv(cv.obj)


# ORIGINAL METHOD: Cross-validation of thrs and npcs always lead to highest npcs ------

library(superpc)
set.seed(464)


x <- matrix(rnorm(1000 * 100), ncol = 100)
v1 <- svd(x[1:80, ])$v[, 1]

y <- 2 + 5 * v1 + .05 * rnorm(100)

xtest <- x
ytest <- 2 + 5 * v1 + .05 * rnorm(100)
censoring.status <- sample(c(rep(1, 80), rep(0, 20)))
censoring.status.test <- sample(c(rep(1, 80), rep(0, 20)))

featurenames <- paste("feature", as.character(1:1000), sep = "")



# create train and test data objects. censoring.status=1 means the event occurred;
#  censoring.status=0 means censored

data <- list(x = x, y = y, censoring.status = censoring.status, featurenames = featurenames)
data.test <- list(x = xtest, y = ytest, censoring.status = censoring.status.test, featurenames = featurenames)


# train  the model. This step just computes the  scores for each feature

train.obj <- superpc.train(data, type = "survival")

# note for regression (non-survival) data, we leave the component "censoring.status"
# out of the data object, and call superpc.train with type="regression".
# otherwise the superpc commands are all the same



# cross-validate the model

cv.obj <- superpc.cv(train.obj, data)

# plot the cross-validation curves. From this plot we see that the 1st
# principal component is significant and the best threshold  is around 0.7

superpc.plotcv(cv.obj)


# See pdf version of the cv plot

# here we have the luxury of  test data, so we can compute the  likelihood ratio statistic
# over the test data and plot them. We see that the threshold of 0.7
# works pretty well
#

lrtest.obj <- superpc.lrtest.curv(train.obj, data, data.test, n.components = 3)

superpc.plot.lrtest(lrtest.obj)
