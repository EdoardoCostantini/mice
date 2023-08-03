# Test 1: Simple use -----------------------------------------------------------

# Define the test context
context("mice.impute.gspcr.logreg: simple use")

# Load example data
data(GSPCRexdata, package = "gspcr")

# Set a seed
set.seed(123)

# define an example predictor set
x <- GSPCRexdata$X$cont

# define an example dependent variable
y <- GSPCRexdata$y$bin

# make missingness
y[sample(1:nrow(iris), nrow(iris) * .3)] <- NA

# Define response indicator
ry <- !is.na(y)

# Define missingness indicator
wy <- !ry

# Use univariate imputation model
imps_t1 <- mice.impute.gspcr.logreg(y, ry, x)

# Right number of imputations
testthat::expect_equal(length(imps_t1), sum(wy))

# Imputations are in a factor
testthat::expect_true(is.factor(imps_t1))

# Correct number of imputations
testthat::expect_equal(length(imps_t1), sum(wy))

# Factor as right levels
testthat::expect_equal(levels(imps_t1), levels(y))

# Test 2: Perfect prediction ---------------------------------------------------

set.seed(20042)

# Generate some predictors
n <- 1e2
p <- 10
Sigma <- matrix(.7, nrow = p, ncol = p)
diag(Sigma) <- 1
x <- data.frame(MASS::mvrnorm(n, rep(0, p), Sigma))

# Create a dv with a perfect predictor
y <- factor(x[, 1] < 0, labels = c("y", "n"))

# Missing values
y[sample(1:n, n * .3)] <- NA

# Define response indicator
ry <- !is.na(y)

# Define missingness indicator
wy <- !ry

# Run imputation 
imps <- suppressMessages(mice.impute.gspcr.logreg(y = y, ry = ry, x = x))

# Returns imputations as expected (perfect prediction handled within gspcr)
testthat::expect_true(length(imps) == sum(wy))

# Test 3: Escape with same input if bootstrap causes constant ------------------

# Force constant y to mimic bootstrap result
y_no <- y[y == "n" | is.na(y)]
ry_no <- !is.na(y_no)
x_no <- x[y == "n" | is.na(y), ]

# Run imputation
imps <- mice.impute.gspcr.logreg(y = y_no, ry = ry_no, x = x_no)

# Test is the same as observed values of input imputations as expected (perfect prediction handled within gspcr)
testthat::expect_equal(imps, y_no[ry_no])