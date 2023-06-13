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

# Factor as right levels
testthat::expect_equal(levels(imps_t1), levels(y))
