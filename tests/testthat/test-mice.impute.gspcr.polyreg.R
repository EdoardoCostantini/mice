# Test 1: Simple use -----------------------------------------------------------

# Define the test context
context("mice.impute.gspcr.polr: simple use")

# Load example data
data(GSPCRexdata, package = "gspcr")

# Set a seed
set.seed(20230614)

# define an example predictor set
x <- GSPCRexdata$X$cont

# define an example dependent variable
y <- GSPCRexdata$y$cat

# make missingness
y[sample(1:length(GSPCRexdata$y$ord), length(GSPCRexdata$y$ord) * .1)] <- NA

# Define response indicator
ry <- !is.na(y)

# Define missingness indicator
wy <- !ry

# Use univariate imputation model
imps_t1 <- mice.impute.gspcr.polyreg(y, ry, x, nthrs = 3)

# Impute with standard polr
imp_polyreg_standard <- mice.impute.polyreg(y, ry, x)

# Test: returns same class as standard polr
testthat::expect_equal(class(imp_polr_standard), class(imps_t1))