# Test 2: Simple use -----------------------------------------------------------

# Define the test context
context("mice.impute.gspcr.norm: simple use")

# Set a seed
set.seed(123)

# define an example predictor set
x <- iris[, -1]

# define an example dependent variable
y <- iris[, 1]

# make missingness
y[sample(1:nrow(iris), nrow(iris) * .3)] <- NA

# Define response indicator
ry <- !is.na(y)

# Define missingness indicator
wy <- !ry

# Use univariate imputation model
imps_t1 <- mice.impute.gspcr.norm(y, ry, x)

