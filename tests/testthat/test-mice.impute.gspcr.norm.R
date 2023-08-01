# Test 1: Simple use -----------------------------------------------------------

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

# Test 2: Categorical predictors -----------------------------------------------

# Define the test context
context("mice.impute.gspcr.norm: categorical predictors")

# Define an example dataset
x <- iris[, -1]

# Create an additional binary variable
x[, 1] <- factor(x[, 1] < mean(x[, 1]), labels = c("Y", "N"))

# define an example dependent variable
y <- iris[, 1]

# make missingness
y[sample(1:nrow(iris), nrow(iris) * .3)] <- NA

# Define response indicator
ry <- !is.na(y)

# Define missingness indicator
wy <- !ry

# Use univariate imputation model
imps_t2 <- mice.impute.gspcr.norm(y, ry, x)

