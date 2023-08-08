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
testthat::expect_equal(class(imp_polyreg_standard), class(imps_t1))

# Test 2: bootstrap leads to empty levels in input y ---------------------------

# define an example predictor set
x <- iris[, -ncol(iris)]

# define an example dependent variable
y <- iris$Species

# add an empty level
levels(y) <- c(levels(iris$Species), "empty")

# make missingness
y[sample(1:nrow(iris), nrow(iris) * .1)] <- NA

# Define response indicator
ry <- !is.na(y)

# Define missingness indicator
wy <- !ry

# Use univariate imputation model
imps <- tryCatch(
    expr = mice.impute.gspcr.polyreg(y, ry, x, nthrs = 3),
    error = function(e) {
        e
    }
)

# Test the outcome is not an error
testthat::expect_false("error" %in% class(imps))

# Test 3: Escape with same input if bootstrap causes constant ------------------

# Force constant y to mimic the bootstrap result
x_no <- iris[1:50, -ncol(iris)]
y_no <- iris[1:50, "Species"]
y_no[sample(1:50, 50 * .1)] <- NA
ry_no <- !is.na(y_no)

# Run imputation
imps <- mice.impute.gspcr.polyreg(y = y_no, ry = ry_no, x = x_no)

# Test is the same as observed values of input imputations as expected (perfect prediction handled within gspcr)
testthat::expect_equal(imps, y_no[ry_no])