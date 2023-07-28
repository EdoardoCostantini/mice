# Define a tolerance for acceptable difference
tol <- 1e-5

# Test 1: Simple use -----------------------------------------------------------

# Define the test context
context("revert.factors: simple use")

# Define an example dataset
preds <- iris[, -1]

# Create an additional binary variable
preds[, 1] <- factor(preds[, 1] < mean(preds[, 1]), labels = c("Y", "N"))

# define an example dependent variable
y <- iris[, 1]

# Create a model matrix as used in mice
x <- obtain.design(preds)

# Run the imputation method with this model matrix as X input
x_mm_reverse <- revert.factors(x)

# Nominal factor is recovered exactly
freq_table <- table(x_mm_reverse$Species, preds$Species)
testthat::expect_true(all(freq_table[upper.tri(freq_table)] == 0))

# Binary factor is recovered exactly
freq_table <- table(x_mm_reverse$Sepal.Width, preds$Sepal.Width)
testthat::expect_true(all(freq_table[upper.tri(freq_table)] == 0))

# Test 2: Nothing happens if there are no contrasts ----------------------------

# Define the test context
context("revert.factors: no contrasts")

# Define an example dataset
x <- iris[, -ncol(iris)]

# Create a model matrix as used in mice
x_mm <- obtain.design(x)

# Run the imputation method with this model matrix as X input
x_mm_reverse <- revert.factors(x_mm)[, -1]

# Check the two matrices are the same up until tolerance
testthat::expect_true(all(x - x_mm_reverse < tol))