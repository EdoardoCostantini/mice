# Define a tolerance for acceptable difference
tol <- 1e-5

# Test 1: Simple use -----------------------------------------------------------

# Define the test context
context("revert.factors: simple use")

# Define an example dataset with shuffled order
preds <- iris[c(51:150, 1:50), -1]

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

# Test 3: Ordered factors are recovered ----------------------------------------

# Define the test context
context("revert.factors: ordered factors")

# Define a variable of interest to be transformed into an ordered factor
v <- iris$Sepal.Length

# Create ordered factors from iris data
K <- 5

# Create lags
lags <- rep(abs(min(v) - max(v)) / K, (K - 1))

# Define the breakpoints for v
breaks <- c(cumsum(c(minimum = min(v), fixed = lags)), maximum = max(v))

# Cut v with the given brakes
x_dis <- as.numeric(cut(x = iris$Sepal.Length, breaks = breaks, include.lowest = TRUE))

# Make an ordered factor
SL_ordered <- factor(x_dis, ordered = TRUE)

# Create a model matrix as used in mice
x <- obtain.design(SL_ordered)

# Run the imputation method with this model matrix as X input
x_mm_reverse <- revert.factors(x)

# Compute frequnecy table
freq_table <- table(x_mm_reverse[, 2], SL_ordered)

# Factor levels order is the same factor is recovered exactly
testthat::expect_true(all(freq_table[upper.tri(freq_table)] == 0))

# Same number of levels before and after
testthat::expect_true(ncol(freq_table) == nrow(freq_table))