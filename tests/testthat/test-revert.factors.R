# Define a tolerance for acceptable difference
tol <- 1e-5

# Test 1: Simple use -----------------------------------------------------------

# Define the test context
context("revert.factors: simple use")

# Define an example dataset with shuffled order
preds <- iris[c(51:150, 1:50), -1]

# Create an additional binary variable
preds[, 1] <- factor(preds[, 1] < mean(preds[, 1]), labels = c("Y", "N"))

# Create a model matrix as used in mice
x <- obtain.design(preds)

# Store contrasts
x_constrasts <- attributes(x)$contrasts

# Remove the intercept
x <- x[, -1]

# Add required variable info
attributes(x)$contrasts <- x_constrasts
attributes(x)$assign_better <- mm.column.variable.map(data = preds)

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
SL_ordered <- data.frame(v1 = factor(x_dis, ordered = TRUE))

# Create a model matrix as used in mice
x <- obtain.design(SL_ordered)

# Store contrasts
x_constrasts <- attributes(x)$contrasts

# Remove the intercept
x <- x[, -1]

# Add required variable info
attributes(x)$contrasts <- x_constrasts
attributes(x)$assign_better <- mm.column.variable.map(data = SL_ordered)

# Run the imputation method with this model matrix as X input
x_mm_reverse <- revert.factors(x)

# Compute frequnecy table
freq_table <- table(x_mm_reverse$v1, SL_ordered$v1)

# Factor levels order is the same factor is recovered exactly
testthat::expect_true(all(freq_table[upper.tri(freq_table)] == 0))

# Same number of levels before and after
testthat::expect_true(ncol(freq_table) == nrow(freq_table))

# Test: variables with nested names --------------------------------------------

# Define an example dataset with "nested" variable names
preds <- data.frame(
    v1 = factor(iris$Species, ordered = TRUE),
    v2 = iris$Sepal.Length,
    v12 = sample(iris$Species, nrow(iris)),
    v3 = iris$Petal.Length,
    v4 = factor(iris$Petal.Length < mean(iris$Petal.Length), labels = c(1, 2)),
    v42 = sample(iris$Sepal.Width, nrow(iris)),
    v5 = iris$Sepal.Width,
    v6 = sample(c(0, 1), nrow(iris), replace = TRUE),
    v7 = cut(
        x = iris$Sepal.Length,
        breaks = boxplot.stats(iris$Sepal.Length)$stats,
        include.lowest = TRUE
    )
)

# Create a better assign attribute
assign_better <- mm.column.variable.map(data = preds)

# Create a model matrix as used in mice
x <- obtain.design(preds)

# Save the contrast attributes
x_constrasts <- attributes(x)$contrasts

# Save the attributes
x_attr <- attributes(x)

# Remove the intercept
x <- x[, -1]

# Add useful attributes
attributes(x)$contrasts <- x_constrasts
attributes(x)$assign_better <- assign_better

# Reverse the model matrix
x_mm_reverse <- revert.factors(x)

# Nominal factor is recovered exactly
freq_table <- table(x_mm_reverse$v12, preds$v12)
testthat::expect_true(all(freq_table[upper.tri(freq_table)] == 0))

# Nominal factor with nested name is recovered exactly
freq_table <- table(x_mm_reverse$v1, preds$v1)
testthat::expect_true(all(freq_table[upper.tri(freq_table)] == 0))

# Binary factor with nested name is recovered exactly
freq_table <- table(x_mm_reverse$v4, preds$v4)
testthat::expect_true(all(freq_table[upper.tri(freq_table)] == 0))

# Binary factor are reconstructed as factors
testthat::expect_true(is.factor(x_mm_reverse$v4))

# Binary variables that are coded as numeric stay as such
testthat::expect_true(is.numeric(x_mm_reverse$v6))

# Ordinal factor is recovered factor
freq_table <- table(x_mm_reverse$v7, preds$v7)
testthat::expect_true(all(freq_table[upper.tri(freq_table)] == 0))
