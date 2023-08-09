# Test 1: Same pmm tests -------------------------------------------------------

# Define the test context
context("mice.impute.gspcr.pmm: same pmm tests")

# Set up data as in the regular pmm tests
xname <- c("age", "hgt", "wgt")
br <- boys[c(1:10, 101:110, 501:510, 601:620, 701:710), ]
r <- stats::complete.cases(br[, xname])
x <- br[r, xname]
y <- br[r, "tv"]
ry <- !is.na(y)

wy1 <- !ry
wy2 <- rep(TRUE, length(y))
wy3 <- rep(FALSE, length(y))
wy4 <- rep(c(TRUE, FALSE), times = c(1, length(y) - 1))

# Function output is as expected
test_that("Returns requested length", {
    expect_equal(length(mice.impute.gspcr.pmm(y, ry, x)), sum(!ry))
    expect_equal(length(mice.impute.gspcr.pmm(y, ry, x, wy = wy1)), sum(wy1))
    expect_equal(length(mice.impute.gspcr.pmm(y, ry, x, wy = wy2)), sum(wy2))
    expect_equal(length(mice.impute.gspcr.pmm(y, ry, x, wy = wy3)), sum(wy3))
    expect_equal(length(mice.impute.gspcr.pmm(y, ry, x, wy = wy4)), sum(wy4))
})

# Donor exclusion works
test_that("Excludes donors", {
    expect_false(all(c(15:25) %in% mice.impute.gspcr.pmm(y, ry, x, exclude = c(15:25))))
})

# Test 3: Works if bootstrap leads to constant variable ------------------------

# Define a constant X
x <- iris[1:50, -1]

# Define some y variable
y <- iris[1:50, 1]

# give it some NAs
y[1:10] <- NA
ry <- !is.na(y)

# Test returns requested length when constant variable is still in
expect_equal(length(mice.impute.gspcr.pmm(y, ry, x)), sum(!ry))

# Test 4: Works in mice function -----------------------------------------------

# Mice call works
mids_gspcr_pmm <- mice(
    data = cbind(fdd[, 8:12], fdd[, 5:7]),
    m = 3,
    maxit = 3,
    method = "gspcr.pmm",
    ridge = 0, 
    eps = 0, # bypasses remove.lindep()
    threshold = 1L,
    printFlag = FALSE,
    seed = 1234
)

# Check resulting object is a mids object
expect_equal(class(mids_gspcr_pmm), "mids")

# Test 5: Escape with same input if bootstrap causes constant ------------------

# Force constant y to mimic the bootstrap result
x_3 <- GSPCRexdata$X$cont[GSPCRexdata$y$ord == 3, ]
y_3 <- GSPCRexdata$y$ord[GSPCRexdata$y$ord == 3]
y_3[sample(1:length(y_3), length(y_3) * .1)] <- NA
ry_3 <- !is.na(y_3)

# Run imputation
imps <- mice.impute.gspcr.pmm(y = y_3, ry = ry_3, x = x_3)

# Test is the same as observed values of input imputations as expected
testthat::expect_equal(imps, y_3[ry_3])