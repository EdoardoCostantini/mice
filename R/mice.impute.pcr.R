#' Imputation by principal component regression
#'
#' Imputes univariate missing data by extracting PCs from all potential auxiliary variables.
#'
#' @aliases mice.impute.pcr pcr
#' @inheritParams mice.impute.norm.boot
#' @param npcs The number of principal components to extract for PC regression.
#' @return Vector with imputed data, same type as \code{y}, and of length
#' \code{sum(wy)}
#' @details
#' Extracts components from x and uses them as predictors in a regular call of
#' norm.boot().
#' @author Edoardo Costantini, 2022
#' @references
#'
#' Jolliffe, I. (2002) Principal Component Analysis (2nd ed). Springer.
#'
#' @family univariate imputation functions
#' @keywords datagen
#' @export
mice.impute.pcr <- function(y, ry, x, wy = NULL, npcs = 1L, ...) {

    # Set up -------------------------------------------------------------------

    if (is.null(wy)) wy <- !ry

    # Extract PCs --------------------------------------------------------------

    # PCA
    pcr_out <- stats::prcomp(x,
        center = TRUE,
        scale = TRUE
    )

    # Compute Explained Variance by each principal component
    pc_var_exp <- prop.table(pcr_out$sdev^2)

    # Keep PCs based on npcs object
    if (npcs >= 1) {
        # npcs as NOT-a-proportion
        pcs_keep <- 1:npcs
    } else {
        # npcs as a proportion
        pcs_keep <- cumsum(pc_var_exp) <= npcs
    }
    x_pcs <- pcr_out$x[, pcs_keep, drop = FALSE]
    pca_exp <- sum(pc_var_exp[pcs_keep])

    # Impute -------------------------------------------------------------------

    # Use traditional norm.boot machinery to obtain replacements
    imputes <- mice.impute.norm.boot(y = y, ry = ry, x = x_pcs, wy = NULL)

    # Return
    return(imputes)
}
