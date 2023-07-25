#' Revert factors from dummy coding
#'
#' Given a model.matrix with contrast attribute, this function returns a data.frame with factors replacing the dummy coded binary and categorical predictors
#'
#' @param x A model.matrix with a contrast attribute
#' @details 
#' This function is used by any mice.impute.gspcr based imputation method to reconstruct categorical predictors as factors so that PCAmix can use them.
#' If the input \code{x} has no contrast attribute, it returns the input.
#' @return A version of the input \code{x} where dummy codes are replaced with factors.
#'
#' data.frame
#'  
#' @author Edoardo Costantini, 2023
#' @export
revert.factors <- function(x){
    # Extract the attributes of the matrix
    mm_attributes <- attributes(x)

    # Extract the attributes that related to contrasts
    mm_contrasts <- mm_attributes$contrasts

    # Identify the name of the factor predictors
    mm_factor <- names(mm_contrasts)

    # Make it a data.frame
    x <- data.frame(x)

    # For every factor location re-create the original factor
    for(i in seq_along(mm_factor)){
        # Define columns containing information
        index_columns <- grepl(mm_factor[i], colnames(x))

        # Collapse values of the model matrix
        x_collapsed <- apply(
            x[, index_columns, drop = FALSE],
            1,
            paste0,
            collapse = ""
        )

        # Create the factor
        x_factor <- factor(
            x_collapsed,
            labels = paste0("lvl_", seq_along(1:length(unique(x_collapsed))))
            )

        # Drop the original model matrix columns
        x <- x[, !index_columns]

        # Add the reconstructed factor variable
        x <- cbind(x, x_factor)

        # Rename it
        colnames(x)[ncol(x)] <- paste0(mm_factor[i])
    }

    # Return result
    x
}
