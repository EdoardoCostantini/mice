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
    x_df <- data.frame(x)

    # For every factor location re-create the original factor
    for(i in seq_along(mm_factor)){

        # Define columns containing information
        index_columns <- grepl(mm_factor[i], colnames(x_df))

        # Collapse values of the model matrix
        x_collapsed <- apply(
            x_df[, index_columns, drop = FALSE],
            1,
            paste0,
            collapse = ""
        )

        # Define number of categories
        ncat <- length(unique(x_collapsed))

        # What contrast is it?
        contr_type <- mm_contrasts[[mm_factor[i]]]

        # Evaluate the log-likelihood of new data under this model
        codes <- do.call(
            what = contr_type,
            args = list(
                n = ncat
            )
        )

        # Collapse codes to make them matchable
        code_collapse <- apply(codes, 1, paste0, collapse = "")

        # Create the factor
        x_factor <- factor(match(x_collapsed, code_collapse))

        # Drop the original model matrix columns
        x_df <- x_df[, !index_columns]

        # Add the reconstructed factor variable
        x_df <- cbind(x_df, x_factor)

        # Rename it
        colnames(x_df)[ncol(x_df)] <- paste0(mm_factor[i])
    }

    # Return result
    x_df
}
