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

    # Original data starting point
    x_original <- data.frame(x, check.names = FALSE)

    # For every factor location re-create the original factor
    for(i in seq_along(mm_factor)){

        # Find column membership
        index_columns <- mm_attributes$assign_better %in% mm_factor[i]

        # Collapse values of the model matrix
        x_collapsed <- apply(
            x[, index_columns, drop = FALSE],
            1,
            paste0,
            collapse = ""
        )

        # Define number of categories
        ncat <- sum(index_columns) + 1

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
        x_factor <- factor(
            match(x_collapsed, code_collapse)
        )

        # Drop the original model matrix columns
        to_drop <- which(colnames(x_original) %in% colnames(x)[index_columns])
        x_original <- x_original[, -to_drop]

        # Add the reconstructed factor variable
        x_original <- cbind(x_original, x_factor)

        # Rename it
        colnames(x_original)[ncol(x_original)] <- paste0(mm_factor[i])
    }

    # Return result
    x_original
}

#' Define a column-variable map for the model.matrix
#'
#' Given a data.frame for which a model.matrix should be produced, this function returns a vector clearly defining which columns of the model.matrix represent which variable of the original data.frame
#'
#' @param data A data.frame
#' @details
#' This function is used to create an attribute object that the revert.factor function can use to reconstruct the original data.frame from a model.matrix.
#' @return A character vector indicating to which variable of the original data.frame the model.matrix columns belong
#'
#' character vector
#'
#' @author Edoardo Costantini, 2023
#' @export
mm.column.variable.map <- function(data) {

    # Define an empty vector to store column-variable assignments
    assign_names <- NULL

    # Loop through the variables to populate the column-variable map
    for (j in 1:ncol(data)) {
        # Store the variable for easy processing
        v <- data[, j]

        # Store the variable name
        varname <- colnames(data)[j]

        if (is.factor(v)) {
            assign_names <- c(assign_names, rep(varname, nlevels(v) - 1))
        } else {
            assign_names <- c(assign_names, varname)
        }
    }

    return(assign_names)
}
