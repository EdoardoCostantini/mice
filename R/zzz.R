#' Echoes the package version number
#'
#' @param pkg A character vector with the package name.
#' @return A character vector containing the package name, version number and
#' installed directory.
#' @author Stef van Buuren, Oct 2010
#' @keywords misc
#' @examples
#' version()
#' version("base")
#' @export
version <- function(pkg = "mice") {
  lib <- dirname(system.file(package = pkg))
  d <- packageDescription(pkg)
  return(paste(d$Package, d$Version, d$Date, lib))
}

.onLoad <- function(libname, pkgname) {
  warning(
    paste0(
      version(), " was loaded.\n\n",
      "This is an experimental version of the `mice` R package
developed by Edoardo Costantini as an unofficial fork to the 
main package. You can find the code base by following this link:

https://github.com/EdoardoCostantini/mice/tree/develop-gspcr

I created this version to demonstrate the potential use of 
generalised supervised principal components as a univariate
imputation method in `mice`.\n\n",
      "IMPORTANT!",
      "\n\n",
      "I created this version in autonomy from maintainers and 
other contributors of the CRAN version of the `mice` R package.
If you encounter any issues using it please contact me at: 

e.costantini@tilburguniversity.edu."
    ),
    call. = FALSE
  )
}