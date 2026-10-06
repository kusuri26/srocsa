#' @title Print the summary of 'srocsa'
#' @description Print the summary of 'srocsa'
#' @param x
#' An object returned by the function \code{summary.srocsa}.
#' @param digits
#' The digits of the results
#' @param ...
#' Other arguments in function \code{\link{print.default}}
#' @method print summary.srocsa
#' @return Print the summary of a srocsa object
#' @export
#'
print.summary.srocsa <- function(x, digits = 3, ...) {

  if (!is.null(x$call)) {
    cat("Input:\n")
    print(x$call)
    cat("\n")
  }

  ## main results ----
  for (b in x$blocks) .srocsa.print.block(b, digits = digits)

  ## all the elements of the fitted models ----
  cat("==================================================\n")
  cat("All elements of the fitted models\n")
  cat("==================================================\n\n")

  for (b in x$blocks) {
    for (i in seq_along(b$fits)) {

      cat("--------------------------------------------------\n")
      cat(if (!is.null(b$time)) paste0("time = ", b$time, ", "), "p = ", b$p[i], "\n", sep = "")
      cat("--------------------------------------------------\n")

      fit <- b$fits[[i]]

      if (!is.list(fit)) {
        cat("Estimation failed\n\n")
        next
      }

      for (nm in names(fit)) {
        cat("$", nm, "\n", sep = "")
        v <- fit[[nm]]
        print(if (is.numeric(v)) round(v, digits) else v, ...)
        cat("\n")
      }
    }
  }

  invisible(x)

}
