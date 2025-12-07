# ==============================================================================
# Internal Helper Functions (Not exported)
# ==============================================================================

cov_my <- function(a, b) {
  t(a) %*% b
}

tr <- function(x) {
  sum(diag(x))
}

corr <- function(mat) {
  mat <- t(mat)
  mat <- mat - rowMeans(mat)
  # Standardize each variable
  mat <- mat / sqrt(rowSums(mat^2))
  # Calculate correlations
  cr <- tcrossprod(mat)
  return(cr)
}

cca_my <- function(a, b, r) {
  # Solving the generalized eigenvalue problem for CCA
  # Note: The formulation corresponds to finding canonical correlations
  mat <- solve(cov_my(a, a)) %*% cov_my(a, b) %*% solve(cov_my(b, b)) %*% cov_my(b, a)
  res <- eigen(mat)

  cor <- sqrt(res$values[1:r])
  xcoef <- res$vectors[, 1:r]
  xscore <- a %*% xcoef

  list(cor = cor, xcoef = xcoef, xscore = xscore)
}

# ==============================================================================
# Exported Functions and Methods
# ==============================================================================

#' Trace Ratio
#'
#' @description Evaluation of the estimated factors by trace ratios. The value is between 0 and 1; higher values indicate better estimation accuracy.
#'
#' @param G The true factors matrix.
#' @param Ghat The estimated factors matrix.
#'
#' @return A numeric value representing the trace ratio, defined as:
#' \eqn{\mathrm{TR} = \mathrm{tr} ( \mathbf{G}' \widehat{\mathbf{G}} (\widehat{\mathbf{G}}'\widehat{\mathbf{G}})^{-1} \widehat{\mathbf{G}}'\mathbf{G})/\mathrm{tr}(\mathbf{G'G})}.
#'
#' @export
#'
#' @examples
#' G <- matrix(rnorm(100 * 2), 100, 2)
#' Ghat <- G + matrix(rnorm(100 * 2, sd = 0.1), 100, 2)
#' TraceRatio(G, Ghat)
#'
TraceRatio <- function(G, Ghat) {
  if (ncol(G) > 0 & !all(is.na(Ghat))) {
    # TR calculation
    TR <- tr(t(G) %*% Ghat %*% solve(t(Ghat) %*% Ghat) %*% t(Ghat) %*% G) / tr(t(G) %*% G)
    return(TR)
  } else {
    return(0)
  }
}

#' Print Method for GFA Objects
#'
#' @description Print the summarized results of the estimated group factor model, such as the estimated global and local factor numbers and reference statistics.
#'
#' @param x The \code{GFA} object returned from the estimation algorithms (e.g., \code{APM}, \code{CCA}, \code{GCC}, \code{CP}).
#' @param ... Additional arguments passed to methods.
#'
#' @return No return value, called for side effects.
#' @export
#'
#' @method print GFA
print.GFA <- function(x, ...) {
  if (!inherits(x, "GFA")) {
    stop("Not a legitimate \"GFA\" object")
  }
  cat("The number of global factors is:", x$r0hat, "\n")
  if (!is.null(x$rhat)) {
    cat("The number of local factors are:", x$rhat, "\n")
  }
  if (!is.null(x$threshold)) {
    cat("The threshold is:", x$threshold, "\n")
  }
  # rho might be NULL in some methods or formatted differently, check before printing
  if (!is.null(x$rho)) {
    cat("The reference statistics are:", round(x$rho, 4), "\n")
  }
}

#' Print Method for GFD Objects
#'
#' @description Print the summary of the generated grouped data.
#'
#' @param x The \code{GFD} object returned from \code{\link{gendata}}.
#' @param ... Additional arguments passed to methods.
#'
#' @return No return value, called for side effects.
#' @export
#'
#' @method print GFD
print.GFD <- function(x, ...) {
  if (!inherits(x, "GFD")) {
    stop("Not a legitimate \"GFD\" object")
  }
  cat("Information of the data:", "\n")
  cat("T (Time points):", x[["T"]], "\n")
  cat("N (Variables per group):", x[["N"]], "\n")
  cat("M (Number of groups):", x[["M"]], "\n")
  cat("r0 (True global factors):", x[["r0"]], "\n")
  cat("r (True local factors):", x[["r"]], "\n")
  cat("Case:", x[["case"]], "\n")
}
