#' Generalised Canonical Correlation
#'
#' @description Generalised Canonical Correlation Estimation for Group Factor Model.
#'
#' @param y A list of the observation data, each element is a data matrix of each group with dimension \eqn{T \times N_m}.
#' @param rmax The maximum factor numbers of all groups. Default is 8.
#' @param r0 The number of global factors. Default is \code{NULL}, the algorithm will automatically estimate the number of global factors.
#'           If you have prior information about the true number of global factors, you can set it manually.
#' @param r The number of local factors in each group. Default is \code{NULL}, the algorithm will automatically estimate the number of local factors.
#'          If you have prior information, set it manually as an integer vector of length \eqn{M} (the number of groups).
#' @param localfactor Logical. If \code{FALSE} (default), local factors are not estimated. If \code{TRUE}, local factors will be estimated.
#' @param type The method used in estimating the factor numbers in each group initially. Default is \code{"IC3"}.
#'
#' @return An object of class \code{"GFA"} containing:
#' \item{r0hat}{The estimated number of global factors.}
#' \item{rhat}{The estimated number of local factors (if \code{localfactor = TRUE}).}
#' \item{rho}{The ratio of the singular values used to estimate the number of global factors.}
#' \item{Ghat}{The estimated global factors.}
#' \item{Fhat}{The estimated local factors (if \code{localfactor = TRUE}).}
#' \item{loading_G}{A list consisting of the estimated global factor loadings.}
#' \item{loading_F}{A list consisting of the estimated local factor loadings (if \code{localfactor = TRUE}).}
#' \item{residual}{A list consisting of the residuals (if \code{localfactor = TRUE}).}
#'
#' @references
#' Lin, R., & Shin, Y. (2023). Generalised Canonical Correlation Estimation of the Multilevel Factor Model. Available at SSRN 4295429.
#'
#' @export
#'
#' @examples
#' dat <- GrFA::gendata()
#' GCC(dat$y, rmax = 8, localfactor = TRUE)
#'
GCC <- function(y, rmax = 8, r0 = NULL, r = NULL, localfactor = FALSE, type = "IC3") {
  M <- length(y)
  T <- nrow(y[[1]])
  Nm <- sapply(y, ncol)
  CNT <- sqrt(min(c(T, Nm)))
  K <- list()

  # Initial Factor Extraction
  for (m in 1:M) {
    # Assuming FA is an internal function
    K[[m]] <- FA(y[[m]], r = rmax)$F
  }

  # Construct the Phi matrix
  Phi <- matrix(0, T * M * (M - 1) / 2, M * rmax)
  row <- 1
  for (i in 1:(M - 1)) {
    for (j in (i + 1):M) {
      Phi[row:(row + T - 1), ((i - 1) * rmax + 1):(i * rmax)] <- K[[i]]
      Phi[row:(row + T - 1), ((j - 1) * rmax + 1):(j * rmax)] <- -K[[j]]
      row <- row + T
    }
  }

  # SVD and Rho Calculation
  svd_phi <- svd(Phi)
  delta2 <- sort(svd_phi$d^2)
  delta2 <- c(sum(delta2) / (CNT * M * rmax), delta2)
  rho <- delta2[2:(rmax + 1)] / delta2[1:rmax]

  # Estimate r0
  if (is.null(r0)) {
    r0hat <- which.max(rho) - 1
  } else {
    if (!(r0 %% 1 == 0) | r0 < 0) {
      stop("invalid 'r0' input")
    }
    r0hat <- r0
  }

  # Estimate Global Factors (G)
  if (r0hat > 0) {
    psi <- matrix(0, T, 0)
    Q <- svd_phi$v[, 1:rmax]
    for (m in 1:M) {
      psi <- cbind(psi, K[[m]] %*% Q[((m - 1) * rmax + 1):(m * rmax), ])
    }
    Ghat <- eigen(psi %*% t(psi))$vectors[, 1:r0hat]
    Ghat <- sqrt(T) * as.matrix(Ghat)
    Proj_G <- Ghat %*% solve(t(Ghat) %*% Ghat) %*% t(Ghat)
    loading_G <- list()
    for (m in 1:M) {
      loading_G[[m]] <- 1 / T * t(y[[m]]) %*% Ghat
    }
  } else {
    Ghat <- NA
    Proj_G <- matrix(0, T, T)
    loading_G <- NA
  }

  # Estimate Local Factors (F) and Residuals
  if (localfactor == FALSE) {
    res <- list(
      r0hat = r0hat,
      rho = rho,
      Ghat = Ghat,
      loading_G = loading_G
    )
  } else {
    Fhat <- list()
    loading_F <- list()
    y_proj_G <- lapply(y, function(x) x - Proj_G %*% x)

    if (is.null(r)) {
      rhat <- rep(0, M)
      for (m in 1:M) {
        rhat[m] <- est_num(y_proj_G[[m]], kmax = rmax - r0hat, type = type)
        fit <- FA(y_proj_G[[m]], r = rhat[m])
        Fhat[[m]] <- fit$F
        loading_F[[m]] <- fit$L
      }
    } else {
      if (!(all(r %% 1 == 0) && all(r >= 0))) {
        stop("invalid 'r' input")
      }
      rhat <- r
      for (m in 1:M) {
        fit <- FA(y_proj_G[[m]], r = rhat[m])
        Fhat[[m]] <- fit$F
        loading_F[[m]] <- fit$L
      }
    }

    # Estimate residuals (e)
    e <- list()
    for (m in 1:M) {
      if (rhat[m] > 0) {
        e[[m]] <- y_proj_G[[m]] - Fhat[[m]] %*% t(loading_F[[m]])
      } else {
        e[[m]] <- y_proj_G[[m]]
      }
    }

    res <- list(
      r0hat = r0hat,
      rhat = rhat,
      rho = rho,
      Ghat = Ghat,
      Fhat = Fhat,
      loading_G = loading_G,
      loading_F = loading_F,
      residual = e
    )
  }

  class(res) <- "GFA"
  return(res)
}
