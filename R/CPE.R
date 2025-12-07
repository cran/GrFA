#' Circularly Projected Estimation
#'
#' @description Circularly Projected Estimation for Group Factor Model.
#'
#' @param y A list of the observation data, each element is a data matrix of each group with dimension \eqn{T \times N_m}.
#' @param rmax The maximum factor numbers of all groups. Default is 8.
#' @param r0 The number of global factors. Default is \code{NULL}, the algorithm will automatically estimate the number of global factors.
#'           If you have prior information about the true number of global factors, you can set it manually.
#' @param r The number of local factors in each group. Default is \code{NULL}, the algorithm will automatically estimate the number of local factors.
#'          If you have prior information, set it manually as an integer vector of length \eqn{M} (the number of groups).
#' @param localfactor Logical. If \code{FALSE} (default), local factors are not estimated. If \code{TRUE}, local factors will be estimated.
#' @param type The method used in estimating the local factor numbers in each group after projecting out the global factors. Default is \code{"IC3"}.
#'
#' @return An object of class \code{"GFA"} containing:
#' \item{r0hat}{The estimated number of global factors.}
#' \item{rhat}{The estimated number of local factors (if \code{localfactor = TRUE}).}
#' \item{rho}{The eigenvalues of the circular projection matrix.}
#' \item{Ghat}{The estimated global factors.}
#' \item{Fhat}{The estimated local factors (if \code{localfactor = TRUE}).}
#' \item{loading_G}{A list consisting of the estimated global factor loadings.}
#' \item{loading_F}{A list consisting of the estimated local factor loadings (if \code{localfactor = TRUE}).}
#' \item{residual}{A list consisting of the residuals (if \code{localfactor = TRUE}).}
#'
#' @references
#' Chen, M. (2023). Circularly Projected Common Factors for Grouped Data. Journal of Business & Economic Statistics, 41(2), 636-649.
#'
#' @export
#'
#' @examples
#' dat <- GrFA::gendata()
#' CP(dat$y, rmax = 8, localfactor = TRUE)
#'
CP <- function(y, rmax = 8, r0 = NULL, r = NULL, localfactor = FALSE, type = "IC3") {
  M <- length(y)
  T <- nrow(y[[1]])
  Nm <- sapply(y, ncol)
  K <- list()
  Proj_Mat <- list()

  # Initial Projection Matrices
  for (m in 1:M) {
    K[[m]] <- FA(y[[m]], r = rmax)$F
    Proj_Mat[[m]] <- K[[m]] %*% solve(t(K[[m]]) %*% K[[m]]) %*% t(K[[m]])
  }

  # Circular Projection Calculation
  Cir_Proj_Mat <- diag(1, T, T)
  for (m in 1:M) {
    Cir_Proj_Mat <- Cir_Proj_Mat %*% Proj_Mat[[m]]
  }
  Cir_Proj_Mat <- t(Cir_Proj_Mat) %*% Cir_Proj_Mat

  eig_deco <- eigen(Cir_Proj_Mat)
  eig_vec <- eig_deco$vectors[, 1:rmax]
  rho <- eig_deco$values[1:rmax]

  # Estimate r0
  if (is.null(r0)) {
    ARSS <- rep(0, rmax)
    for (i in 1:rmax) {
      for (m in 1:M) {
        ARSS[i] <- ARSS[i] + t(eig_vec[, i]) %*% (diag(T) - Proj_Mat[[m]]) %*% eig_vec[, i]
      }
    }
    ARSS <- ARSS / M
    nmin <- min(c(T, Nm))

    logistic <- function(x) {
      P1 <- 1
      P0 <- 10^-3
      A <- P1 / P0 - 1
      tau <- 14
      P1 / (1 + A * exp(-tau * x))
    }

    # Determine r0hat based on the gap in logistic transformed ARSS
    val_diff <- logistic(log(log(nmin)) * ARSS[2:rmax]) - logistic(log(log(nmin)) * ARSS[1:(rmax - 1)])
    r0hat <- which.max(val_diff)
  } else {
    if (!(r0 %% 1 == 0) | r0 <= 0) {
      stop("invalid 'r0' input")
    }
    r0hat <- r0
  }

  # Estimate Global Factors (G)
  Ghat <- as.matrix(eig_deco$vectors[, 1:r0hat]) * sqrt(T)
  Proj_G <- Ghat %*% solve(t(Ghat) %*% Ghat) %*% t(Ghat)
  loading_G <- list()
  for (m in 1:M) {
    loading_G[[m]] <- 1 / T * t(y[[m]]) %*% Ghat
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
        stop("invalid r input")
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
