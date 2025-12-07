#' Canonical Correlation Estimation
#'
#' @description Canonical Correlation Estimation for Group Factor Model.
#'
#' @param y A list of the observation data, each element is a data matrix of each group with dimension \eqn{T \times N_m}.
#' @param rmax The maximum factor numbers of all groups. Default is 8.
#' @param r0 The number of global factors. Default is \code{NULL}, the algorithm will automatically estimate the number of global factors.
#'           If you have prior information about the true number of global factors, you can set it manually.
#' @param r The number of local factors in each group. Default is \code{NULL}, the algorithm will automatically estimate the number of local factors.
#'          If you have prior information, set it manually as an integer vector of length \eqn{M} (the number of groups).
#' @param localfactor Logical. If \code{FALSE} (default), local factors are not estimated. If \code{TRUE}, local factors will be estimated.
#' @param method The method used in the algorithm. Default is \code{"CCD"}, can also be \code{"MCC"}.
#' @param type The method used in estimating the factor numbers in each group initially. Default is \code{"IC3"}.
#'
#' @return An object of class \code{"GFA"} containing:
#' \item{r0hat}{The estimated number of global factors.}
#' \item{rhat}{The estimated number of local factors (if \code{localfactor = TRUE}).}
#' \item{rho}{The vector of average canonical correlations (eigenvalues).}
#' \item{Ghat}{The estimated global factors.}
#' \item{Fhat}{The estimated local factors (if \code{localfactor = TRUE}).}
#' \item{loading_G}{A list consisting of the estimated global factor loadings.}
#' \item{loading_F}{A list consisting of the estimated local factor loadings (if \code{localfactor = TRUE}).}
#' \item{residual}{A list consisting of the residuals (if \code{localfactor = TRUE}).}
#' \item{threshold}{The threshold used in determining the number of global factors (only for \code{method = "MCC"}).}
#'
#' @references
#' Choi, I., Lin, R., & Shin, Y. (2021). Canonical correlation-based model selection for the multilevel factors. Journal of Econometrics.
#'
#' @export
#'
#' @examples
#' dat <- GrFA::gendata()
#' CCA(dat$y, rmax = 8, localfactor = TRUE, method = "CCD")
#' CCA(dat$y, rmax = 8, localfactor = TRUE, method = "MCC")
#'
CCA <- function(y, rmax = 8, r0 = NULL, r = NULL, localfactor = FALSE, method = "CCD", type = "IC3") {
  if (is.na(match(method, c("CCD", "MCC")))) {
    stop("invalid 'method' input")
  }
  M <- length(y)
  T <- nrow(y[[1]])
  Nm <- sapply(y, ncol)
  rhat <- rep(0, M)

  # Initial estimation of factor numbers
  for (m in 1:M) {
    # Assuming est_num is an internal function
    rhat[m] <- est_num(y[[m]], kmax = rmax, type = "BIC3")
  }
  rmaxstar <- max(rhat)

  K <- list()
  for (m in 1:M) {
    K[[m]] <- FA(y[[m]], r = rmaxstar)$F
  }

  lmh <- array(0, dim = c(M, M, rmaxstar))
  for (i in 1:(M - 1)) {
    for (j in (i + 1):M) {
      a <- K[[i]]
      b <- K[[j]]
      lmh[i, j, ] <- eigen(solve(cov_my(a, a)) %*% cov_my(a, b) %*% solve(cov_my(b, b)) %*% cov_my(b, a))$values
    }
  }

  rho <- apply(lmh, 3, sum) * 2 / (M * (M - 1))
  rho <- c(1, rho)

  # Estimate r0
  if (is.null(r0)) {
    if (method == "CCD") {
      CCD <- rep(0, rmaxstar)
      for (i in 1:rmaxstar) {
        CCD[i] <- rho[i] - rho[i + 1]
      }
      r0hat <- which.max(CCD) - 1
      threshold <- NULL
    } else {
      # MCC Method
      Nmin <- min(Nm)
      PNT <- (log(Nmin) + log(T)) / sqrt(Nmin * T) * log(log(Nmin * T))
      sigma2_y <- sum(sapply(y, function(x) sum(x^2))) / (T * sum(Nm))
      f <- function(y, F) {
        e <- y - F %*% solve(t(F) %*% F) %*% t(F) %*% y
        sum(e^2)
      }
      sigma2_e <- sum(mapply(f, y, K) / (T * sum(Nm)))
      threshold <- 1 - exp(sigma2_e / sigma2_y) * PNT
      r0hat <- sum(rho > threshold) - 1
    }
  } else {
    if (!(r0 %% 1 == 0) | r0 < 0) {
      stop("invalid 'r0' input")
    }
    r0hat <- r0
    threshold <- NULL
  }

  # Estimate G
  if (r0hat > 0) {
    max_idx <- which.max(lmh[, , 1])
    arr_idx <- arrayInd(max_idx, .dim = c(M, M))
    i <- arr_idx[1]
    j <- arr_idx[2]

    Ghat <- cca_my(K[[i]], K[[j]], r0hat)$xscore
    Proj_G <- Ghat %*% solve(t(Ghat) %*% Ghat) %*% t(Ghat)
    y_proj <- lapply(y, function(x) x - Proj_G %*% x)

    Fhat <- list()
    y_proj_F <- list()
    rhat <- rep(0, M)
    for (m in 1:M) {
      rhat[m] <- est_num(y_proj[[m]], kmax = rmaxstar - r0hat, type = type)
      if (rhat[m] == 0) {
        Fhat[[m]] <- matrix(0, T, 0)
        y_proj_F[[m]] <- y[[m]]
      } else {
        Fhat[[m]] <- FA(y_proj[[m]], rhat[m])$F
        y_proj_F[[m]] <- y[[m]] - Fhat[[m]] %*% solve(t(Fhat[[m]]) %*% Fhat[[m]]) %*% t(Fhat[[m]]) %*% y[[m]]
      }
    }
    y_proj_F_all <- do.call("cbind", y_proj_F)
    Ghat <- FA(y_proj_F_all, r0hat)$F
    Proj_G <- Ghat %*% solve(t(Ghat) %*% Ghat) %*% t(Ghat)
    loading_G <- list()
    for (m in 1:M) {
      loading_G[[m]] <- 1 / T * t(y[[m]]) %*% Ghat
    }
  } else {
    Ghat <- NA
    Proj_G <- diag(0, T, T)
    loading_G <- NA
  }

  # Estimate F (Local Factors)
  if (localfactor == FALSE) {
    res <- list(
      r0hat = r0hat,
      rho = rho[-1],
      Ghat = Ghat,
      loading_G = loading_G,
      threshold = threshold
    )
  } else {
    Fhat <- list()
    loading_F <- list()
    y_proj_G <- lapply(y, function(x) x - Proj_G %*% x)

    if (is.null(r)) {
      rhat <- rep(0, M)
      for (m in 1:M) {
        rhat[m] <- est_num(y_proj_G[[m]], kmax = rmaxstar - r0hat, type = type)
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
      rho = rho[-1],
      Ghat = Ghat,
      Fhat = Fhat,
      loading_G = loading_G,
      loading_F = loading_F,
      residual = e,
      threshold = threshold
    )
  }

  class(res) <- "GFA"
  return(res)
}
