#' Aggregated Projection Method
#'
#' @description Aggregated Projection Method for Group Factor Model.
#'
#' @param y A list of the observation data, each element is a data matrix of each group with dimension \eqn{T \times N_m}.
#' @param rmax The maximum factor numbers of all groups. Default is 8.
#' @param r0 The number of global factors. Default is \code{NULL}, the algorithm will automatically estimate the number of global factors.
#'           If you have prior information about the true number of global factors, you can set it manually.
#' @param r The number of local factors in each group. Default is \code{NULL}, the algorithm will automatically estimate the number of local factors.
#'          If you have prior information, set it manually as an integer vector of length \eqn{M} (the number of groups).
#' @param localfactor Logical. If \code{FALSE} (default), local factors are not estimated. If \code{TRUE}, local factors will be estimated.
#' @param weight The weight of each projection matrix.
#'               If \code{TRUE} (default), weights are \eqn{w_m = N_m/N}.
#'               If \code{FALSE}, the mean of all projection matrices is calculated (equal weights).
#'               Can also be a numeric vector of length \eqn{M} specifying custom weights.
#' @param method The method used in the algorithm. Default is \code{"ic"}, can also be \code{"gap"}.
#' @param type The method used in estimating the factor numbers in each group initially. Default is \code{"IC3"}.
#'
#' @return An object of class \code{"GFA"} containing:
#' \item{r0hat}{The estimated number of global factors.}
#' \item{rhat}{The estimated number of local factors (if \code{localfactor = TRUE}).}
#' \item{rho}{The first \code{rmax} eigenvalues of the weighted projection matrix.}
#' \item{Ghat}{The estimated global factors.}
#' \item{loading_G}{A list consisting of the estimated global factor loadings.}
#' \item{Fhat}{The estimated local factors (if \code{localfactor = TRUE}).}
#' \item{loading_F}{A list consisting of the estimated local factor loadings (if \code{localfactor = TRUE}).}
#' \item{residual}{A list consisting of the residuals (if \code{localfactor = TRUE}).}
#' \item{threshold}{The threshold used in determining the number of global factors (only for \code{method = "ic"}).}
#'
#' @references
#' Aggregated Projection Method: A New Approach for Group Factor Model. Jiaqi Hu, Ting Li, Xueqin Wang (2025). Journal of the American Statistical Association, doi:10.1080/01621459.2025.2491154
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- GrFA::gendata()
#' APM(dat$y, rmax = 8, localfactor = TRUE, method = "ic")
#' APM(dat$y, rmax = 8, localfactor = TRUE, method = "gap")
#' }
APM <- function(y, rmax = 8, r0 = NULL, r = NULL, localfactor = FALSE, weight = TRUE, method = "ic", type = "IC3") {
  if (is.na(match(method, c("ic", "gap")))) {
    stop("invalid 'method' input")
  }
  M <- length(y)
  T <- nrow(y[[1]])
  Nm <- sapply(y, ncol)
  K <- list()
  Proj_Mat <- list()
  rhat <- rep(0, M)

  # Initial projection matrix estimation
  for (m in 1:M) {
    # Assuming est_num and FA are internal functions or imported from another package
    rhat[m] <- est_num(y[[m]], kmax = rmax, type = type)
    K[[m]] <- FA(y[[m]], r = rhat[m])$F
    if (rhat[m] == 0) {
      Proj_Mat[[m]] <- matrix(0, T, T)
    } else {
      Proj_Mat[[m]] <- K[[m]] %*% solve(t(K[[m]]) %*% K[[m]]) %*% t(K[[m]])
    }
  }

  # Weighting the projection matrices
  W_Proj_Mat <- matrix(0, T, T)
  if (isTRUE(weight)) {
    weights <- Nm / sum(Nm)
    for (m in 1:M) {
      W_Proj_Mat <- W_Proj_Mat + Proj_Mat[[m]] * weights[m]
    }
  } else if (isFALSE(weight)) {
    for (m in 1:M) {
      W_Proj_Mat <- W_Proj_Mat + Proj_Mat[[m]] / M
    }
  } else if (length(weight) == M && all(weight > 0)) {
    weight <- weight / sum(weight)
    for (m in 1:M) {
      W_Proj_Mat <- W_Proj_Mat + Proj_Mat[[m]] * weight[m]
    }
  } else {
    stop("invalid 'weight' input")
  }

  # Estimate r0 (Global Factor Number)
  eig_deco <- eigen(W_Proj_Mat)
  rho <- eig_deco$values[1:rmax]

  if (is.null(r0)) {
    if (method == "gap") {
      rho_gap_calc <- c(1 - 1 / sqrt(min(c(Nm, T))), rho)
      rho_gap <- rep(0, rmax)
      for (i in 1:rmax) {
        rho_gap[i] <- rho_gap_calc[i] - rho_gap_calc[i + 1]
      }
      r0hat <- which.max(rho_gap) - 1
      threshold <- NULL
      # rho stays as the original eigenvalues for return
    } else if (method == "ic") {
      # Re-calculate residuals for sigma estimation
      # Note: This uses the initial Proj_Mat, not W_Proj_Mat
      e_list <- list()
      for (m in 1:M) {
        e_list[[m]] <- y[[m]] - Proj_Mat[[m]] %*% y[[m]]
      }
      sigma2_e <- sum(sapply(e_list, function(x) sum(x^2))) / (T * sum(Nm))
      sigma2_y <- sum(sapply(y, function(x) sum(x^2))) / (T * sum(Nm))

      if (isTRUE(weight)) {
        Nmin <- mean(Nm) # Logic from original code kept, though usually Nmin implies min()
        Nmin <- min(Nm)  # Overwrites previous line immediately
      } else {
        Nmin <- min(Nm)
      }

      threshold <- 1 - (sqrt(Nmin) + sqrt(T)) / (sqrt(Nmin * T)) * log(log(sqrt(min(Nmin, T)))) * exp(sigma2_e / sigma2_y) * (1 - 1 / M)
      r0hat <- sum(rho >= threshold)
    }
  } else {
    threshold <- NULL
    r0hat <- r0
  }

  # Estimate G (Global Factors)
  if (r0hat > 0) {
    Ghat <- eig_deco$vectors[, 1:r0hat]
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

  # Estimate F (Local Factors) and Construct Result
  if (localfactor == FALSE) {
    res <- list(
      r0hat = r0hat,
      rho = rho,
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
      residual = e,
      threshold = threshold
    )
  }

  class(res) <- "GFA"
  return(res)
}
