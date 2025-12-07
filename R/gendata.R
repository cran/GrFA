#' Generate the grouped data
#'
#' @description Generate the grouped data for simulation studies.
#'
#' @param seed The seed used in \code{set.seed}. Default is 1.
#' @param T The number of time points. Default is 50.
#' @param N A vector representing the number of variables in each group. Default is \code{rep(30, 5)}.
#' @param r0 The number of global factors. Default is 2.
#' @param r A vector representing the number of the local factors. Notice, the length of \eqn{r} is the same as the length of \eqn{N} (which implies the number of groups \eqn{M}). Default is \code{rep(2, 5)}.
#' @param Phi_G Hyperparameter of the global factors (AR(1) coefficient). Default is 0.5. The value should be between 0 and 1.
#' @param Phi_F Hyperparameter of the local factors (AR(1) coefficient). Default is 0.5. The value should be between 0 and 1.
#' @param Phi_e Hyperparameter of the errors. Default is 0.5. The value should be between 0 and 1.
#' @param W_F Hyperparameter of the correlation of local factors. Only applicable when \code{case = 3}. The value should be between 0 and 1. Default is 0.5.
#' @param beta Hyperparameter of the errors (spatial correlation). Default is 0.2.
#' @param kappa Hyperparameter of signal to noise ratio. Default is 1.
#' @param case The case of the data-generating process. Default is 1. It can also be 2 or 3.
#'
#' @return An object of class \code{"GFD"} containing:
#' \item{y}{A list of the generated data matrices.}
#' \item{G}{The global factors matrix.}
#' \item{F}{A list of the local factors.}
#' \item{loading_G}{A list of the global factor loadings.}
#' \item{loading_F}{A list of the local factor loadings.}
#' \item{T}{The number of time points.}
#' \item{N}{The vector of variables per group.}
#' \item{M}{The number of groups.}
#' \item{r0}{The number of global factors.}
#' \item{r}{The vector of local factors.}
#' \item{case}{The generation case used.}
#'
#' @references
#' Aggregated Projection Method: A New Approach for Group Factor Model. Jiaqi Hu, Ting Li, Xueqin Wang (2025). Journal of the American Statistical Association, doi:10.1080/01621459.2025.2491154
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm
#' @export
#'
#' @examples
#' dat <- gendata()
#' print(dat)
gendata <- function(seed = 1, T = 50, N = rep(30, 5),
                    r0 = 2, r = rep(2, 5),
                    Phi_G = 0.5, Phi_F = 0.5, Phi_e = 0.5,
                    W_F = 0.5, beta = 0.2, kappa = 1, case = 1) {
  # Check the validation of input parameters
  if (r0 %% 1 != 0 | r0 < 0) stop("invalid r0 input")
  if (!(all(r %% 1 == 0) && all(r > 0))) stop("invalid r input")
  if (length(r) < 2) stop("length of r (number of groups) must be at least 2") # Fixed error message
  if (Phi_G < 0 | Phi_G >= 1) stop("invalid Phi_G input")
  if (Phi_F < 0 | Phi_F >= 1) stop("invalid Phi_F input")
  if (Phi_e < 0 | Phi_e >= 1) stop("invalid Phi_e input")
  if (W_F < 0 | W_F >= 1) stop("invalid W_F input")
  if (beta < 0 | beta >= 1) stop("invalid beta input")
  # if(kappa < 0) stop("invalid kappa input")
  if (!all(case %in% c(1, 2, 3))) stop("invalid case input")
  if (case == 2) {
    if (length(unique(r)) > 1) stop("invalid r input for case 2: all groups must have same r")
  }

  set.seed(seed)
  M <- length(r)
  if (length(kappa) == 1) kappa <- rep(kappa, M)

  F_list <- list() # Renamed to avoid confusion with F (FALSE)
  y <- list()
  loading_F <- list()
  loading_G <- list()
  e <- list()

  # Generate global factors
  G <- matrix(0, T, r0)
  G[1, ] <- rnorm(r0)
  for (t in 2:T) {
    G[t, ] <- Phi_G * G[t - 1, ] + rnorm(r0)
  }

  # Generate local factors
  if (case == 1) {
    for (m in 1:M) {
      Fm <- matrix(0, T, r[m])
      Fm[1, ] <- rnorm(r[m])
      for (t in 2:T) {
        Fm[t, ] <- Phi_F * Fm[t - 1, ] + rnorm(r[m])
      }
      F_list[[m]] <- Fm
    }
  } else if (case == 2) {
    F1 <- matrix(0, T, r[1])
    F1[1, ] <- rnorm(r[1])
    for (t in 2:T) {
      F1[t, ] <- Phi_F * F1[t - 1, ] + rnorm(r[1])
    }

    F2 <- matrix(0, T, r[1])
    F2[1, ] <- rnorm(r[1])
    for (t in 2:T) {
      F2[t, ] <- Phi_F * F2[t - 1, ] + rnorm(r[1])
    }

    for (m in 1:floor(M / 2)) {
      F_list[[m]] <- F1
    }
    for (m in (1 + floor(M / 2)):M) {
      F_list[[m]] <- F2
    }
  } else {
    # Case 3: Correlated local factors
    r_all <- sum(r)
    Omega <- matrix(W_F, r_all, r_all)
    diag(Omega) <- 1
    # Note: Requires 'mvtnorm' package
    F_all <- mvtnorm::rmvnorm(n = T, sigma = Omega)
    dim_1 <- c(1, 1 + cumsum(r))[1:M]
    dim_2 <- cumsum(r)
    for (m in 1:M) {
      F_list[[m]] <- F_all[, dim_1[m]:dim_2[m]]
    }
  }

  # Generate global and local factor loadings
  for (m in 1:M) {
    loading_F[[m]] <- matrix(rnorm(N[m] * r[m]), N[m], r[m])
    loading_G[[m]] <- matrix(rnorm(N[m] * r0), N[m], r0)
  }

  # Scale factors
  theta1 <- rep(0, M)
  theta2 <- rep(0, M)
  if (r0 == 0) {
    for (m in 1:M) {
      theta1[m] <- 1
      theta2[m] <- r[m] / (1 - Phi_G^2) / ((1 + 16 * beta^2) / (1 - Phi_e^2))
    }
  } else {
    for (m in 1:M) {
      theta1[m] <- r0 / (1 - Phi_G^2) / (r[m] / (1 - Phi_F^2))
      theta2[m] <- r0 / (1 - Phi_G^2) / ((1 + 16 * beta^2) / (1 - Phi_e^2))
    }
  }

  # Idiosyncratic error
  for (m in 1:M) {
    em <- matrix(0, T, N[m])
    em[1, ] <- rnorm(N[m])
    v <- matrix(rnorm(T * (N[m] + 2 * 8)), T, N[m] + 2 * 8)
    for (t in 2:T) {
      for (i in 1:N[m]) {
        em[t, i] <- Phi_e * em[t - 1, i] + v[t, i + 8] +
          beta * sum(v[t, i:(i + 8 - 1)]) +
          beta * sum(v[t, (i + 8 + 1):(i + 2 * 8)])
      }
    }
    e[[m]] <- em
  }

  # Construct Y
  for (m in 1:M) {
    y[[m]] <- G %*% t(loading_G[[m]]) +
      sqrt(theta1[m]) * F_list[[m]] %*% t(loading_F[[m]]) +
      sqrt(kappa[m] * theta2[m]) * e[[m]]
  }

  res <- list(
    y = y,
    G = G,
    F = F_list,
    loading_G = loading_G,
    loading_F = loading_F,
    T = T,
    N = N,
    M = M,
    r0 = r0,
    r = r,
    case = case
  )
  class(res) <- "GFD"
  return(res)
}
