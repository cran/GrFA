#' Factor Analysis
#'
#' @description Performs Factor Analysis using Principal Component Analysis (PCA) to extract factors and loadings.
#'
#' @param X The observation data matrix of dimension \eqn{T \times N}.
#' @param r The number of factors to estimate.
#'
#' @return A list containing:
#' \item{F}{The estimated factors matrix of dimension \eqn{T \times r}.}
#' \item{L}{The estimated factor loadings matrix of dimension \eqn{N \times r}.}
#'
#' @references
#' Bai, J., & Ng, S. (2002). Determining the number of factors in approximate factor models. Econometrica, 70(1), 191-221.
#'
#' @author Jiaqi Hu
#' @export
#'
#' @examples
#' X <- matrix(rnorm(100*20), 100, 20)
#' res <- FA(X, r = 2)
#' head(res$F)
#' head(res$L)
#'
FA <- function(X, r) {
  T <- nrow(X)
  N <- ncol(X)

  if (r == 0) {
    return(list(F = NA, L = NA))
  }

  if (T > N) {
    S <- t(X) %*% X / (N * T)
    eig <- eigen(S)
    e <- eig$vectors
    L <- e[, 1:r] * sqrt(N)
    L <- as.matrix(L)
    F <- X %*% L / N
    F <- as.matrix(F)
    # Normalize factors
    F <- F %*% diag(1 / sqrt(diag(t(F) %*% F / T)), nrow = r, ncol = r)
    L <- t(X) %*% F / T
    L <- as.matrix(L)
  } else {
    S <- X %*% t(X) / (N * T)
    eig <- eigen(S)
    e <- eig$vectors
    F <- e[, 1:r] * sqrt(T)
    F <- as.matrix(F)
    L <- t(X) %*% F / T
    L <- as.matrix(L)
  }
  return(list(F = F, L = L))
}


#' Estimate Factor Numbers
#'
#' @description Estimates the number of factors using various Information Criteria (IC) and Eigenvalue Ratio tests.
#'
#' @param X The observation data matrix of dimension \eqn{T \times N}.
#' @param kmax The maximum number of factors to consider. Default is 8.
#' @param type The criterion used in determining the number of factors. Default is \code{"BIC3"}.
#'             Options: \code{"PC1"}, \code{"PC2"}, \code{"PC3"}, \code{"IC1"}, \code{"IC2"}, \code{"IC3"},
#'             \code{"AIC3"}, \code{"BIC3"}, \code{"ER"}, \code{"GR"}.
#'
#' @return \item{rhat}{The estimated number of factors (an integer).}
#'
#' @references
#' Bai, J., & Ng, S. (2002). Determining the number of factors in approximate factor models. Econometrica, 70(1), 191-221.
#'
#' Ahn, S. C., & Horenstein, A. R. (2013). Eigenvalue ratio test for the number of factors. Econometrica, 81(3), 1203-1227.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100*20), 100, 20)
#' est_num(X, kmax = 8, type = "BIC3")
#' est_num(X, kmax = 8, type = "ER")
#' }
est_num <- function(X, kmax = 8, type = "BIC3") {
  valid_types <- c("PC1", "PC2", "PC3", "IC1", "IC2", "IC3", "AIC3", "BIC3", "ER", "GR")
  if (is.na(match(type, valid_types))) {
    stop("Invalid input 'type'. Must be one of: ", paste(valid_types, collapse = ", "))
  }
  if (kmax == 0) {
    return(0)
  }

  T <- nrow(X)
  N <- ncol(X)

  # Efficiently compute eigenvalues based on dimensions
  if (T > N) {
    eig_val <- eigen(t(X) %*% X, only.values = TRUE)$values
  } else {
    eig_val <- eigen(X %*% t(X), only.values = TRUE)$values
  }

  # Calculate residual variance V(k)
  V <- rep(0, kmax + 1)
  V[1] <- sum(eig_val)
  for (i in 1:kmax) {
    V[i + 1] <- V[i] - eig_val[i]
  }
  V <- V / (N * T)

  CNT <- min(sqrt(N), sqrt(T))

  # Bai and Ng (2002) Criteria
  PC1 <- V + (0:kmax) * V[kmax + 1] * (N + T) / (N * T) * log(N * T / (N + T))
  PC2 <- V + (0:kmax) * V[kmax + 1] * (N + T) / (N * T) * log(CNT^2)
  PC3 <- V + (0:kmax) * V[kmax + 1] * log(CNT^2) / CNT^2
  IC1 <- log(V) + (0:kmax) * (N + T) / (N * T) * log(N * T / (N + T))
  IC2 <- log(V) + (0:kmax) * (N + T) / (N * T) * log(CNT^2)
  IC3 <- log(V) + (0:kmax) * log(CNT^2) / CNT^2
  AIC3 <- V + (0:kmax) * V[kmax + 1] * 2 * (N + T - (0:kmax)) / (N * T)
  BIC3 <- V + (0:kmax) * V[kmax + 1] * (N + T - (0:kmax)) * log(N * T) / (N * T)

  # Ahn and Horenstein (2013) Eigenvalue Ratio (ER)
  ER <- eig_val[1:kmax] / eig_val[2:(kmax + 1)]
  ER_k <- which.max(ER)

  # Ahn and Horenstein (2013) Growth Ratio (GR)
  V_GR <- rep(0, kmax + 2)
  V_GR[1] <- sum(eig_val)
  for (i in 1:(kmax + 1)) {
    V_GR[i + 1] <- V_GR[i] - eig_val[i]
  }
  V_GR <- V_GR / (N * T)

  GR <- rep(0, kmax)
  for (k in 1:kmax) {
    GR[k] <- log(V_GR[k] / V_GR[k + 1]) / log(V_GR[k + 1] / V_GR[k + 2])
  }
  GR_k <- which.max(GR)

  opt <- c(
    which.min(PC1) - 1, which.min(PC2) - 1, which.min(PC3) - 1,
    which.min(IC1) - 1, which.min(IC2) - 1, which.min(IC3) - 1,
    which.min(AIC3) - 1, which.min(BIC3) - 1, ER_k, GR_k
  )
  names(opt) <- valid_types

  return(opt[[type]])
}
