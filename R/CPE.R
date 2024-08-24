CP <- function(y, rmax = 8, r0 = NULL, r = NULL, localfactor = FALSE, type = "IC3") {
  M = length(y)
  T = nrow(y[[1]])
  Nm = sapply(y, ncol)
  K = list()
  Proj_Mat = list()
  for (m in 1:M) {
    K[[m]] = FA(y[[m]], r = rmax)$F
    Proj_Mat[[m]] = K[[m]] %*% solve(t(K[[m]]) %*% K[[m]]) %*% t(K[[m]])
  }
  Cir_Proj_Mat = diag(1, T, T)
  for (m in 1:M) {
    Cir_Proj_Mat = Cir_Proj_Mat %*% Proj_Mat[[m]]
  }
  Cir_Proj_Mat = t(Cir_Proj_Mat) %*% Cir_Proj_Mat
  eig_deco = eigen(Cir_Proj_Mat)
  eig_vec = eig_deco$vectors[, 1:rmax]
  rho = eig_deco$values[1:rmax]
  if (is.null(r0)) {
    ARSS = rep(0, rmax)
    for (i in 1:rmax) {
      for (m in 1:M) {
        ARSS[i] = ARSS[i] + t(eig_vec[, i]) %*% (diag(T) - Proj_Mat[[m]]) %*%
          eig_vec[, i]
      }
    }
    ARSS = ARSS/M
    nmin = min(c(T, Nm))
    logistic <- function(x) {
      P1 = 1
      P0 = 10^-3
      A = P1/P0 - 1
      tau = 14
      P1/(1 + A * exp(-tau * x))
    }
    r0hat = which.max(logistic(log(log(nmin)) * ARSS[2:rmax]) - logistic(log(log(nmin)) *
                                                                           ARSS[1:(rmax - 1)]))
  } else {
    if (!(r0%%1 == 0) | r0 <= 0) {
      stop("invalid 'r0' input")
    }
    r0hat = r0
  }
  Ghat = as.matrix(eig_deco$vectors[, 1:r0hat])*sqrt(T)
  Proj_G = Ghat %*% solve(t(Ghat) %*% Ghat) %*% t(Ghat)
  loading_G = list()
  for(m in 1:M){
    loading_G[[m]] = 1/T*t(y[[m]]) %*% Ghat
  }

  # estimate F
  if(localfactor == FALSE){
    res = list(r0hat = r0hat, rho = rho, Ghat = Ghat, loading_G = loading_G)
  }else{
    Fhat = list()
    loading_F = list()
    y_proj_G = lapply(y, function(x) x - Proj_G %*% x)
    if (is.null(r)) {
      rhat = rep(0, M)
      for (m in 1:M) {
        rhat[m] = est_num(y_proj_G[[m]], kmax = rmax - r0hat, type = type)
        fit = FA(y_proj_G[[m]], r = rhat[m])
        Fhat[[m]] = fit$F
        loading_F[[m]] = fit$L
      }
    } else {
      if (!(all(r%%1 == 0) && all(r >= 0))){
        stop("invalid r input")
      }
      rhat = r
      for (m in 1:M) {
        fit = FA(y_proj_G[[m]], r = rhat[m])
        Fhat[[m]] = fit$F
        loading_F[[m]] = fit$L
      }
    }
    # estimate e
    e = list()
    for(m in 1:M){
      if(rhat[m] > 0){
        e[[m]] = y_proj_G[[m]] - Fhat[[m]] %*% t(loading_F[[m]])
      }else{
        e[[m]] = y_proj_G[[m]]
      }
    }
    res = list(r0hat = r0hat, rhat = rhat, rho = rho, Ghat = Ghat, Fhat = Fhat,
               loading_G = loading_G, loading_F = loading_F, residual = e)
  }

  class(res) = "GFA"
  return(res)
}
