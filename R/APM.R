APM <- function(y, rmax = 8, r0 = NULL, r = NULL, weight = TRUE, method = "ic1", type = "BIC3") {
  if (is.na(match(method, c("ic1", "ic2", "gap")))) {
    stop("invalid 'method' input")
  }
  M = length(y)
  T = nrow(y[[1]])
  Nm = sapply(y, ncol)
  K = list()
  Proj_Mat = list()
  for (m in 1:M) {
    K[[m]] = FA(y[[m]], r = rmax)$F
    Proj_Mat[[m]] = K[[m]] %*% solve(t(K[[m]]) %*% K[[m]]) %*% t(K[[m]])
  }
  W_Proj_Mat = matrix(0, T, T)
  if (weight == TRUE) {
    weights = Nm/sum(Nm)
    for (m in 1:M) {
      W_Proj_Mat = W_Proj_Mat + Proj_Mat[[m]] * weights[m]
    }
  } else if (weight == FALSE) {
    for (m in 1:M) {
      W_Proj_Mat = W_Proj_Mat + Proj_Mat[[m]]/M
    }
  } else if (length(weight) == M & all(weight > 0)) {
    weight = weight/sum(weight)
    W_Proj_Mat = W_Proj_Mat + Proj_Mat[[m]] * weight[m]
  } else {
    stop("invalid 'weight' input")
  }
  eig_deco = eigen(W_Proj_Mat)
  rho = eig_deco$values[1:rmax]
  if (is.null(r0)) {
    if (method == "gap") {
      rho = c(1 - 1/sqrt(min(c(Nm, T))), rho)
      rho_gap = rep(0, rmax)
      for (i in 1:rmax) {
        rho_gap[i] = rho[i] - rho[i + 1]
      }
      r0hat = which.max(rho_gap) - 1
      threshold = NULL
    } else {
      rho = c(1, rho)
      rhat = rep(0, M)
      K = list()
      for (m in 1:M) {
        rhat[m] = est_num(y[[m]], kmax = rmax, type = type)
        K[[m]] = FA(y[[m]], r = rhat[m])$F
      }
      sigma2_y = sum(sapply(y, function(x) sum(x^2)))/(T * sum(Nm))
      f = function(y, F) {
        e = y - F %*% solve(t(F) %*% F) %*% t(F) %*% y
        sum(e^2)
      }
      sigma2_e = sum(mapply(f, y, K)/(T * sum(Nm)))
      Nmin = min(Nm)
      if (method == "ic1") {
        threshold = 1 - (Nmin * T)^(-1/4) * log(log(Nmin * T/(Nmin + T))) * exp(sigma2_e/sigma2_y) * (1 - 1/M)
      } else {
        threshold = 1 - (sqrt(Nmin) + sqrt(T))/(sqrt(Nmin * T)) * log(log(sqrt(min(Nmin, T)))) * exp(sigma2_e/sigma2_y) * (1 - 1/M)
      }
      r0hat = sum(rho >= threshold) - 1
    }
  } else {
    if (!(r0%%1 == 0) | r0 < 0) {
      stop("invalid 'r0' input")
    }
    threshold = NULL
    r0hat = r0
    rho = c(1, rho)
  }
  if (r0hat > 0) {
    Ghat = eig_deco$vectors[, 1:r0hat]
    Ghat = sqrt(T) * as.matrix(Ghat)
    Proj_G = Ghat %*% solve(t(Ghat) %*% Ghat) %*% t(Ghat)
  } else {
    Ghat = NULL
    Proj_G = matrix(0, T, T)
  }
  y_proj_G = lapply(y, function(x) x - Proj_G %*% x)
  Fhat = list()
  if (is.null(r)) {
    rhat = rep(0, M)
    for (m in 1:M) {
      rhat[m] = est_num(y_proj_G[[m]], kmax = rmax - r0hat, type = type)
      Fhat[[m]] = FA(y_proj_G[[m]], r = rhat[m])$F
    }
  } else {
    if (!(all(r%%1 == 0) && all(r >= 0))) {
      stop("invalid 'r' input")
    }
    rhat = r
    for (m in 1:M) {
      Fhat[[m]] = FA(y_proj_G[[m]], r[m])$F
    }
  }
  loading_G = list()
  loading_F = list()
  e = list()
  for(m in 1:M){
    if(r0hat == 0 & rhat[m] == 0){
      loading_F[[m]] = NA
      loading_G[[m]] = NA
      e[[m]] = y[[m]]
    }else if(r0hat == 0 & rhat[m] > 0){
      loading_G[[m]] = NA
      loading_F[[m]] = 1/T*t(y[[m]]) %*% Fhat[[m]]
      e[[m]] = y[[m]] - Fhat[[m]] %*% t(loading_F[[m]])
    }else if(r0hat > 0 & rhat[m] == 0){
      loading_G[[m]] = 1/T*t(y[[m]]) %*% Ghat
      loading_F[[m]] = NA
      e[[m]] = y[[m]] - Ghat %*% t(loading_G[[m]])
    }else{
      loading_G[[m]] = 1/T*t(y[[m]]) %*% Ghat
      loading_F[[m]] = 1/T*t(y[[m]]) %*% Fhat[[m]]
      e[[m]] = y[[m]] - Ghat %*% t(loading_G[[m]]) - Fhat[[m]] %*% t(loading_F[[m]])
    }
  }
  res = list(r0hat = r0hat, rhat = rhat, rho = rho[-1], Ghat = Ghat, Fhat = Fhat,
             loading_G = loading_G, loading_F = loading_F, residual = e, threshold = threshold)
  class(res) = "GFA"
  return(res)
}
