APM <- function(y, rmax = 8, r0 = NULL, r = NULL, localfactor = FALSE, weight = TRUE, method = "ic", type = "IC3") {
  if (is.na(match(method, c("ic", "gap")))) {
    stop("invalid 'method' input")
  }
  M = length(y)
  T = nrow(y[[1]])
  Nm = sapply(y, ncol)
  K = list()
  Proj_Mat = list()
  rhat = rep(0, M)
  for (m in 1:M) {
    rhat[m] = est_num(y[[m]], kmax = rmax, type = type)
    K[[m]] = FA(y[[m]], r = rhat[m])$F
    if(rhat[m] == 0){
      Proj_Mat[[m]] = matrix(0, T, T)
    } else{
      Proj_Mat[[m]] = K[[m]] %*% solve(t(K[[m]]) %*% K[[m]]) %*% t(K[[m]])
    }
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

  # estimate r0
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
      rho = rho[-1]
    } else if (method == "ic"){
      rhat = rep(0, M)
      K = list()
      e = list()
      for (m in 1:M) {
        e[[m]] = y[[m]] - Proj_Mat[[m]] %*% y[[m]]
      }
      sigma2_e = sum(sapply(e, function(x) sum(x^2)))/(T * sum(Nm))
      sigma2_y = sum(sapply(y, function(x) sum(x^2)))/(T * sum(Nm))
      if(weight == TRUE){
        Nmin = mean(Nm)
        Nmin = min(Nm)
      }else{
        Nmin = min(Nm)
      }
      threshold = 1 - (sqrt(Nmin) + sqrt(T))/(sqrt(Nmin * T)) * log(log(sqrt(min(Nmin, T)))) * exp(sigma2_e/sigma2_y) * (1 - 1/M)
      r0hat = sum(rho >= threshold)
    }
  } else {
    threshold = NULL
    r0hat = r0
  }

  # estimate G
  if (r0hat > 0) {
    Ghat = eig_deco$vectors[, 1:r0hat]
    Ghat = sqrt(T) * as.matrix(Ghat)
    Proj_G = Ghat %*% solve(t(Ghat) %*% Ghat) %*% t(Ghat)
    loading_G = list()
    for(m in 1:M){
      loading_G[[m]] = 1/T*t(y[[m]]) %*% Ghat
    }
  } else {
    Ghat = NA
    Proj_G = matrix(0, T, T)
    loading_G = NA
  }

  # estimate F
  if(localfactor == FALSE){
    res = list(r0hat = r0hat, rho = rho, Ghat = Ghat, loading_G = loading_G, threshold = threshold)
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
      if (!(all(r%%1 == 0) && all(r >= 0))) {
        stop("invalid 'r' input")
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
               loading_G = loading_G, loading_F = loading_F, residual = e, threshold = threshold)
  }
  class(res) = "GFA"
  return(res)
}
