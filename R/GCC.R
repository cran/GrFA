GCC <- function(y, rmax = 8, r0 = NULL, r = NULL, type = "BIC3"){
  M = length(y)
  T = nrow(y[[1]])
  Nm = sapply(y, ncol)
  CNT = min(c(T, Nm))
  K = list()
  for(m in 1:M){
    K[[m]] = FA(y[[m]], r = rmax)$F
  }
  Phi = matrix(0, T*M*(M-1)/2, M*rmax)
  row = 1
  for(i in 1:(M-1)){
    for(j in (i+1):M){
      Phi[row:(row + T-1), ((i-1)*rmax+1):(i*rmax)] = K[[i]]
      Phi[row:(row + T-1), ((j-1)*rmax+1):(j*rmax)] = -K[[j]]
      row = row + T
    }
  }
  svd_phi = svd(Phi)
  delta2 = sort(svd_phi$d^2)
  delta2 = c(sum(delta2)/(CNT*M*rmax), delta2)
  rho = delta2[2:(rmax + 1)]/delta2[1:rmax]
  if(is.null(r0)){
    r0hat = which.max(rho) - 1
  }else{
    if (!(r0%%1 == 0) | r0 < 0) {
      stop("invalid 'r0' input")
    }
    r0hat = r0
  }
  if(r0hat > 0){
    psi = matrix(0, T, 0)
    Q = svd_phi$v[, 1:rmax]
    for(m in 1:M){
      psi = cbind(psi, K[[m]] %*% Q[((m-1)*rmax+1):(m*rmax),])
    }
    Ghat = eigen(psi %*% t(psi))$vectors[, 1:r0hat]
    Ghat = sqrt(T)*as.matrix(Ghat)
    Proj_G = Ghat %*% solve(t(Ghat) %*% Ghat) %*% t(Ghat)
  }else{
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
    if (!(all(r%%1 == 0) && all(r >= 0))){
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

  res = list(r0hat = r0hat, rhat = rhat, rho = rho, Ghat = Ghat, Fhat = Fhat,
             loading_G = loading_G, loading_F = loading_F, residual = e)
  class(res) = "GFA"
  return(res)
}





