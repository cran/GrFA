FA <- function(X, r){
  T = nrow(X)
  N = ncol(X)
  if(r == 0){
    return(list(F = NA, L =  NA))
  }
  if(T > N) {
    S = t(X) %*% X/(N*T)
    eig = eigen(S)
    e = eig$vectors
    L = e[, 1:r] * sqrt(N)
    L = as.matrix(L)
    F = X %*% L/N
    F = as.matrix(F)
    F = F %*% diag(1/sqrt(diag(t(F) %*% F/T)), nrow = r, ncol = r)
    L = t(X) %*% F/T
    L = as.matrix(L)
  } else {
    S = X %*% t(X)/(N*T)
    eig = eigen(S)
    e = eig$vectors
    F = e[, 1:r] * sqrt(T)
    F = as.matrix(F)
    L = t(X) %*% F/T
    L = as.matrix(L)
  }
  return(list(F = F, L = L))
}

est_num <- function(X, kmax = 8, type = "BIC3"){
  if(is.na(match(type, c("PC1", "PC2", "PC3", "IC1", "IC2","IC3", "AIC3", "BIC3", "ER", "GR")))){
    stop("Wrong input type")
  }
  if(kmax == 0){
    return(0)
  }
  T = nrow(X)
  N = ncol(X)
  if(T > N){
    eig_val = eigen(t(X) %*% X, only.values = TRUE)$values
  } else{
    eig_val = eigen(X %*% t(X), only.values = TRUE)$values
  }
  V = rep(0, kmax + 1)
  V[1] = sum(eig_val)
  for(i in 1:kmax){
    V[i+1] = V[i] - eig_val[i]
  }
  V = V/(N*T)
  CNT = min(sqrt(N), sqrt(T))
  PC1 = V + (0:kmax)*V[kmax+1]*(N+T)/(N*T)*log(N*T/(N+T))
  PC2 = V + (0:kmax)*V[kmax+1]*(N+T)/(N*T)*log(CNT^2)
  PC3 = V + (0:kmax)*V[kmax+1]*log(CNT^2)/CNT^2
  IC1 = log(V) + (0:kmax)*(N+T)/(N*T)*log(N*T/(N+T))
  IC2 = log(V) + (0:kmax)*(N+T)/(N*T)*log(CNT^2)
  IC3 = log(V) + (0:kmax)*log(CNT^2)/CNT^2
  AIC3 = V + (0:kmax)*V[kmax+1]*2*(N+T-(0:kmax))/(N*T)
  BIC3 = V + (0:kmax)*V[kmax+1]*(N+T-(0:kmax))*log(N*T)/(N*T)
  ER = eig_val[1:kmax]/eig_val[2:(kmax + 1)]
  ER_k = which.max(ER)
  # For GR
  V = rep(0, kmax + 2)
  V[1] = sum(eig_val)
  for(i in 1:(kmax+1)){
    V[i+1] = V[i] - eig_val[i]
  }
  V = V/(N*T)
  GR = rep(0, kmax)
  for(k in 1:kmax){
    GR[k] = log(V[k]/V[k+1])/log(V[k+1]/V[k+2])
  }
  GR_k = which.max(GR)
  opt = c(which.min(PC1)-1, which.min(PC2)-1, which.min(PC3)-1,
          which.min(IC1)-1, which.min(IC2)-1, which.min(IC3)-1,
          which.min(AIC3)-1, which.min(BIC3)-1, ER_k, GR_k)
  names(opt) = c("PC1", "PC2", "PC3", "IC1", "IC2","IC3", "AIC3", "BIC3", "ER", "GR")

  # res = list(opt = opt, PC1 = PC1, PC2 = PC2, PC3 = PC3, IC1 = IC1, IC2 = IC2, IC3 = IC3,
  #            AIC3 = AIC3, BIC3 = BIC3, ER_k = ER_k, GR_k = GR_k)
  return(opt[[type]])
}


