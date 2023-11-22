gendata <- function(seed = 1, T = 50, N = rep(30, 5),
                    r0 = 2, r = rep(2, 5),
                    Phi_G = 0.5, Phi_F = 0.5, Phi_e = 0.5,
                    W_F = 0.5, beta = 0.1, kappa = 1, case = 1){
  # check the validation of input parameters
  if(r0 %% 1 != 0 | r0 < 0) stop("invalid r0 input")
  if(!(all(r %% 1 == 0) && all(r > 0))) stop("invalid r input")
  if(length(r) < 2) stop("invalid r0 input")
  if(Phi_G < 0 | Phi_G >= 1) stop("invalid Phi_G input")
  if(Phi_F < 0 | Phi_F >= 1) stop("invalid Phi_F input")
  if(Phi_e < 0 | Phi_e >= 1) stop("invalid Phi_e input")
  if(W_F < 0 | W_F >= 1) stop("invalid W_F input")
  if(beta < 0 | beta >= 1) stop("invalid beta input")
  if(kappa < 0) stop("invalid kappa input")
  if(!all(case %in% c(1, 2, 3))) stop("invalid case input")
  if(case == 2){
    if(length(unique(r)) > 1) stop("invalid r input")
  }
  set.seed(seed)
  M = length(r)
  F = list()
  y = list()
  loading_F = list()
  loading_G = list()
  e = list()
  # global factors
  G = matrix(0, T, r0)
  G[1, ] = rnorm(r0)
  for(t in 2:T){
    G[t, ] = Phi_G*G[t-1, ] + rnorm(r0)
  }
  # local factors
  if(case == 1){
    for(m in 1:M){
      Fm = matrix(0, T, r[m])
      Fm[1, ] = rnorm(r[m])
      for(t in 2:T){
        Fm[t, ] = Phi_F*Fm[t-1, ] + rnorm(r[m])
      }
      F[[m]] = Fm
    }
  } else if(case == 2){
    F1 = matrix(0, T, r[1])
    F1[1, ] = rnorm(r[1])
    for(t in 2:T){
      F1[t, ] = Phi_F*F1[t-1, ] + rnorm(r[1])
    }
    for(m in 1:floor(M/2)){
      F[[m]] = F1
    }
    F2 = matrix(0, T, r[1])
    F2[1, ] = rnorm(r[1])
    for(t in 2:T){
      F2[t, ] = Phi_F*F2[t-1, ] + rnorm(r[1])
    }
    for(m in (1+floor(M/2)):M){
      F[[m]] = F2
    }
  } else{
    r_all = sum(r)
    Omega = matrix(W_F, r_all, r_all)
    diag(Omega) = 1
    F_all = rmvnorm(n = T, sigma = Omega)
    dim_1 = c(1, 1 + cumsum(r))[1:M]
    dim_2 = cumsum(r)
    for(m in 1:M){
      F[[m]] = F_all[, dim_1[m]:dim_2[m]]
    }
  }
  # global and local factor loadings
  for(m in 1:M){
    loading_F[[m]] = matrix(rnorm(N[m]*r[m]), N[m], r[m])
    loading_G[[m]] = matrix(rnorm(N[m]*r0), N[m], r0)
  }
  # scale factor
  theta1 = rep(0, M)
  theta2 = rep(0, M)
  if(r0 == 0){
    for(m in 1:M){
      theta1[m] = 1
      theta2[m] = r[m]/(1-Phi_G^2)/((1+16*beta^2)/(1-Phi_e^2))
    }
  } else {
    for(m in 1:M){
      theta1[m] = r0/(1-Phi_G^2)/(r[m]/(1-Phi_F^2))
      theta2[m] = r0/(1-Phi_G^2)/((1+16*beta^2)/(1-Phi_e^2))
    }
  }
  # idiosyncratic error
  for(m in 1:M){
    em = matrix(0, T, N[m])
    em[1, ] = rnorm(N[m])
    v = matrix(rnorm(T*(N[m]+2*8)), T, N[m]+2*8)
    for(t in 2:T){
      for(i in 1:N[m]){
        em[t, i] = Phi_e*em[t-1, i] + v[t, i + 8] + beta*sum(v[t, i:(i+8-1)]) + beta*sum(v[t, (i+8+1):(i+2*8)])
      }
    }
    e[[m]] = em
  }
  for(m in 1:M){
    y[[m]] = G %*% t(loading_G[[m]]) + sqrt(theta1[m])*F[[m]] %*% t(loading_F[[m]]) + sqrt(kappa*theta2[m])*e[[m]]
  }
  res = list(y = y, G = G, F = F, loading_G = loading_G, loading_F = loading_F, T = T, N = N, M = M, r0 = r0, r = r, case = case)
  class(res) = "GFD"
  return(res)
}


print.GFD <- function(x, ...){
  if (!inherits(x, "GFD"))
    stop("Not a legitimate \"GFD\" object")
  cat("Information of the data:", "\n")
  cat("T:", x[["T"]], "\n")
  cat("N:", x[["N"]], "\n")
  cat("groups numbers M:", x[["M"]], "\n")
  cat("global factors r0:", x[["r0"]], "\n")
  cat("local factors r:", x[["r"]], "\n")
  cat("case:", x[["case"]], "\n")
}
