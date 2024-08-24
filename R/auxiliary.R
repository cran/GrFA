
print.GFA <- function(x, ...){
  if (!inherits(x, "GFA"))
    stop("Not a legitimate \"GFA\" object")
  cat("The number of global factors is:", x$r0hat, "\n")
  if(!is.null(x$rhat)) cat("The number of local factors are:", x$rhat, "\n")
  if(!is.null(x$threshold)) cat("The threshold is:", x$threshold, "\n")
  cat("The reference statistics are:", round(x$rho, 4), "\n")
}

corr <- function(mat){
  mat = t(mat)
  mat = mat - rowMeans(mat)
  # Standardize each variable
  mat = mat / sqrt(rowSums(mat^2));
  # Calculate correlations
  cr = tcrossprod(mat)
  return(cr)
}

cov_my <- function(a, b){
  t(a) %*% b
}

cca_my <- function(a, b, r){
  res = eigen(solve(cov_my(a, a)) %*% cov_my(a, b) %*% solve(cov_my(b, b)) %*% cov_my(b, a))
  cor = sqrt(res$values[1:r])
  xcoef = res$vectors[, 1:r]
  xscore = a %*% xcoef
  list(cor = cor, xcoef = xcoef, xscore = xscore)
}

tr <- function(x){
  sum(diag(x))
}

TraceRatio <- function(G, Ghat){
  if(ncol(G) > 0 & !all(is.na(Ghat))){
    # TR:
    TR = tr(t(G) %*% Ghat %*% solve(t(Ghat) %*% Ghat) %*% t(Ghat) %*% G) / tr(t(G) %*% G)
    return(TR)
  } else{
    return(0)
  }
}

