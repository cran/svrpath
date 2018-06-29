svr.sic <- function(obj){

  sic <- rep(0, length(obj$lambda))
  n <- length(obj$y)
  y <- obj$y
  svr.eps <- obj$eps

  for(i in 1:length(obj$lambda)){
    fx <- ( 1/obj$lambda[i]) * (obj$theta[,i]%*%obj$Kscript + obj$theta0[i])
    loss_sic <- loss(t(fx), y, svr.eps)
    sic[i] <- log( loss_sic / n ) + ( log(n) / (2*n) )*( length(obj$Elbow.L[[i]]) + length(obj$Elbow.R[[i]]) )
  }

  lambda1 <- obj$lambda[which.min(sic)]
  theta1  <- obj$theta[,which.min(sic)]
  theta1.0    <- obj$theta0[which.min(sic)]
  elbow.l <- obj$Elbow.L[[which.min(sic)]]
  elbow.r <- obj$Elbow.R[[which.min(sic)]]

  return(list(SIC = sic, optimal.lambda = lambda1,
              Elbow.L = elbow.l, Elbow.R = elbow.r,
              theta = theta1, theta0 = theta1.0))
}
