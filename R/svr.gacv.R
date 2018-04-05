#' @rdname svrpath-internal
#' @export
svr.gacv <- function(obj){
  gacv <- rep(0, length(obj$lambda))
  n <- length(obj$y)
  y <- obj$y
  svr.eps <- obj$eps

  for(i in 1:length(obj$lambda)){
    fx <- ( 1/obj$lambda[i]) * (obj$theta[,i]%*%obj$Kscript + obj$theta0[i])
    loss_gacv <- loss(t(fx), y, svr.eps)
    gacv[i] <- loss_gacv / ( n - (length(obj$Elbow.L[[i]]) + length(obj$Elbow.R[[i]])) )
  }
  lambda3 <- obj$lambda[which.min(gacv)]
  theta3  <- obj$theta[,which.min(gacv)]
  theta3.0    <- obj$theta0[which.min(gacv)]
  elbow.l <- obj$Elbow.L[[which.min(gacv)]]
  elbow.r <- obj$Elbow.R[[which.min(gacv)]]

  return(list(GACV = gacv, optimal.lambda = lambda3,
              Elbow.L = elbow.l, Elbow.R = elbow.r,
              theta = theta3, theta0 = theta3.0))
}
