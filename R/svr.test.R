svr.test <- function(newy, newx, obj, kernel.function = radial.kernel, param.kernel = 1){
  avg.loss <- rep(0, length(obj$lambda))
  n <- length(newy)
  x.train <- obj$x
  K <- kernel.function(x.train, newx, param.kernel = param.kernel)
  svr.eps <- obj$eps

  for(i in 1:length(obj$lambda)){
    lambda <- obj$lambda[i]
    fx <- ( 1/lambda ) * (obj$theta[,i]%*%obj$Kscript + obj$theta0[i])
    loss_test_avg <- loss(t(fx), newy, eps = svr.eps)
    avg.loss[i] <- loss_test_avg / n
  }

  lambda2 <- obj$lambda[which.min(avg.loss)]
  theta2 <- obj$theta[,which.min(avg.loss)]
  theta2.0 <- obj$theta0[which.min(avg.loss)]
  elbow.l <- obj$Elbow.L[[which.min(avg.loss)]]
  elbow.r <- obj$Elbow.R[[which.min(avg.loss)]]

  return(list(Average.Loss = avg.loss, optimal.lambda = lambda2,
              Elbow.L = elbow.l, Elbow.R = elbow.r,
              theta = theta2, theta0 = theta2.0))
}
