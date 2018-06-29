#' Make predictions from an "epspath" object
#'
#' Provides a prediction value at a given \code{epsilon} from \code{epspath} object.
#' @param object The epspath object
#' @param newx Values of x to be predicted. This is a matrix with observations per row. Default is x in the epspath object.
#' @param svr.eps The value of the "epsilon-insensitive loss" paramter, epsilon.
#' @param ... Generic compatibility
#' @return In each case, the desired prediction.
#' @author Do Hyun Kim, Seung Jun Shin
#' @examples
#' \donttest{
#' # The 'eobj' is given by examples description of epspath().
#' predict(eobj, svr.eps = .1) }
#' @export predict.epspath
#' @export
predict.epspath <- function(object, newx, svr.eps = 1,...){
  kernel.function <- object$kernel.function
  param.kernel <- object$param.kernel

  if(missing(newx)){
    newx <- object$x
  }

  K <- kernel.function(newx, object$x, param.kernel)
  lambda <- object$lambda
  otheta <- object$theta
  otheta0 <- object$theta0
  osvr.eps <- object$svr.eps
  minl <- min(osvr.eps)
  maxl <- max(osvr.eps)

  if(length(object$svr.eps) == 1){
    warning('This model has only one epsilon value; prediction values at this epsilon are returned')
    nlam <- length(svr.eps)
    svr.eps <- rep(object$svr.eps, nlam)
    theta <- otheta[rep(1, nlam), ,drop=FALSE]
    theta0 <- rep(otheta0, nlam)
  }else if(length(svr.eps) == 1){

    if(svr.eps >= maxl){
      theta <- object$theta[,1]
      theta0 <- object$theta0[1]
    }else if(svr.eps <= minl){
      theta <- object$theta[,ncol(object$theta)]
      theta0 <- object$theta0[ncol(object$theta)]
      warning("The epsilon value is too small to predict values.")
    }else{
      theta <- rep(1, times = nrow(otheta))
      for(i in 1:nrow(otheta)){
        theta[i] <-  approx(osvr.eps, otheta[i,], svr.eps)$y
      }
      theta0 <- approx(osvr.eps, otheta0, svr.eps)$y
    }

    fx <- (K%*%(theta) + theta0) / lambda

  }else if(length(svr.eps) > 1){

    theta <- otheta[,1:length(svr.eps)]
    theta0 <- otheta0[1:length(svr.eps)]

    for(i in 1:length(svr.eps)){
      if(svr.eps[i] >= maxl){
        theta[,i] <- object$theta[,1]
        i <- i+1
      }else if(svr.eps[i] <= minl){
        theta[,i] <- rep(0, nrow(otheta))
        i <- i+1
      }
      if(i > length(svr.eps))break

      for(j in 1:nrow(otheta)){
        theta[j,i]  <-  approx(osvr.eps, otheta[j,], svr.eps[i])$y
      }
      theta0[i] <- approx(osvr.eps, otheta0, svr.eps[i])$y
    }

    if(nrow(newx) > 1){
      theta0 <- matrix(rep(theta0,nrow(newx)), ncol = length(svr.eps), byrow =T)
      lambda <- matrix(rep(lambda,2*nrow(newx)), ncol = length(svr.eps), byrow = T)

      fx <- (K%*%(theta) + theta0) / lambda
      theta0 <- theta0[1,]
      lambda <- lambda[1,1]
    }else{
      fx <- (K%*%(theta) + theta0) / lambda
    }
  }

  obj <- list(fx = fx, theta = theta, theta0 = theta0, svr.eps = svr.eps, lambda = lambda)
  class(obj) <- "predict.epspath"
  obj
}
