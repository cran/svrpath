#' Make predictions from a "svrpath" object
#'
#' Provides a prediction value at a given \code{lambda} from \code{svrpath} object.
#' @param object The svrpath object
#' @param newx Values of x to be predicted. This is a matrix with observations per row. Default is x in the epspath object.
#' @param lambda The value of the regularization paramter, lambda.
#' @param criterion It provides predictions at an optimal \code{lambda} selected by SIC or GACV. \code{"sic"} or \code{"gacv"}.
#' @param ... Generic compatibility
#' @return In each case, the desired prediction.
#' @author Do Hyun Kim, Seung Jun Shin
#'
#' @examples
#' \donttest{
#' # The 'eobj' is given by examples description of epspath().
#' predict.svrpath(obj, lambda = 10) # or
#' predict(obj, criterion = 'sic') }
#'
#' @importFrom stats approx
#'
#' @export
#' @export predict.svrpath
predict.svrpath <- function(object, newx , lambda = NULL , criterion = 'sic', ...){

kernel.function <- object$kernel.function
param.kernel <- object$param.kernel

if(missing(newx)){
newx <- object$x
}

K <- kernel.function(newx, object$x, param.kernel)
otheta <- object$theta
otheta0 <- object$theta0
olambda <- object$lambda
minl <- min(olambda)
maxl <- max(olambda)

if(!is.null(lambda)){
  criterion <- NULL
}

if(is.null(lambda) && criterion == 'sic' ){
  res <- svr.sic(object)
  lambda <- res$optimal.lambda
  theta <- res$theta
  theta0 <- res$theta0
}else if(is.null(lambda) && criterion == 'gacv'){
  res <- svr.gacv(object)
  lambda <- res$optimal.lambda
  theta <- res$theta
  theta0 <- res$theta0
}

if(length(object$lambda) == 1){
warning('This model has only one lambda value; prediction values at this lambda are returned')
nlam <- length(lambda)
lambda <- rep(object$lambda, nlam)
theta <- otheta[rep(1, nlam), ,drop=FALSE]
theta0 <- rep(otheta0, nlam)
}else if(length(lambda) == 1){

if(lambda >= maxl){
  theta <- object$theta[,1]
  theta0 <- object$theta0[1]
  }else if(lambda <= minl){
    theta <- theta0 <- rep(0, nrow(otheta))
    }else{
      theta <- rep(1, times = nrow(otheta))
      for(i in 1:nrow(otheta)){
        theta[i] <-  approx(olambda, otheta[i,], lambda)$y
      }
      theta0 <- approx(olambda, otheta0, lambda)$y
    }

  fx <- (K%*%(theta) + theta0) / lambda

  }else if(length(lambda) > 1){

    theta <- otheta[,1:length(lambda)]
    theta0 <- otheta0[1:length(lambda)]

    for(i in 1:length(lambda)){
      if(lambda[i] >= maxl){
        theta[,i] <- object$theta[,1]
        i <- i+1
        }else if(lambda[i] <= minl){
          theta[,i] <- rep(0, nrow(otheta))
          i <- i+1
        }
      if(i > length(lambda))break

      for(j in 1:nrow(otheta)){
        theta[j,i]  <-  approx(olambda, otheta[j,], lambda[i])$y
        }
      theta0[i] <- approx(olambda, otheta0, lambda[i])$y
    }

    if(nrow(newx) > 1){
      theta0 <- matrix(rep(theta0,nrow(newx)), ncol = length(lambda), byrow =T)
      lambda <- matrix(rep(lambda,nrow(newx)), ncol = length(lambda), byrow = T)

      fx <- (K%*%(theta) + theta0) / lambda
      theta0 <- theta0[1,]
      lambda <- lambda[1,]
    }else{
      fx <- (K%*%(theta) + theta0) / lambda
    }

}

obj <- list(fx = fx, lambda = lambda, theta = theta, theta0 = theta0, criterion = criterion)
class(obj) <- "predict.svrpath"
obj
}

