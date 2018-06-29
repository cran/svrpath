#' plot the svrpath, solution paths of SVR as a function of lambda
#'
#' produces a plot of the SVR \code{lambda} path.
#'
#' @param x The svrpath object
#' @param intercept If it is \code{TRUE}, then an intercept path plot is given.
#' @param ... Generic compatibility
#' @return The entire solution path of SVR solution as a function of \code{lambda}.
#' @author Do Hyun Kim, Seung Jun Shin
#' @examples
#' \donttest{
#' # The 'obj' is given by examples description of svrpath().
#' plot(obj, lty = 2, lwd = 2, col = 2, cex.lab = 1.5) }
#' @importFrom graphics abline lines plot title
#'
#' @export
#' @export plot.svrpath
#'
plot.svrpath <- function(x, intercept = FALSE, ...){

  lambda <- x$lambda
  result <- x$theta
  theta0 <- x$theta0
  len <- length(lambda)

  result <- cbind(result, rep(0, nrow(result)))
  lambda <- c(lambda,0)
  new <- TRUE
  if(intercept == T){
    par(mfrow = c(1,1), mar = c(5,5,4,2) + 0.1)
    plot(x$lambda, theta0, type = "l", ylab = "intercept", xlab = expression(lambda),
         main = "Regularization Path for SVR (intercept)",...)
    abline(v = x$lambda, col = "gray", lty=2)
    new <- FALSE
  }
  if(new == T){
  par(mfrow = c(1,1), mar = c(5,5,4,2) + 0.1)
  plot(lambda, result[1,], type="n", xlim = c(0, max(lambda))
                , ylim =c(-1,1), ylab=expression(theta),xlab = expression(lambda), ...)
  for(i in 2:(length(y))){
    lines(lambda, result[i,],...)
  }
  title("Regularization Path for SVR")
  abline(v = lambda, col = "gray", lty=2)
  }
  y <- NULL
}
