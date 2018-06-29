#' plot the epspath, solution paths of SVR as a function of epsilon
#'
#' produces a plot of the SVR \code{epsilon} path.
#'
#' @param x The epspath object
#' @param intercept If it is \code{TRUE}, then an intercept path plot is given.
#' @param ... Generic compatibility
#' @return The entire solution path of SVR solution as a function of \code{epsilon}.
#' @author Do Hyun Kim, Seung Jun Shin
#' @examples \donttest{
#' # The 'eobj' is given by examples description of epspath().
#' plot(eobj, lty = 2, lwd = 2, col = 2, cex.lab = 1.5) }
#' @importFrom graphics abline lines par plot
#'
#' @export plot.epspath
#' @export
plot.epspath <- function(x, intercept = FALSE, ...){

  svr.eps <- x$svr.eps
  theta <- x$theta
  theta0 <- x$theta0
  param.kernel <- x$param.kernel
  kernel.function <- x$kernel.function
  len <- length(svr.eps)
  new <- TRUE
  if(intercept == T){
    par(mfrow = c(1,1), mar = c(5,5,4,2) + 0.1)
    plot(x$svr.eps, theta0, type = "l", ylab = "intercept", xlab = expression(epsilon),
         main = "Epsilon Path for SVR (intercept)",...)
    abline(v = x$svr.eps, col = "gray", lty=2)
    new <- FALSE
  }
  if(new == T){
  par(mfrow = c(1,1), mar = c(5,5,4,2) + 0.1)
  plot(svr.eps, theta[1,], type="n", xlim = c(0, max(svr.eps))
                , ylim =c(-1,1), ylab=expression(theta), xlab = expression(epsilon),...)
  for(i in 2:length(y)){
    lines(svr.eps, theta[i,],...)
  }
  title("Epsilon Path for SVR")
  abline(v=svr.eps, col = "gray", lty=2)
  }
  y <- NULL
}
