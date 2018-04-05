#' plot the epspath, solution paths of SVR as a function of epsilon
#'
#' produces a plot of the SVR \code{epsilon} path.
#'
#' @param x The epspath object
#' @param intercept if it is \code{TRUE}, then provides a intercept path plot.
#' @param ... Generic compatibility
#' @return The entire solution path of SVR solution as a function of \code{epsilon}.
#' @author Dohyun Kim, Seung Jun Shin
#' @examples
#' set.seed(1)
#' n <- 30
#' p <- 50
#'
#' x <- matrix(rnorm(n*p), n, p)
#' e <- rnorm(n, 0, 1)
#' beta <- c(1, 1, rep(0, p-2))
#' y <- x %*% beta + e
#' lambda <- 1
#' eobj <- epspath(x, y, lambda = lambda)
#' plot(eobj)
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
    plot(x$svr.eps, theta0, type = "l", ylab = "intercept", xlab = "epsilon",
         main = "Epsilon Path for SVR (intercept)",...)
    abline(v = x$svr.eps, col = "gray", lty=2)
    new <- FALSE
  }
  if(new == T){
  plot(svr.eps, theta[1,], type="n", xlim = c(0, max(svr.eps))
                , ylim =c(-1,1), ylab="theta", xlab = "epsilon",...)
  for(i in 2:length(y)){
    lines(svr.eps, theta[i,],...)
  }
  title("Epsilon Path for SVR")
  abline(v=svr.eps, col = "gray", lty=2)
  }
  y <- NULL
}
