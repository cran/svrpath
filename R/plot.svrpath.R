#' plot the svrpath, solution paths of SVR as a function of lambda
#'
#' produces a plot of the SVR \code{lambda} path.
#'
#' @param x The svrpath object
#' @param intercept if it is \code{TRUE}, then provides a intercept path plot.
#' @param ... Generic compatibility
#' @return The entire solution path of SVR solution as a function of \code{lambda}.
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
#' svr.eps <- 1
#' obj <- svrpath(x, y, svr.eps = svr.eps)
#' plot(obj)
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
    plot(x$lambda, theta0, type = "l", ylab = "intercept", xlab = "lambda",
         main = "Regularization Path for SVR (intercept)",...)
    abline(v = x$lambda, col = "gray", lty=2)
    new <- FALSE
  }
  if(new == T){
  plot(lambda, result[1,], type="n", xlim = c(0, max(lambda))
                , ylim =c(-1,1), ylab="theta",...)
  for(i in 2:(length(y))){
    lines(lambda, result[i,],...)
  }
  title("Regularization Path for SVR")
  abline(v = lambda, col = "gray", lty=2)
  }
  y <- NULL
}
