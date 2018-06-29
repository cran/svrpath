#' QP solver for SVR
#'
#' solves quadratic programming(QP) for SVR.
#'
#' @param a The data matrix (n x p) with n rows (observations) on p variables (columns)
#' @param b The real number valued response variable
#' @param lambda The regularization parameter
#' @param svr.eps Epsilon in epsion-insensitive loss function
#' @param kernel.function User defined kernel function. See \code{svmpath}.
#' @param param.kernel Parameter(s) of the kernels. See \code{svmpath}.
#' @param ... Generic compatibility
#'
#' @return SVR solution at a given \code{lambda} and \code{epsilon}
#' @author Dohyun Kim, Seung Jun Shin
#' @examples
#' \donttest{
#' # set.seed(1)
#' n <- 30
#' p <- 50
#'
#' x <- matrix(rnorm(n*p), n, p)
#' e <- rnorm(n, 0, 1)
#' beta <- c(1, 1, rep(0, p-2))
#' y <- x %*% beta + e
#' solve.svr(x, y) }
#'
#' @importFrom quadprog solve.QP
#' @importFrom svmpath poly.kernel radial.kernel UpdateKstar SolveKstar
#'
#' @export solve.svr
#' @export
solve.svr <- function(a, b, lambda = 1, svr.eps = 1,
                      kernel.function = radial.kernel, param.kernel = 1,...){
  a <- as.matrix(a)
  K <- kernel.function(a, a, param.kernel = param.kernel)
  eps <- svr.eps
  digits <- 6

  K <- rbind(cbind(K,-K), cbind(-K, K))
  Dmat <- K
  diag(Dmat) <- diag(Dmat) + 1e-5
  dvec <- c(-eps+b, -eps-b)
  Amat <- cbind(c(rep(1,length(b)), rep(-1, length(b))),
                diag(2*length(b)), -diag(2*length(b)))
  bvec <- c(0,  rep(0, 2*length(b)), rep(-1, 2*length(b)))
  Dmati  <- Dmat/lambda
  result <- quadprog::solve.QP(Dmat=Dmati,dvec=dvec, Amat=Amat,bvec=bvec, meq=1)
  alpha <- result$solution
  alpha1 <- alpha[seq(length(alpha)/2)]
  alpha2 <- alpha[seq(from = length(alpha)/2 + 1, length(alpha))]

  w <- alpha1 - alpha2
  w <- w%*%a/lambda

  # alpha0
  e.right.index <- which(0 + 1e-8  < alpha1 & alpha1 < 1 - 1e-8)
  e.right.index <- e.right.index[round(alpha2[e.right.index],3) == 0 ]

  e.left.index  <- which(0 + 1e-8 < alpha2 & alpha2 < 1 - 1e-8)
  e.left.index  <- e.left.index[round(alpha1[e.left.index],3) == 0 ]

  alpha0.right <- rep(0, length(e.right.index))
  alpha0.left  <- rep(0, length(e.left.index))

  if(length(e.right.index) ==0){
    alpha0.right <- NA
  }else if(length(e.right.index) == 1 | ncol(x) == 1){
    alpha0.right <- b[e.right.index] - w%*%a[e.right.index,] - eps
  }else{
    for(i in 1:length(e.right.index)){
      alpha0.right[i] <- b[e.right.index][i] - w%*%a[e.right.index,][i,] - eps
    }
  }

  if(length(e.left.index) == 0){
    alpha0.left <- NA
  }else if(length(e.left.index) == 1 | ncol(a) == 1){
    alpha0.left  <- b[e.left.index] - w%*%a[e.left.index,] + eps
  }else{
    for(i in 1:length(e.left.index)){
      alpha0.left[i] <- b[e.left.index][i] - w%*%a[e.left.index,][i,] + eps
    }
  }
  k <- x <- y <- NULL
  theta0 <- lambda*mean(c(alpha0.left, alpha0.right), na.rm=T) # lambda * b
  theta  <- alpha1 - alpha2
  result <- list(theta = round(theta, 3), theta0 = theta0, lambda = lambda, svr.eps=svr.eps)
  class(result) <- "solve.svr"
  result
}
