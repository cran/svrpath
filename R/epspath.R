#' Fit the entire \code{epsilon} path for Support Vector Regression
#'
#' The Suport Vector Regression (SVR) employs epsilon-intensive loss which ignores errors smaller than epsilon. This algorithm computes the entire paths for SVR solution as a function of \code{epsilon} at a given regularization parameter \code{lambda}, which we call \code{epsilon} path.
#' @param x The data matrix (n x p) with n rows (observations) on p variables (columns)
#' @param y The real number valued response variable
#' @param lambda The regularization parameter value.
#' @param kernel.function User defined kernel function. See \code{svmpath}.
#' @param param.kernel parameter(s) of the kernels. See \code{svmpath}.
#' @param eps.min The smallest value of epsilon for termination of the algorithm. Default is \code{eps.min = 1e-8}
#' @param ridge Sometimes the algorithm encounters singularities; in this case a small value of ridge can help, default is \code{ridge = 1e-8}
#' @param eps a small machine number which is used to identify minimal step sizes
#' @param ... Generic compatibility
#' @return a 'epspath' object is returned.
#' @author Dohyun Kim, Seung Jun Shin
#' @seealso \code{\link{predict.epspath}}, \code{\link{plot.epspath}}, \code{\link{svrpath}}
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
#' obj <- epspath(x, y, lambda = lambda)
#' @importFrom svmpath poly.kernel radial.kernel UpdateKstar SolveKstar
#' @export
epspath <- function(x, y, lambda = 1,
                     kernel.function = radial.kernel,
                     param.kernel = 1,
                     ridge = 1e-8,
                     eps = 1e-7,
                     eps.min = 1e-8,...){

  n <- length(y)
  Nmoves <- 5*n
  options(digits=6)
  Kscript <- kernel.function(x, x, param.kernel = param.kernel)
  theta <- matrix(1, n, Nmoves)
  theta0 <- double(Nmoves)
  Elbow.L.list <- as.list(seq(Nmoves))
  Elbow.R.list <- as.list(seq(Nmoves))
  Elbow.Length <- as.list(seq(Nmoves))
  svr.eps <- double(Nmoves)

  # Initialization
  init  <- e.Initialization(x, y, lambda = lambda, kernel.function = kernel.function, param.kernel = param.kernel)
  Elbow.R <- init$Elbow.R
  Elbow.L <- init$Elbow.L
  Right <- init$Right
  Left  <- init$Left
  Center <- init$Center
  svr.eps[1] <- init$eps0
  theta0 <- init$theta0 * lambda
  theta[,1] <- init$theta

  Elbow.L.list[[1]] <- Elbow.L
  Elbow.R.list[[1]] <- Elbow.R
  Elbow.Length[[1]] <- length(Elbow.L) + length(Elbow.R)

  if(length(init$Elbow.L) == 0 & length(init$Elbow.R) ==0) {
    stop("all points are already in the center")
  }

  ystar <- c(0, rep(1,length(Elbow.R)), rep(-1, length(Elbow.L)))
  updated.index <- c(Elbow.R.list[[1]], Elbow.L.list[[1]])
  Kstar <- UpdateKstar(0, Kscript[updated.index, updated.index],
                       NULL, rep(1, length(updated.index)))

  fl <- (Kscript %*% theta[,1] + theta0)/lambda
  k <- 1

  while ((k < Nmoves) && (svr.eps[k] > eps.min)) {
    if (length(Elbow.L) == 0 & length(Elbow.R) == 0) {
     if(length(Center) == 0)break
      re.init <- re.e.Initialization(x,y,Left,Right,Center,svr.eps=svr.eps[k], lambda = lambda,
                                     kernel.function = kernel.function, param.kernel = param.kernel)
      Elbow.L <- re.init$Elbow.L
      Elbow.R <- re.init$Elbow.R
      Left <- re.init$Left
      Right <- re.init$Right
      Center <- re.init$Center
      theta[,k+1] <- re.init$theta
      svr.eps[k+1] <- re.init$eps0
      theta0[k+1] <- re.init$theta0 * lambda
      fl <- (Kscript %*% theta[,k+1] + theta0[k+1])/lambda

    }else{
      bstar <- solve(Kstar)%*%ystar
      b0  <- bstar[1]
      b   <- bstar[-1]

      if(length(Elbow.L) == 0){
        b.R <- b[seq(from=1, to=length(b))]
        b.L <- integer(0)
      }else if(length(Elbow.R) ==0){
        b.R <- integer(0)
        b.L <- b[seq(from=1, to=length(b))]
      }else{
        b.R <- b[seq(from=1, to=length(Elbow.R))]
        b.L <- b[seq(from=length(Elbow.R) + 1, to=length(b))]
      }

      gl <- apply(Kscript[Elbow.R,,drop = FALSE]*b.R, 2,sum) + apply(Kscript[Elbow.L,,drop = FALSE]*b.L, 2,sum) + b0
      immobile <- sum(abs(gl))/n < eps

      temp.L <- theta[Elbow.L, k] + svr.eps[k]*lambda*b.L
      eps.Elbow.L.to.L <- (temp.L + 1)/(lambda*b.L)
      eps.Elbow.L.to.L[abs(b.L) < eps] <- -2

      eps.Elbow.L.to.C <- temp.L/(lambda*b.L)
      eps.Elbow.L.to.C[abs(b.L) < eps] <- -2

      temp.R <- theta[Elbow.R, k] + svr.eps[k]*lambda*b.R
      eps.Elbow.R.to.R <- (temp.R-1)/(lambda*b.R)
      eps.Elbow.R.to.R[abs(b.R) < eps] <- -2

      eps.Elbow.R.to.C <- temp.R/(lambda*b.R)
      eps.Elbow.R.to.C[abs(b.R) < eps]  <- -2

      eps01 <- c(eps.Elbow.L.to.L, eps.Elbow.L.to.C, eps.Elbow.R.to.R, eps.Elbow.R.to.C)
      eps.exit <- eps01[eps01 < svr.eps[k] - eps]

      if(length(eps.exit) == 0){
        eps.exit <- -2
      }else{
        eps.exit <- max(eps01[eps01 < svr.eps[k] - eps])
      }
      eps.exit

      if(immobile & (abs(eps.exit) < 1e-15)) break
      if(!immobile){
        eps.hit.Elbow.L <- (svr.eps[k] * gl + fl - y)/(gl + 1)
        eps.hit.Elbow.R <- (svr.eps[k] * gl + fl - y)/(gl - 1)
        eps.hit.Elbow.L[abs(gl + 1) < eps] <- -Inf
        eps.hit.Elbow.R[abs(gl - 1) < eps] <- -Inf

        eps.entry.L <- max(eps.hit.Elbow.L[eps.hit.Elbow.L < svr.eps[k]-eps])
        eps.entry.R <- max(eps.hit.Elbow.R[eps.hit.Elbow.R < svr.eps[k]-eps])

        c(eps.entry.L, eps.entry.R)
        eps.entry <- max(eps.entry.L, eps.entry.R)

      }else eps.entry <- -Inf

      eps.max <- max(eps.entry, eps.exit)
      if(eps.max < - 1.001) break
      svr.eps[k+1]   <- eps.max
      theta0[k + 1] <- theta0[k] + ((svr.eps[k] - svr.eps[k+1]) * b0 * lambda)
      theta[,k+1]   <- theta[,k]
      theta[Elbow.R, k + 1] <- theta[Elbow.R, k] + ((svr.eps[k] - svr.eps[k+1])*b.R*lambda)
      theta[Elbow.L, k + 1] <- theta[Elbow.L, k] + ((svr.eps[k] - svr.eps[k+1])*b.L*lambda)
      fl <- (svr.eps[k] - svr.eps[k + 1]) * gl + fl

      if(eps.entry > eps.exit){
        if(eps.entry.R > eps.entry.L){
          i <- match(eps.entry, eps.hit.Elbow.R, 0)[1]
          obs <- i
          if(match(i, Right, FALSE)){
            Right <- setdiff(Right, i)
          }else{
            Center <- setdiff(Center, i)
          }
          Elbow.R <- c(Elbow.R, i)
          i <- NULL
        }else{
          i <- match(eps.entry, eps.hit.Elbow.L, 0)[1]
          obs <- i
          if(match(i, Left, FALSE)){
            Left <- setdiff(Left, i)
          }else{
            Center <- setdiff(Center, i)
          }
          Elbow.L <- c(Elbow.L, i)
        }
      }else{
        i <- idrop <- Leaveright <- Leaveleft <- NULL
        i <- Elbow.L[abs(eps.Elbow.L.to.L - eps.exit) < eps]
        if (length(i) > 0) {
          Leaveleft <- rep(TRUE, length(i))
          idrop <- i
        }
        i <- Elbow.L[abs(eps.Elbow.L.to.C - eps.exit) < eps]
        if (length(i) > 0) {
          Leaveleft <- c(Leaveleft, rep(FALSE, length(i)))
          idrop <- c(idrop, i)
        }
        obs <- idrop
        for (j in seq(along = idrop)) {
          if (Leaveleft[j]) {
            Left <- c(Left, idrop[j])
          }
          else {
            Center <- c(Center, idrop[j])
          }
          mi <- match(idrop[j], Elbow.L)
          Elbow.L <- setdiff(Elbow.L , Elbow.L[mi])
        }

        i <- idrop <- NULL
        i <- Elbow.R[abs(eps.Elbow.R.to.R - eps.exit) < eps]
        if (length(i) > 0) {
          Leaveright <- rep(TRUE, length(i))
          idrop <- i
        }
        i <- Elbow.R[abs(eps.Elbow.R.to.C - eps.exit) < eps]
        if (length(i) > 0) {
          Leaveright <- c(Leaveright, rep(FALSE, length(i)))
          idrop <- c(idrop, i)
        }
        obs <- idrop
        for (j in seq(along = idrop)) {
          if (Leaveright[j]) {
            Right <- c(Right, idrop[j])
          }
          else {
            Center <- c(Center, idrop[j])
          }
          mi <- match(idrop[j], Elbow.R)
          Elbow.R <- setdiff(Elbow.R, Elbow.R[mi])
        }
      }}
    k <- k+1
    Elbow.L.list[[k]] <- Elbow.L
    Elbow.R.list[[k]] <- Elbow.R
    Elbow.Length[[k]] <- (length(Elbow.L) + length(Elbow.R))
    ystar <- c(0, rep(1, length(Elbow.R)), rep(-1, length(Elbow.L)))
    updated.index <- c(Elbow.R, Elbow.L)
    if(length(updated.index) != 0){
      Kstar <- UpdateKstar(0, Kscript[updated.index, updated.index],
                           NULL, rep(1, length(updated.index)))
      Kstar <- Kstar + diag(length(ystar)) * ridge
    }
  }
  obj <- list(theta = round(theta[ ,seq(k), drop = FALSE],5),
              theta0 = round(theta0,5),
              svr.eps = round(svr.eps[seq(k)],5),
              Elbow.L = Elbow.L.list[seq(k)], Elbow.R = Elbow.R.list[seq(k)],
              Elbow.Length = Elbow.Length[seq(k)],
              Left = Left, Center = Center, Right = Right,
              Kscript = Kscript, x = x, y = y, lambda = lambda,
              kernel.function = kernel.function, param.kernel = param.kernel)
  class(obj) <- "epspath"
  obj
}
