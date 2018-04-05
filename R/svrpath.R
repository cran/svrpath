#' Fit the entire regularization path for Support Vector Regression
#'
#' This algorithm computes the entire regularization path for the support vector regression  with a relatively low cost compared to quadratic programming problem.
#'
#' @param x The data matrix (n x p) with n rows (observations) on p variables (columns)
#' @param y The real number valued response variable
#' @param svr.eps epsilon in epsion-insensitive loss function
#' @param kernel.function This is a user-defined function. Provided are \code{poly.kernel} (the default, with parameter set to default to a linear kernel) and \code{radial.kernel}
#' @param param.kernel the parameter(s) for the kernel. For this radial kernel, the parameter is known in the fields as "gamma". For the polynomial kernel, it is the "degree"
#' @param lambda.min The smallest value of lambda for termination of the algorithm. Default is \code{lambda = 1e-8}
#' @param ridge Sometimes the algorithm encounters singularities; in this case a small value of ridge can help, default is \code{ridge = 1e-8}#'
#' @param eps a small machine number which is used to identify minimal step sizes
#' @param ... Generic compatibility
#' @return a 'svrpath' object is returned, for which there are print, \code{lambda} values, corresponding values of \code{theta} for each data point.
#' @author Dohyun Kim, Seung Jun Shin
#' @seealso \code{\link{predict.svrpath}}, \code{\link{plot.svrpath}}
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
#' @importFrom svmpath poly.kernel radial.kernel UpdateKstar SolveKstar
#' @export
svrpath <- function(x, y, svr.eps = 1,
                    kernel.function = radial.kernel,
                    param.kernel = 1,
                    ridge = 1e-8,
                    eps = 1e-8,
                    lambda.min = 1e-8,...){
  n <- length(y)
  Nmoves <- 5*n
  options(digits=6)
  Kscript <- kernel.function(x, x, param.kernel = param.kernel)
  theta <- matrix(1, n, Nmoves)
  theta0 <- double(Nmoves)
  Elbow.L.list <- as.list(seq(Nmoves))
  Elbow.R.list <- as.list(seq(Nmoves))
  Elbow.Length <- as.list(seq(Nmoves))
  lambda <- double(Nmoves)

  # Initialization
  init  <- Initialization(x, y, svr.eps = svr.eps, kernel.function = kernel.function, param.kernel = param.kernel)
  Elbow.R <- init$Elbow.R
  Elbow.L <- init$Elbow.L
  Right <- init$Right
  Left  <- init$Left
  Center <- init$Center
  lambda0 <- init$lambda
  theta0 <- init$theta0 * lambda0
  theta[,1] <- init$theta

  Elbow.L.list[[1]] <- Elbow.L
  Elbow.R.list[[1]] <- Elbow.R
  Elbow.Length[[1]] <- length(Elbow.L) + length(Elbow.R)

  if(length(init$Elbow.L) == 0 & length(init$Elbow.R) ==0) {
    stop("all points are already in the center")
  }

  ystar <- c(0, y[Elbow.R] - svr.eps, y[Elbow.L] + svr.eps)
  updated.index <- c(Elbow.R.list[[1]], Elbow.L.list[[1]])
  Kstar <- UpdateKstar(0, Kscript[updated.index, updated.index],
                       NULL, rep(1, length(updated.index)))
  lambda[1] <- lambda0

  fl <- (Kscript %*% theta[,1] + theta0)/lambda0
  k <- 1

  while ((k < Nmoves) && (lambda[k] > lambda.min)) {
    if (length(Elbow.L) == 0 & length(Elbow.R) == 0) {
      re.init <- re.Initialization(x,y,Left,Right,Center,lambda=lambda[k], svr.eps = svr.eps,
                                    kernel.function = kernel.function, param.kernel = param.kernel)
      Elbow.L <- re.init$Elbow.L
      Elbow.R <- re.init$Elbow.R
      Left <- re.init$Left
      Right <- re.init$Right
      Center <- re.init$Center
      theta[,k+1] <- re.init$theta
      lambda[k+1] <- re.init$lambda0
      theta0[k+1] <- re.init$theta0 * lambda[k+1]
      fl <- (Kscript %*% (theta[,k+1]) + theta0[k+1])/lambda[k+1]

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
      dl <- fl - gl

      immobile <- sum(abs(dl))/n < eps
      temp.L <- -theta[Elbow.L, k] + lambda[k]*b.L

      lambda.Elbow.L.to.L <- (temp.L-1)/b.L
      lambda.Elbow.L.to.L[abs(b.L) < eps] <- -1

      lambda.Elbow.L.to.C <- temp.L/b.L
      lambda.Elbow.L.to.C[abs(b.L) < eps] <- -1

      temp.R <- -theta[Elbow.R, k] + lambda[k]*b.R

      lambda.Elbow.R.to.R <- (temp.R+1)/b.R
      lambda.Elbow.R.to.R[abs(b.R) < eps] <- -1

      lambda.Elbow.R.to.C <- temp.R/b.R
      lambda.Elbow.R.to.C[abs(b.R) < eps]  <- -1

      lambda01 <- c(lambda.Elbow.L.to.L, lambda.Elbow.L.to.C, lambda.Elbow.R.to.R, lambda.Elbow.R.to.C)
      lambda.exit <- lambda01[lambda01 < lambda[k] - eps]

      if(length(lambda.exit) == 0){
        lambda.exit <- -1
      }else{
        lambda.exit <- max(lambda01[lambda01 < lambda[k] - eps])
      }
      lambda.exit

      if(immobile & (abs(lambda.exit) < 1e-15)) break
      if(!immobile){
        lambda.hit.Elbow.L <- (lambda[k] * (dl))/(y - gl + svr.eps)
        lambda.hit.Elbow.R <- (lambda[k] * (dl))/(y - gl - svr.eps)
        lambda.hit.Elbow.L[abs(y-gl+svr.eps) < eps] <- -Inf
        lambda.hit.Elbow.R[abs(y-gl-svr.eps) < eps] <- -Inf

        lambda.entry.L <- max(lambda.hit.Elbow.L[lambda.hit.Elbow.L < lambda[k]-eps])
        lambda.entry.R <- max(lambda.hit.Elbow.R[lambda.hit.Elbow.R < lambda[k]-eps])

        c(lambda.entry.L, lambda.entry.R)
        lambda.entry <- max(lambda.entry.L, lambda.entry.R)

      }else lambda.entry <- -1

      lambda.max <- max(lambda.entry, lambda.exit)
      if(lambda.max < - 1.0e-5 ) break
      lambda[k+1]   <- lambda.max
      theta0[k + 1] <- theta0[k] + ((lambda[k+1] - lambda[k]) * b0)
      theta[,k+1]   <- theta[,k]
      theta[Elbow.R, k + 1] <- theta[Elbow.R, k] + ((lambda[k+1] - lambda[k])*b.R)
      theta[Elbow.L, k + 1] <- theta[Elbow.L, k] + ((lambda[k+1] - lambda[k])*b.L)
      fl <- (lambda[k]/lambda[k + 1]) * (dl) + gl

      if(lambda.entry > lambda.exit){
        if(lambda.entry.R >lambda.entry.L){
          i <- match(lambda.entry, lambda.hit.Elbow.R, 0)[1]
          obs <- i
          if(match(i, Right, FALSE)){
            Right <- setdiff(Right, i)
          }else{
            Center <- setdiff(Center, i)
          }
          Elbow.R <- c(Elbow.R, i)
          i <- NULL
        }else{
          i <- match(lambda.entry, lambda.hit.Elbow.L, 0)[1]
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
        i <- Elbow.L[abs(lambda.Elbow.L.to.L - lambda.exit) < eps]
        if (length(i) > 0) {
          Leaveleft <- rep(TRUE, length(i))
          idrop <- i
        }
        i <- Elbow.L[abs(lambda.Elbow.L.to.C - lambda.exit) < eps]
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
        i <- Elbow.R[abs(lambda.Elbow.R.to.R - lambda.exit) < eps]
        if (length(i) > 0) {
          Leaveright <- rep(TRUE, length(i))
          idrop <- i
        }
        i <- Elbow.R[abs(lambda.Elbow.R.to.C - lambda.exit) < eps]
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
    ystar <- c(0, y[Elbow.R] - svr.eps, y[Elbow.L] + svr.eps)
    updated.index <- c(Elbow.R, Elbow.L)
    if(length(updated.index) != 0){
      Kstar <- UpdateKstar(0, Kscript[updated.index, updated.index],
                           NULL, rep(1, length(updated.index)))
      Kstar <- Kstar + diag(length(ystar)) * ridge
    }
  }
  obj <- list(theta = round(theta[ ,seq(k), drop = FALSE],5),
              theta0 = round(theta0,5),
              lambda = round(lambda[seq(k)],5),
              Elbow.L = Elbow.L.list[seq(k)], Elbow.R = Elbow.R.list[seq(k)],
              Elbow.Length = Elbow.Length[seq(k)],
              Left = Left, Center = Center, Right = Right,
              eps = svr.eps, Kscript = Kscript, x = x, y = y,
              kernel.function = kernel.function, param.kernel = param.kernel)
  class(obj) <- 'svrpath'
  obj
}
