re.e.Initialization <- function(x,y,Left, Right, Center, svr.eps = svr.eps[k], lambda = 1,
                                kernel.function = radial.kernel, param.kernel = 1){

  Left.cp <- Left
  Right.cp <- Right
  lambdak <- lambda
  Elbow.L <- Elbow.R <- NULL
  eps <- svr.eps
  mod <- 1e-5
  K <- kernel.function(x, x, param.kernel = param.kernel)

  gx <- apply(K[,Right, drop = F], 1, sum) - apply(K[,Left, drop = F], 1, sum)
  lambda1 <- (outer(y[Center], y[Center], "-") / 2)
  temp <- lambda1
  l1 <- lambda1[lambda1 < lambdak - mod]

  if(length(l1) == 0){ l1 <- 0 }
  lambda0 <- max(l1)

  if(ncol(as.matrix(temp)) == 1 && nrow(as.matrix(temp)) == 1){
    i1 <- 1
    i2 <- 1
  }else{
    i1 <- which(temp == lambda0, arr.ind = T)[,1]
    i2 <- which(temp == lambda0, arr.ind = T)[,2]
  }

  Elbow.R <- Center[i1]
  Elbow.L <- Center[i2]
  Center <- setdiff(Center, c(Elbow.R, Elbow.L))

  eps0 <- lambda0

  theta0 <- (y[Elbow.R] + y[Elbow.L]) / 2
  theta <- rep(0, length(y))
  theta[c(Right)] <- 1
  theta[c(Left)] <- -1

  k <- NULL
  list(Elbow.L = Elbow.L, Elbow.R=Elbow.R, Center=Center,
       Right=Right, Left=Left, theta0=theta0,
       theta=theta, eps0 = eps0)
}
