#' @rdname svrpath-internal
#' @export
e.Initialization <- function(x, y, lambda = 1,
                             kernel.function = radial.kernel, param.kernel = 1){

  Center <- c(1:length(y))
  Right <- integer(0)
  Left  <- integer(0)
  Elbow.R <- Elbow.L <- NULL
  K <- kernel.function(x, x, param.kernel = param.kernel)

  lambda1 <- (outer(y[Center], y[Center], "-") / 2)
  max.lambda1 <- lambda1[which.max(lambda1)]
  lambda0 <- max.lambda1
  temp <- lambda1

  i1 <- row(matrix(0, nrow(temp), ncol(temp)))[which.max(temp)]
  i2 <- col(matrix(0, nrow(temp), ncol(temp)))[which.max(temp)]

  Elbow.R <- Center[i1]
  Elbow.L <- Center[i2]
  Center <- setdiff(Center, c(Elbow.R, Elbow.L))

  if(length(Elbow.R)==0){
    Elbow.R <- integer(0)
  }else if(length(Elbow.L) == 0){
    Elbow.L <- integer(0)
  }

  eps <- lambda0
  theta0 <- (y[Elbow.R] + y[Elbow.L]) / 2
  theta <- rep(0, length(y))

  list(Elbow.L = Elbow.L, Elbow.R=Elbow.R, Center=Center,
       Right=Right, Left=Left, eps0 = eps,
       theta0 = theta0, theta=theta)
}
