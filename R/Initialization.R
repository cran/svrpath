#' @importFrom stats optimize quantile
Initialization <- function(x, y, svr.eps= 1,kernel.function = radial.kernel, param.kernel = 1){

  eps <- svr.eps
  w0 <- optimize(loss, quantile(y, c(0,1)), y = y, eps = eps)$minimum

  Center <- (which(abs(y - w0) < eps))
  Right <- (which((y - w0) > eps))
  Left <- (which((y - w0) < -eps))

  Left.cp <- Left
  Right.cp <- Right

  Elbow.R <- Elbow.L <- NULL

  K <- kernel.function(x, x, param.kernel = param.kernel)

  gx <- apply(K[,Right, drop = F], 1, sum) - apply(K[,Left, drop = F], 1, sum)
  lambda1 <- (outer(gx[Left], gx[Right], "-")) / (outer(y[Left], y[Right], "-") + 2 * eps)
  lambda2 <- (outer(gx[Left], gx[Center], "-")) / (outer(y[Left], y[Center], "-"))
  lambda3 <- (outer(gx[Right], gx[Center], "-")) / (outer(y[Right], y[Center], "-"))
  lambda4 <- (outer(gx[Center], gx[Center], "-")) / (outer(y[Center], y[Center], "-") + 2 * eps)

  max.lambda1 <- lambda1[which.max(lambda1)]
  max.lambda2 <- lambda2[which.max(lambda2)]
  max.lambda3 <- lambda3[which.max(lambda3)]
  max.lambda4 <- lambda4[which.max(lambda4)]
  max.lambda <- c(max.lambda1,
                  max.lambda2,
                  max.lambda3,
                  max.lambda4)
  sel.case <- which.max(max.lambda)

  lambda0 <- max.lambda[sel.case]

  if (sel.case == 1) {
    temp <- lambda1
  } else if (sel.case == 2) {
    temp <- lambda2
  } else if (sel.case == 3) {
    temp <- lambda3
  } else if (sel.case == 4) {
    temp <- lambda4
  } else step()

  i1 <- row(matrix(0, nrow(temp), ncol(temp)))[which.max(temp)]
  i2 <- col(matrix(0, nrow(temp), ncol(temp)))[which.max(temp)]

  # set update
  if (sel.case == 1) {
    Elbow.L <- Left[i1]
    Elbow.R <- Right[i2]
    Left <- setdiff(Left, Elbow.L)
    Right <- setdiff(Right, Elbow.R)
  } else if (sel.case == 2) {
    Elbow.L <- Left[i1]
    Elbow.L <- c(Elbow.L, Center[i2])
    Center <- setdiff(Center, Center[i2])
    Left <- setdiff(Left, Left[i1])
  } else if (sel.case == 3) {
    Elbow.R <- Right[i1]
    Elbow.R <- c(Elbow.R, Center[i2])
    Center <- setdiff(Center, Center[i2])
    Right <- setdiff(Right, Right[i1])
  } else if (sel.case == 4) {
    Elbow.L <- Center[i1]
    Elbow.R <- Center[i2]
    Center <- setdiff(Center, c(Elbow.L, Elbow.R))
  } else step()

  if(length(Elbow.R)==0){
    Elbow.R <- integer(0)
  }else if(length(Elbow.L) == 0){
    Elbow.L <- integer(0)
  }

  theta0 <- c(y[Elbow.L] - gx[Elbow.L]/lambda0 + eps,  y[Elbow.R] - gx[Elbow.R]/lambda0 - eps)
  theta0 <- mean(theta0)
  theta <- rep(0, length(y))

  if(sel.case==1) {
    theta[c(Right,Elbow.R)] <- 1
    theta[c(Left, Elbow.L)] <- -1
  }else if(sel.case==2) {
    theta[c(Right)] <- 1
    theta[c(Left, Left.cp[i1])] <- -1
  }else if(sel.case==3) {
    theta[c(Right, Right.cp[i1])] <- 1
    theta[c(Left)] <- -1
  }else if(sel.case==4) {
    theta[c(Right)] <- 1
    theta[c(Left)] <- -1
  }else step()

  list(Elbow.L = Elbow.L, Elbow.R=Elbow.R, Center=Center,
       Right=Right, Left=Left, lambda0=lambda0,
       case.of.lambda = paste("case :",sel.case),
       theta0=theta0, theta = theta)
}
