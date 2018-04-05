#' @rdname svrpath-internal
#' @export
re.Initialization <- function(x,y,Left, Right, Center, svr.eps= 1, lambda = lambda[k],
                               kernel.function = poly.kernel,param.kernel = 1){
  Left.cp <- Left
  Right.cp <- Right
  lambdak <- lambda
  Elbow.L <- Elbow.R <- NULL
  eps <- svr.eps
  mod <- 1e-5
  K <- kernel.function(x, x, param.kernel = param.kernel)

  gx <- apply(K[,Right, drop = F], 1, sum) - apply(K[,Left, drop = F], 1, sum)
  lambda1 <- (outer(gx[Left], gx[Right], "-")) / (outer(y[Left], y[Right], "-") + 2 * eps)
  lambda2 <- (outer(gx[Left], gx[Center], "-")) / (outer(y[Left], y[Center], "-"))
  lambda3 <- (outer(gx[Right], gx[Center], "-")) / (outer(y[Right], y[Center], "-"))
  lambda4 <- (outer(gx[Center], gx[Center], "-")) / (outer(y[Center], y[Center], "-") + 2 * eps)

  l1 <- lambda1[lambda1 < lambdak - mod]
  l2 <- lambda2[lambda2 < lambdak - mod]
  l3 <- lambda3[lambda3 < lambdak - mod]
  l4 <- lambda4[lambda4 < lambdak - mod]

  if(length(l1) == 0){ l1 <- 0 }
  if(length(l2) == 0){ l2 <- 0 }
  if(length(l3) == 0){ l3 <- 0 }
  if(length(l4) == 0){ l4 <- 0 }

  max.lambda1 <- max(l1)
  max.lambda2 <- max(l2)
  max.lambda3 <- max(l3)
  max.lambda4 <- max(l4)
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

  if(ncol(as.matrix(temp)) == 1 && nrow(as.matrix(temp)) == 1){
    i1 <- 1
    i2 <- 1
  }else{
    i1 <- which(temp == lambda0, arr.ind = T)[,1]
    i2 <- which(temp == lambda0, arr.ind = T)[,2]
  }

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

  k <- NULL
  list(Elbow.L = Elbow.L, Elbow.R=Elbow.R, Center=Center,
       Right=Right, Left=Left, lambda0=lambda0,
       case.of.lambda = paste("case :",sel.case),
       theta0 = theta0, theta = theta)
}
