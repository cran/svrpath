loss <- function(w, y, eps = .1) {
  r <- (y - w)
  v <- (r-eps) * (r > eps) + (-eps - r) * (r < -eps)
  sum(v)
}
