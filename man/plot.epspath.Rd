% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.epspath.R
\name{plot.epspath}
\alias{plot.epspath}
\title{plot the epspath, solution paths of SVR as a function of epsilon}
\usage{
\method{plot}{epspath}(x, intercept = FALSE, ...)
}
\arguments{
\item{x}{The epspath object}

\item{intercept}{If it is \code{TRUE}, then an intercept path plot is given.}

\item{...}{Generic compatibility}
}
\value{
The entire solution path of SVR solution as a function of \code{epsilon}.
}
\description{
produces a plot of the SVR \code{epsilon} path.
}
\examples{
\donttest{
# The 'eobj' is given by examples description of epspath().
plot(eobj, lty = 2, lwd = 2, col = 2, cex.lab = 1.5) }
}
\author{
Do Hyun Kim, Seung Jun Shin
}
