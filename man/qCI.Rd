% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/procFunctions.R
\name{qCI}
\alias{qCI}
\title{(internal) Reformatting bootstrapped quantile estimates into the CIs}
\usage{
qCI(lenQ, Qmat, alpha, nB, n)
}
\arguments{
\item{lenQ}{an integer, length of Q}

\item{Qmat}{a matrix, matrix of bootstrapped quantile estimates}

\item{alpha}{a float, the size of the test}

\item{nB}{an integer, number of bootstrap iterations}

\item{n}{an integer, sample size}
}
\value{
a matrix of three columns of CIs with mean estimates in the middle.
}
\description{
a function reformat matrix of quantile estimates into the uniform CIs
}
