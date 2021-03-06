% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/int_mat.R
\name{int_mat}
\alias{int_mat}
\title{Numerical Integration Matrix}
\usage{
int_mat(
  distr = dnorm,
  args = list(mean = 0, sd = 1),
  lb = -4,
  ub = 4,
  npts = 10000
)
}
\arguments{
\item{distr}{A density function with two user-specified parameters. Defaults
to the normal distribution (dnorm), but any density function is permitted.}

\item{args}{Named list of arguments to distr.}

\item{lb}{Lower bound of range over which to numerically integrate.}

\item{ub}{Upper bound of range over which to numerically integrate.}

\item{npts}{Number of integration points.}
}
\value{
Matrix of two columns. Column 1 is a sequence of x-coordinates, and
column 2 is a sequence of y-coordinates from a normalized distribution.
}
\description{
Create a matrix for numerical integration.
}
\seealso{
\link{rimse} \link{th_est_ml} \link{th_est_eap} \link{sl_link}
\link{hb_link}

 @importFrom stats dnorm
}
