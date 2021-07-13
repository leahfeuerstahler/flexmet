#' Numerical Integration Matrix
#'
#' Create a matrix for numerical integration.
#'
#' @param distr A density function with two user-specified parameters. Defaults
#' to the normal distribution (dnorm), but any density function is permitted.
#' @param args Named list of arguments to distr.
#' @param lb Lower bound of range over which to numerically integrate.
#' @param ub Upper bound of range over which to numerically integrate.
#' @param npts Number of integration points.
#'
#' @return Matrix of two columns. Column 1 is a sequence of x-coordinates, and
#' column 2 is a sequence of y-coordinates from a normalized distribution.
#'
#' @seealso \link{rimse} \link{th_est_ml} \link{th_est_eap} \link{sl_link}
#' \link{hb_link}
#'
#'  @importFrom stats dnorm
#' @export


int_mat <- function(distr = dnorm, args = list(mean = 0, sd = 1),
                    lb = -4, ub = 4, npts = 10000) {

  # take a uniform sequence of points over the given range
  xvals <- seq(lb, ub, length = npts)

  # find the height of the density function at each point
  
  yvals <- do.call(distr, c(list(x = xvals), args))

  # normalize the y values
  yvals <- yvals / sum(yvals)

  # output
  cbind(xvals, yvals)
}
