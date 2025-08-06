#' Random number generation for Correlated Survival and Survival Outcomes
#'
#' This function generates random number vectors of correlated survival and
#' survival outcomes. Gaussian or Clayton copula is used to specify the
#' correlation structure. Marginal distributions are based on a exponential
#' distribution.
#'
#' @param target_tau the target value of the Kendall's \eqn{\tau}
#' @param num_samples the number of simulation samples
#' @param censoring_rate the censoring rate for \code{y}
#' @param censoring_ratex the censoring rate for \code{x}
#' @param censtype the type of censoring model for survival and survival
#'   outcomes (default = "univariate"; see Lakhal et al. (2009)).
#' \itemize{
#' \item \code{univariate}: univariate censoring model:
#'   \eqn{Y = \min(T, C)}, \eqn{X = \min(T_X, C)}
#' \item \code{independent}: independent censoring model:
#'   \eqn{Y = \min(T, C)}, \eqn{X = \min(T_X, C_X)}
#' }
#' @param copula_type copula type (default = "Gaussian")
#' \itemize{
#' \item \code{Gaussian}: Gaussian copula
#' \item \code{Clayton}: Clayton copula
#' }
#' @param hr the hazard rate for survival time \code{y} (default = 1.0)
#' @param hrx the hazard rate for survival time \code{x} (default = 1.0)
#' @return
#' \itemize{
#' \item \code{y}: the survival time or censoring time outcome vector \code{y}
#' \item \code{event}: the event indicator outcome vector \code{y}
#' \item \code{x}: the survival time or censoring time outcome vector \code{x}
#' \item \code{eventx}: the event indicator outcome vector \code{x}
#' \item \code{t}: the true survival time vector \code{y} (for simulation)
#' \item \code{c}: the true censoring time vector \code{y} (for simulation)
#' \item \code{tx}: the true survival time vector \code{x} (for simulation)
#' \item \code{cx}: the true censoring time vector \code{x} (for simulation)
#' }
#' @examples
#' library(surrosurvo)
#' set.seed(1234)
#' data <- generate_srv(0.7, 500, 0.1, 0.1)
#' surrosurvo(data$y, data$event, data$x, data$eventx)
#' @importFrom copula normalCopula claytonCopula rCopula
#' @importFrom stats qexp quantile rexp
#' @export
generate_srv <- function(target_tau, num_samples,
                         censoring_rate, censoring_ratex = NULL,
                         censtype = c("univariate", "independent"),
                         copula_type = c("Gaussian", "Clayton"),
                         hr = 1, hrx = 1) {

  # initial check
  util_check_inrange(target_tau, -1.0, 1.0)
  util_check_gt(num_samples, 1)
  util_check_inrange(censoring_rate, 0, 1.0)
  util_check_gt(hr, 0)
  lstc <- c("Gaussian", "Clayton")
  copula_type <- match.arg(copula_type)
  if (!is.element(copula_type, lstc)) {
    stop("Unknown 'copula_type' specified.")
  }
  lstc <- c("univariate", "independent")
  censtype <- match.arg(censtype)
  if (!is.element(censtype, lstc)) {
    stop("Unknown 'censtype' specified.")
  }
  if (!is.null(censoring_ratex) | censtype == "independent") {
    util_check_inrange(censoring_ratex, 0, 1.0)
  }

  # generate random number vectors
  if (copula_type == "Gaussian") {
    cop <- normalCopula(param = sin(target_tau * pi / 2), dim = 2)
  } else if (copula_type == "Clayton") {
    cop <- claytonCopula(param = 2 * target_tau / (1 - target_tau), dim = 2)
  }

  u <- rCopula(num_samples, cop)

  # y
  hrc <- (censoring_rate * hr)/(1 - censoring_rate)
  t <- qexp(u[, 1])/hr
  c <- rexp(num_samples)/hrc
  y <- pmin(t, c)
  event <- ifelse(y == t, 1, 0)

  # x
  tx <- qexp(u[, 2])/hrx
  # tx <- pmin(t, tx)
  if (censtype == "univariate") {
    cx <- c
  } else {
    hrcx <- (censoring_ratex * hrx)/(1 - censoring_ratex)
    cx <- rexp(num_samples)/hrcx
  }
  x <- pmin(tx, cx)
  eventx <- ifelse(x == tx, 1, 0)

  # output
  data.frame(y, event, x, eventx, t, c, tx, cx)

}
