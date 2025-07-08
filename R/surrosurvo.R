#' Evaluation of Individual-Level Surrogacy for Survival and Ordinal,
#' Continuous, or Survival Endpoints
#'
#' This function can estimate Kendall's \eqn{\tau} as follows:
#' Oakes (1982)'s non-parametric estimator \eqn{\hat{\tau}_o},
#' Lakhal et al. (2009)'s inverse probability weighting (IPCW)
#' estimators \eqn{\hat{\tau}_{mo1}} and \eqn{\hat{\tau}_{mo2}},
#' and Okui et al. (2025)'s modified IPCW estimator
#' \eqn{\hat{\tau}_{so}}. \eqn{\hat{\tau}_{so}} can take
#' into account tie data. Confidence intervals were based on
#' jackknife standard errors with a normal approximation.
#' Moreover, theses extensions of stratified estimators,
#' \eqn{\hat{\tau}_o^S}, \eqn{\hat{\tau}_{mo1}^S},
#' \eqn{\hat{\tau}_{mo1}^S}, \eqn{\hat{\tau}_{so}^S}, are also
#' available. Stratified estimators are based on (random-effects)
#' meta-analysis models.
#'
#' @param y the survival time or censoring time outcome vector
#' @param event the event indicator outcome vector
#' @param x the ordinal, continuous, or survival outcome vector
#' @param eventx the event indicator outcome vector
#' @param strata the stratified variable vector
#' @param level the confidence level of the confidence interval
#' (default = 0.95)
#' @param B the number of bootstrap samples
#' @param parallel the number of threads used in parallel
#' computing, or FALSE that means single threading
#' (default = FALSE)
#' @param seed set the value of random seed
#' @return
#' \itemize{
#' \item \code{method}: names for estimators of \eqn{\hat{\tau}}
#' \item \code{tau}: estimates \eqn{\hat{\tau}}
#' \item \code{se}: standard errors for \eqn{\hat{\tau}}
#' \item \code{lcl}: lower confidence limits \eqn{\hat{\tau}_l}
#' \item \code{ucl}: uper confidence limits \eqn{\hat{\tau}_u}
#' }
#' @references
#' Oakes, D. (1982).
#' A model for association in bivariate survival data.
#' \emph{J R Stat Soc Ser B Methodol.}
#' \strong{44}(3): 414-422.
#' \url{https://doi.org/10.1111/j.2517-6161.1982.tb01222.x}
#'
#' Lakhal, L, Rivest, L-P., and Beaudoin, D. (2009).
#' IPCW estimator for Kendall's tau under bivariate censoring.
#' \emph{Int J Biostat.}
#' \strong{5}(1): Article 8.
#' \url{https://doi.org/10.2202/1557-4679.1121}
#'
#' Okui, J., Nagashima, K., Matsuda, S., et al (2025).
#' Evaluation of pathological complete response as a surrogate
#' endpoint for overall survival in resectable esophageal cancer.
#' \emph{Br J Surg.}
#' \strong{112}(6)): znaf131.
#' \url{https://doi.org/10.1093/bjs/znaf131}
#' @seealso
#' \code{\link[=generate_rv]{generate_rv}}.
#' @examples
#' library(surrosurvo)
#' set.seed(1234)
#' data <- generate_rv(0.7, 300, 0.3)
#' surrosurvo(data$y, data$event, data$x)
#'
#' # Parallel computation
#' \donttest{surrosurvo(data$y, data$event, data$x, parallel = 2)}
#' @importFrom survival Surv survfit
#' @importFrom stats qnorm
#' @export
surrosurvo <- function(y, event, x, eventx = NULL,
                       censtype = c("univariate", "independent"),
                       strata = NULL, level = 0.95, B = 25000,
                       parallel = FALSE, seed = NULL) {

  # initial check
  util_check_ge(y, 0)
  util_check_inrange(event, 0, 1)
  util_check_num(x, 0)
  if (!is.null(eventx)) {
    util_check_inrange(eventx, 0, 1)
  }
  util_check_inrange(level, 0.0, 1.0)
  if (!is.null(strata)) {
    strata <- as.numeric(factor(strata))
    L <- max(strata)
    util_check_gt(L, 1)
  }

  lstc <- c("univariate", "independent")
  censtype <- match.arg(censtype)
  if (!is.element(censtype, lstc)) {
    stop("Unknown 'censtype' specified.")
  }

  # analysis
  if (!is.null(strata)) {
    # stratified
    res <- surrosurvo_strata(y = y, event = event, x = x, eventx = eventx,
                             strata = strata, censtype = censtype,
                             level = level, B = B, parallel = parallel,
                             seed = seed)
  } else {
    # not stratified
    res <- surrosurvo_base(y = y, event = event, x = x, eventx = eventx,
                           censtype = censtype, level = level,
                           parallel = parallel)
  }

  return(res)

}
