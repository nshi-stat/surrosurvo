#' Evaluation of Individual-Level Surrogacy for Survival and Ordinal Endpoints
#'
#' This function can estimate Kendall's \eqn{\tau} as follows:
#' Oakes (1982)'s non-parametric estimator \eqn{\hat{\tau}_o},
#' Lakhal et al. (2009)'s inverse probability weighting (IPCW)
#' estimators \eqn{\hat{\tau}_{mo1}} and \eqn{\hat{\tau}_{mo2}},
#' and Okui et al. (2024)'s modified IPCW estimators
#' \eqn{\hat{\tau}_{so1}} and \eqn{\hat{\tau}_{so2}}.
#' \eqn{\hat{\tau}_{so1}} and \eqn{\hat{\tau}_{so2}} can take
#' into account tie data. Confidence intervals were based on
#' jackknife standard errors with a normal approximation.
#'
#' @name surrosurvo
#' @rdname surrosurvo
#' @param y the survival time or censoring time outcome vector
#' @param event the event indicator outcome vector
#' @param x the ordinal (or continuous) outcome vector
#' @param level the confidence level of the confidence interval
#' (default = 0.95)
#' @param confint whether the calculation of confidence
#' intervals is performed (default = TRUE)
#' @param parallel the number of threads used in parallel
#' computing, or FALSE that means single threading
#' (default = FALSE)
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
#' Okui, J., Nagashima, K., Matsuda, S., et al (2024).
#' Surrogacy between pathological complete response and overall survival:
#' An individual patient data analysis of ten RCTs on neoadjuvant treatment
#' plus surgery for esophageal cancer.
#' \emph{In preparation}.
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
surrosurvo <- function(y, event, x, level = 0.95, confint = TRUE, parallel = FALSE) {

  # initial check
  util_check_ge(y, 0)
  util_check_inrange(event, 0, 1)
  util_check_gt(x, 0)
  util_check_inrange(level, 0.0, 1.0)

  # estimating S_c
  n <- length(y)
  time <- y
  censor <- ifelse(event == 0, 1, 0)
  ctime <- ifelse(censor == 1, time, Inf)
  df <- data.frame(id = 1:n, time, ctime, censor, x)
  fit <- survfit(Surv(time, censor) ~ 1, df)
  sc <- data.frame(time = fit$time, surv = fit$surv)
  df <- merge(df, sc, by = "time")
  df <- df[order(df$id),]
  df[, c("id", "censor")] <- NULL

  # point estimates
  res <- surrosurvoCpp(df, n)
  out <- data.frame(method = colnames(res), tau = t(res))
  rownames(out) <- NULL

  # confidence intervals
  if (confint) {
    if (parallel == FALSE) {
      parallel <- 1
    }
    ci <- confintCpp(df, n, parallel)
    # jackknife se
    se <- data.frame(method = colnames(res), se = c(
      sqrt((n - 1)/n*sum((ci$jktauo - res$tauo)^2)),
      sqrt((n - 1)/n*sum((ci$jktaumo1 - res$taumo1)^2)),
      sqrt((n - 1)/n*sum((ci$jktaumo2 - res$taumo2)^2)),
      sqrt((n - 1)/n*sum((ci$jktauso1 - res$tauso1)^2)),
      sqrt((n - 1)/n*sum((ci$jktauso2 - res$tauso2)^2))
    ))
    alpha <- 1 - level
    out <- merge(out, se, by = "method")
    out$lcl <- out$tau + qnorm(alpha/2)*out$se
    out$ucl <- out$tau + qnorm(1 - alpha/2)*out$se
    out$lcl[out$lcl < -1] <- -1
    out$lcl[out$lcl > 1] <- 1
    out$ucl[out$ucl < -1] <- -1
    out$ucl[out$ucl > 1] <- 1
  }

  return(out)

}
