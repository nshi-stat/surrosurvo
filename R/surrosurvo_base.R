#' @importFrom survival Surv survfit
#' @importFrom stats qnorm
surrosurvo_base <- function(y, event, x, eventx = NULL, level = 0.95,
                            parallel = FALSE) {

  # initial check
  util_check_ge(y, 0)
  util_check_inrange(event, 0, 1)
  util_check_num(x, 0)
  if (!is.null(eventx)) {
    util_check_inrange(eventx, 0, 1)
  }
  util_check_inrange(level, 0.0, 1.0)

  # estimating S_c & Sx_c
  n <- length(y)
  time <- y
  censor <- ifelse(event == 0, 1, 0)
  ctime <- ifelse(censor == 1, time, Inf)
  if (!is.null(eventx)) {
    censorx <- ifelse(eventx == 0, 1, 0)
    ctimex <- ifelse(censorx == 1, x, Inf)
    df <- data.frame(id = 1:n, time, ctime, censor, x, ctimex, censorx)
    fit1 <- survfit(Surv(time, censor) ~ 1, df)
    fit2 <- survfit(Surv(x, censorx) ~ 1, df)
    sc1 <- data.frame(time = fit1$time, surv = fit1$surv)
    sc2 <- data.frame(x = fit2$time, survx = fit2$surv)
    df <- merge(df, sc1, by = "time")
    df <- merge(df, sc2, by = "x")
    df <- df[order(df$id),]
  } else {
    df <- data.frame(id = 1:n, time, ctime, censor, x, ctimex = Inf)
    fit <- survfit(Surv(time, censor) ~ 1, df)
    sc <- data.frame(time = fit$time, surv = fit$surv)
    df <- merge(df, sc, by = "time")
    df <- df[order(df$id),]
    df$survx <- 1
  }
  df <- df[,c("time", "ctime", "x", "ctimex", "surv", "survx")]

  # point estimates
  res <- surrosurvoCpp(df, n)
  out <- data.frame(method = colnames(res), tau = t(res))
  rownames(out) <- NULL

  # confidence intervals
  if (parallel == FALSE) {
    parallel <- 1
  }
  ci <- confintCpp(df, n, parallel)
  # jackknife se
  se <- data.frame(method = colnames(res), se = c(
    sqrt((n - 1)/n*sum((ci$jktauo - res$tauo)^2)),
    sqrt((n - 1)/n*sum((ci$jktaumo1 - res$taumo1)^2)),
    sqrt((n - 1)/n*sum((ci$jktaumo2 - res$taumo2)^2)),
    sqrt((n - 1)/n*sum((ci$jktauso - res$tauso)^2))
  ))
  alpha <- 1 - level
  out <- merge(out, se, by = "method")
  out$lcl <- out$tau + qnorm(alpha/2)*out$se
  out$ucl <- out$tau + qnorm(1 - alpha/2)*out$se
  out$lcl[out$lcl < -1] <- -1
  out$lcl[out$lcl > 1] <- 1
  out$ucl[out$ucl < -1] <- -1
  out$ucl[out$ucl > 1] <- 1

  return(out)

}
