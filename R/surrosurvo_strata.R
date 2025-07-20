#' @importFrom pimeta cima
surrosurvo_strata <- function(y, event, x, eventx, strata,
                              censtype = c("univariate", "independent"),
                              level = 0.95, B = 25000, parallel = FALSE,
                              seed = NULL) {

  out <- NULL
  L <- max(strata)
  for (i in 1:L) {

    # by strata
    y2 <- y[strata == i]
    event2 <- event[strata == i]
    x2 <- x[strata == i]
    if (is.null(eventx)) {
      eventx2 <- NULL
    } else {
      eventx2 <- eventx[strata == i]
    }
    outs <- surrosurvo_base(y = y2, event = event2, x = x2, eventx = eventx2,
                            censtype = censtype, level = level, parallel = parallel)
    out <- rbind(
      out,
      data.frame(strata = i, outs)
    )

  }

  out1 <- out[out$method == "taumo1",]
  out2 <- out[out$method == "taumo2",]
  out3 <- out[out$method == "tauo",]
  out4 <- out[out$method == "tauso",]

  cfit1 <- cima(out1$tau, out1$se, alpha = 1 - level, method = "boot",
                B = B, seed = seed, parallel = parallel)
  cfit2 <- cima(out2$tau, out2$se, alpha = 1 - level, method = "boot",
                B = B, seed = seed, parallel = parallel)
  cfit3 <- cima(out3$tau, out3$se, alpha = 1 - level, method = "boot",
                B = B, seed = seed, parallel = parallel)
  cfit4 <- cima(out4$tau, out4$se, alpha = 1 - level, method = "boot",
                B = B, seed = seed, parallel = parallel)

  cout <- rbind(
    data.frame(method = "taumo1", tau = cfit1$muhat,
               lcl = cfit1$lci, ucl = cfit1$uci),
    data.frame(method = "taumo2", tau = cfit2$muhat,
               lcl = cfit2$lci, ucl = cfit2$uci),
    data.frame(method = "tauo", tau = cfit3$muhat,
               lcl = cfit3$lci, ucl = cfit3$uci),
    data.frame(method = "tauso", tau = cfit4$muhat,
               lcl = cfit4$lci, ucl = cfit4$uci)
  )

  return(cout)

}
