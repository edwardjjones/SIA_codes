goof2 = function (Observed, Predicted, 
                  coefficient = c("R2", "concordance","MSE", "RMSE", 
                                  "bias", "MSEc", "RMSEc", "RPD", 
                                  "RPIQ"),
                  plot = TRUE, ...) 
{
  if (any(!coefficient %in% c("R2", "concordance", "MSE", "RMSE", 
                              "bias", "MSEc", "RMSEc", "RPD", "RPIQ"))) 
    stop("Please choose a valid coefficient")
  rLM <- lm(Predicted ~ Observed)
  R2 <- as.matrix(summary(rLM)$adj.r.squared)
  SEP2 <- mean((Observed - Predicted)^2)
  SEP <- sqrt(SEP2)
  bias <- mean(Predicted) - mean(Observed)
  SEP2c <- sum(((Predicted - bias - Observed)^2)/length(Observed))
  SEPc <- sqrt(SEP2c)
  RPD <- sd(Observed)/SEP
  IQ <- c(quantile(Observed))[3] - c(quantile(Observed))[2]
  RPIQ <- IQ/SEP
  mx <- mean(Observed)
  my <- mean(Predicted)
  s2x <- var(Observed)
  s2y <- var(Predicted)
  sxy <- cov(Observed,Predicted)
  ccc <- 2 * sxy/(s2x + s2y + (mx - my)^2)
  if (plot) {
    op = par("pty")
    on.exit(par(pty=op))
    par(pty="s")
    lims= c(min(c(Observed, Predicted)),max(c(Observed, Predicted)))
    plot(Observed, Predicted, ylim = lims, xlim = lims, asp = 1,...)
    abline(a = 0, b = 1, col = "grey")
  }
  coefs_tmp <- data.frame(R2 = R2, concordance = ccc, MSE = SEP2, 
                          RMSE = SEP, bias = bias, MSEc = SEP2c, RMSEc = SEPc, 
                          RPD = RPD, RPIQ = RPIQ, row.names = NULL)
  gf <- data.frame(coefs_tmp[, coefficient])
  return(gf)
}
