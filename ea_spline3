##For rasters
##Convert to matrix structure (x,y,d0_d1,d1_d2,d2_d3,d3_d4)
##Output = matrix (x,y,0-1,1-2,...,)

ea_spline3=function (obj, lam = 0.1, d = t(0:90), vlow = 0, vhigh = 1000){
  
  
  
  ##maxdepth
  mxd <- max(d)
  ##variance of input data
  s2 <- (0.05 * mean(obj[, 3:ncol(obj)]))^2
  
  nl <- length(lam)
  yave <- rep(NA, length(d)+1)
  
  message("Fitting mass preserving splines per pixel...")
  
  n <- length(3:ncol(obj))
  
  np1 <- n + 1
  nm1 <- n - 1
  
  
  ##get upper and lower depths
  depths=str_split(colnames(obj)[3:ncol(obj)],"_")
  ##upper depths
  u=as.numeric(lapply(depths,FUN=function(x)x[1]))
  ##lower depths
  v=as.numeric(lapply(depths,FUN=function(x)x[2]))
  
  delta <- t(v - u)
  
  del <- c(u[2:n], u[n]) - v
  
  r <- matrix(0, ncol = nm1, nrow = nm1)
  for (dig in 1:nm1) {
    r[dig, dig] <- 1
  }
  for (udig in 1:nm1 - 1) {
    r[udig, udig + 1] <- 1
  }
  d2 <- matrix(0, ncol = nm1, nrow = nm1)
  diag(d2) <- delta[2:n]
  r <- d2 %*% r
  r <- r + t(r)
  d1 <- matrix(0, ncol = nm1, nrow = nm1)
  diag(d1) <- delta[1:nm1]
  d3 <- matrix(0, ncol = nm1, nrow = nm1)
  diag(d3) <- del[1:nm1]
  r <- r + 2 * d1 + 6 * d3
  q <- matrix(0, ncol = n, nrow = n)
  for (dig in 1:n) {
    q[dig, dig] <- -1
  }
  for (udig in 1:n - 1) {
    q[udig, udig + 1] <- 1
  }
  q <- q[1:nm1, 1:n]
  dim.mat <- matrix(q[], ncol = length(1:n), nrow = length(1:nm1))
  rinv <- try(solve(r), TRUE)
  ind <- diag(n)
  pr.mat <- matrix(0, ncol = length(1:nm1), nrow = length(1:n))
  pr.mat[] <- 6 * n * lam
  fdub <- pr.mat * t(dim.mat) %*% rinv
  z <- fdub %*% dim.mat + ind
  nj <- max(v)
  if (nj > mxd) {
    nj <- mxd
  }
  xfit <- as.matrix(t(c(1:mxd)))
  yfit <- xfit
  
  
  ea_spline_single=function(x, rinv,dim.mat, sbar, delta, mxd, nj, xfit, yfit, u, v) {
    
    # if (is.matrix(rinv)) {
    sbar <- solve(z, (x))
    b <- 6 * rinv %*% dim.mat %*% sbar
    b0 <- rbind(0, b)
    b1 <- rbind(b, 0)
    gamma <- (b1 - b0)/t(2 * delta)
    alfa <- sbar - b0 * t(delta)/2 - gamma * t(delta)^2/3
    for (k in 1:nj) {
      xd <- xfit[k]
      if (xd < u[1]) {
        p <- alfa[1]
      }
      else {
        for (its in 1:n) {
          if (its < n) {
            tf2 = as.numeric(xd > v[its] & xd < u[its + 
                                                    1])
          }
          else {
            tf2 <- 0
          }
          if (xd >= u[its] & xd <= v[its]) {
            p = alfa[its] + b0[its] * (xd - u[its]) + 
              gamma[its] * (xd - u[its])^2
          }
          else if (tf2) {
            phi = alfa[its + 1] - b1[its] * (u[its + 
                                                 1] - v[its])
            p = phi + b1[its] * (xd - v[its])
          }
        }
      }
      yfit[k] = p
    }
    if (nj < mxd) {
      yfit[, (nj + 1):mxd] = NA
    }
    yfit[which(yfit > vhigh)] <- vhigh
    yfit[which(yfit < vlow)] <- vlow
    
    if(min(u)!=0){
      yfit[,1:min(u)]= NA
    }
    
    # m_fyfit[st, ] <- yfit
    nd <- length(d) - 1
    dl <- d + 1
    for (cj in 1:nd) {
      xd1 <- dl[cj]
      xd2 <- dl[cj + 1] - 1
      if (nj >= xd1 & nj <= xd2) {
        xd2 <- nj - 1
        yave[cj] <- mean(yfit[xd1:xd2])
      }
      else {
        yave[ cj] <- mean(yfit[xd1:xd2])
      }
      # yave[cj + 1] <- min(u)
      # yave[cj + 2] <- max(v)
    }
    
    
    # ssq <- sum((t(y) - sbar)^2)
    # g <- solve(z)
    # ei <- eigen(g)
    # ei <- ei$values
    # df <- n - sum(ei)
    # sig2w <- ssq/df
    # dfc <- n - 2 * sum(ei) + sum(ei^2)
    # sig2c <- ssq/dfc
    # tmse <- ssq/n - 2 * s2 * df/n + s2
    # sset <- tmse
    # }
    
    return(yave)
  }
  # if (show.progress) {
  #   setTxtProgressBar(pb, st)
  # }
  # }
  
  out = apply(obj[,3:ncol(obj)],1,function(x) ea_spline_single(x, rinv,dim.mat, sbar, delta, mxd, nj, xfit, yfit, u, v))
  
  out=t(out[1:(length(d)-1),])
  
  colname_v = paste0("d",d[-length(d)],"_",d[-1])
  
  colnames(out)=colname_v
  
  out=cbind(obj[,1:2],out)
  
  
  # for(m in 1:10){
  #   plot(y=0:89,x=dave[1:90,m],ylim=c(89,0))
  #   for(n in 1:4){
  #     lines(x=c(obj[m,n+2],obj[m,n+2]),y=c(u[n],v[n]),lwd=2,col="red")
  #   }
  # }
  
  
 
  
  # if (show.progress) {
  #   close(pb)
  # }
  # yave <- as.data.frame(yave)
  # jmat <- matrix(NA, ncol = 1, nrow = length(d))
  # for (i in 1:length(d) - 1) {
  #   a1 <- paste(d[i], d[i + 1], sep = "-")
  #   a1 <- paste(a1, "cm", sep = " ")
  #   jmat[i] <- a1
  # }
  # jmat[length(d)] <- "min depth"
  # jmat[length(d)+1] <- "max depth"
  # for (jj in 1:length(jmat)) {
  #   names(yave)[jj] <- jmat[jj]
  # }
  # yave <- cbind(mat_id[, 1], yave)
  # names(yave)[1] <- "id"
  # retval <- list(harmonised = yave, obs.preds = dave, var.1cm = t(m_fyfit), 
  #                tmse = sset)
  message("...Finished processing.")
  return(out)
}
