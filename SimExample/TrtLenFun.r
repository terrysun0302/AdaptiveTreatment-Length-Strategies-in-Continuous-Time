############## created by Hao Sun at 3/24/2018

####################################################################################################################
####################################################################################################################
###########################################  function for treatment length #########################################
####################################################################################################################
####################################################################################################################

surv.cal <- function(lambda, lp, low, t, sp){
  lambda/exp(lp) - 3*pnorm(t, mean = 24, sd = sp) + 3*pnorm(low, mean = 24, sd = sp) - 0.02*(t - low)
}

trt.time.bell <- function(u, x, beta, new.L, sp){ # x is the stacked data, and bh is the baseline hazard
  x <- as.matrix(x)
  
  step <- c()
  N <- nrow(x)
  base.haz <- 3* pnorm(new.L, mean = 24, sd = sp) + 0.02 * new.L 
  base.haz <- base.haz[-length(base.haz)]
  lp <- x %*% beta
  step <- cumsum(exp(lp[-length(lp)]) * c(diff(base.haz),0.20*6))
  
  lambda <- -log(u)
  count <- sum(lambda >= step) + 1
  if(count == 1){
    t <- uniroot(surv.cal, lower = 0, upper = new.L[2], lambda = lambda, lp = lp[1], low = 0, sp  =sp)$root
  }
  else if(count < (length(new.L)-1)){
    t <- uniroot(surv.cal, lower = new.L[count], upper = new.L[count+1], lambda = lambda - step[count-1], lp = lp[count], low = new.L[count], sp = sp)$root
  }
  else{
    t <- (lambda - step[count-1])/0.17/exp(lp[count]) + new.L[count]
  }
  return(t)
}


trt.time.even <- function(u, x, beta, new.L, baset1, baset2, baset3, baset4, baset5){ # x is the stacked data, and bh is the baseline hazard
  x <- as.matrix(x)
  base <- c(baset1, baset2, baset3, baset4, baset5)
  step <- c()
  N <- nrow(x)
  base.haz <- base *6
  base.haz <- base.haz[-length(base.haz)]
  lp <- x %*% beta
  step <- cumsum(exp(lp[-length(lp)]) * base.haz)

  lambda <- -log(u)
  count <- sum(lambda >= step) + 1
  if(count == 1){
    t <- lambda/base[count]/exp(lp[1])
  }
  else if(count < (length(new.L))){
    t <- (lambda - step[count-1])/base[count]/exp(lp[count]) + new.L[count]
  }
  else{
    t <- (lambda - step[count-1])/base[count]/exp(lp[count]) + new.L[count]
  }
  return(t)
}


cen.time <- function(u, x, bh, beta, new.L){ # x is the stacked data, and bh is the baseline hazard
  x <- as.matrix(x)
  
  step <- c()
  N <- nrow(x)
  time.diff <- c(diff(new.L),0) # we don't consider the last inear predictor
  lp <- x %*% beta
  step <- bh * cumsum(exp(lp) * time.diff)[-length(lp)]
  
  lambda <- -log(u)
  count <- sum(lambda >= step) + 1
  if(count == 1){
    t <- lambda/bh/exp(lp[1])
  }
  else if(count < length(new.L)) {
    t <- (lambda - step[count-1])/bh/exp(lp[count]) + new.L[count]
  }
  else {
    t <- (lambda - step[count-1])/bh/exp(lp[count]) + new.L[count]
  }
}

cen.cal <- function(lambda, lp, low, t, sp){
  lambda/exp(lp) - pnorm(t, mean = 24, sd = sp) + pnorm(low, mean = 24, sd = sp) - 0*(t - low)
}

cen.time.bell <- function(u, x, beta, new.L, sp){ # x is the stacked data, and bh is the baseline hazard
  x <- as.matrix(x)
  
  step <- c()
  N <- nrow(x)
  base.haz <- pnorm(new.L, mean = 24, sd = sp) + 0 * new.L
  lp <- x %*% beta
  step <- cumsum(exp(lp[-length(lp)]) * diff(base.haz))
  
  lambda <- -log(u)
  count <- sum(lambda >= step) + 1
  if(count == 1){
    t <- uniroot(cen.cal, lower = 0, upper = new.L[2], lambda = lambda, lp = lp[1], low = 0, sp  =sp)$root
  }
  else if(count < length(new.L)){
    t <- uniroot(cen.cal, lower = new.L[count], upper = new.L[count+1], lambda = lambda - step[count-1], lp = lp[count], low = new.L[count], sp = sp)$root
  }
  else{
    t <- (lambda - step[count-1])/0.15/exp(lp[count]) + new.L[count]
  }
  return(t)
}

true.weight <- function(new.L, full.data){
  ################### parameter for cox models ###########################
  #beta.trt <- c(-0.5, 0.5)
  #gamma.cen <- c(1, -0.3)
  #beta.trt <- c(-0.3, 0.5) # the one used for well defined 45% and ill defined 20%
  beta.trt <- c(-0.1, 0.5) # the one used for well defined 20%
  gamma.cen <- c(0.3, -0.5)
 
    ### new.L is the endpoints, X is the full data
    L <- length(new.L)
    N <- nrow(full.data)/(length(new.L))
    weightmat <- matrix(numeric(0), ncol = L-1, nrow = N)
    for(i in 1:N){
      X <- as.matrix(full.data[((i-1)*5+1) : (i*5),4:5])
      baset <- c(baset1, baset2, baset3, baset4)
      basec <- c(basec1, basec2, basec3, basec4)
      cumhazt <- cumsum(exp(beta.trt %*% t(X[-L,])) %*% diag(baset) * 6)
      cumhazc <- cumsum(exp(gamma.cen %*% t(X[-L,])) %*% diag(basec) * 6)
      surv.t <- exp(-cumhazt)
      surv.c <- exp(-cumhazc) 
      weight <- surv.t*surv.c
      weightmat[i,] <- weight
    }
    weightmat <- cbind(1, weightmat)
  return(weightmat)
}
