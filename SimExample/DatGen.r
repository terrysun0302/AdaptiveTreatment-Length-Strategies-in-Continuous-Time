############## created by Hao Sun at 3/24/2018
############## functions to generate the data and estimate the true potential outcomes over time

####################################################################################################################
####################################################################################################################
#########################  function to generate the true potential outcome #########################################
####################################################################################################################
####################################################################################################################
true.potential <- function(time, base.cen, sp, repl){
  
  #################  firstly generate the full dat
  u.cen <- runif(repl)
  
  #beta.trt <- c(-0.5, 0.5)
  #gamma.cen <- c(1, -0.3)
  #beta.trt <- c(-0.3, 0.5) # the one used for 45% and ill defined 20%
  #beta.trt <- c(-0.1, 0.5) # the one used for well defined 20%
  gamma.cen <- c(0.3, -0.5)

  
  
  ################## parameter for outcome model #########################
  inter <- 1
  #alpha <- c(0.5, -0.2)
  alpha <- c(0.5, -0.5)

  upp <- 30
  
  new.L <- c(0,6,12,18,24)
  stop.L <- c(6,12,18,24,30)
  
  
  ################## generate the full data ##############################
  t.cov0 <- rnorm(repl, mean = 0, sd = 1)
  t.cov6 <- rnorm(repl, mean = t.cov0, sd = 1)
  t.cov12 <- rnorm(repl, mean = t.cov6, sd = 1)
  t.cov18 <- rnorm(repl, mean =  t.cov12, sd = 1)
  t.cov24 <- rnorm(repl, mean = t.cov18, sd = 1)
  t.cov <- cbind(t.cov0,t.cov6,t.cov12,t.cov18,t.cov24)
  t.cov <- as.vector(t(t.cov))
  
  id <- rep(1:repl, each = 5)
  tstart <- rep(new.L, time = repl)
  tstop <- rep(stop.L, time = repl)
  
  #dbts <- rbinom(n, 1, prob = 0.4)
  dbts <- rnorm(repl)
  # ptca <- rbinom(n, 1, prob = 0.25)
  # angina <- rbinom(n,1,prob = 0.4)
  # weight <- rgamma(n, shape = 28, scale = 3)
  # heparin <- rbinom(n, 1, prob = 0.15)
  
  full.data <- data.frame(id = id, tstart = tstart, tstop = tstop, time.cov = t.cov, dbts = rep(dbts, each = 5))
  
  tlength <- length(time) ### the number of test time points
  
  Delta <- matrix(numeric(0), nrow = tlength, ncol = repl) ### can be sued to calculate P(C > t)
  y <- c()
  
  for(i in 1:repl){
    
    time.data <- full.data[full.data$id == i,]
    x <- full.data[full.data$id == i,-c(1:3)]
    
    #ad.event <- cen.time(u.cen[i], x, base.cen, gamma.cen, new.L)
    ad.event <- trt.time.even(u.cen[i], x, gamma.cen, new.L, basec1, basec2, basec3, basec4, basec5) ### the bell shape dist
    
    delta <- time < ad.event ## a vector with lengtht
    
    #error <- rnorm(sum(delta), sd = 0.7)
    #error <- c(error, rep(rnorm(1), length(test.time) - length(error)))
    dur <- pmin(time, ad.event)
    dur <- pmin(dur, upp)
    
    Delta[,i] <- delta
    outcome.idx <- lapply(1:tlength, function(x) tail(which(dur[x] > time.data[,2]), n = 1))
    outcome.data <- time.data[unlist(outcome.idx),]
    
    outcome <- inter + delta*sqrt(dur) + alpha[1] * 2 * exp(outcome.data[4]) + 
      alpha[2] * 5 * (outcome.data[4] + outcome.data[5])^2 - 
      delta + log(outcome.data[3]) * (1-delta)
    
    y <- cbind(y, as.matrix(outcome))
  }
  return(cbind(test.time = time, mean.res = apply(y,1,sum), 
              part1 = apply((y*Delta), 1, sum), p2 = apply((y*(1-Delta)), 1, sum)))
}

######################################################################################################################
######################################################################################################################
########################################### function to simulate the data ############################################
######################################################################################################################
######################################################################################################################

dat.gen <- function(n, base.cen, sp, type = "gaussian"){
  
  ################### parameter for cox models ###########################
  #beta.trt <- c(-0.5, 0.5)
  #gamma.cen <- c(1, -0.3)
  #beta.trt <- c(-0.3, 0.5) # the one used for 45% and ill defined 20%
  beta.trt <- c(-0.1, 0.5) # the one used for well defined 20%
  gamma.cen <- c(0.3, -0.5)
  
  ################## parameter for outcome model #########################
  inter <- 1
  #alpha <- c(0.5, -0.2)
  alpha <- c(1, -0.5)

  
  upp <- 30
  
  new.L <- c(0,6,12,18,24)
  stop.L <- c(6,12,18,24,30)
  
  
  ################## generate the full data ##############################
  t.cov0 <- rnorm(n, mean = 0, sd = 1)
  t.cov6 <- rnorm(n, mean = t.cov0, sd = 1)
  t.cov12 <- rnorm(n, mean = t.cov6, sd = 1)
  t.cov18 <- rnorm(n, mean =  t.cov12, sd = 1)
  t.cov24 <- rnorm(n, mean = t.cov18, sd = 1)
  t.cov <- cbind(t.cov0,t.cov6,t.cov12,t.cov18,t.cov24)
  t.cov <- as.vector(t(t.cov))
  
  id <- rep(1:n, each = 5)
  tstart <- rep(new.L, time = 5)
  tstop <- rep(stop.L, time = 5)
  
  #dbts <- rbinom(n, 1, prob = 0.4)
  dbts <- rnorm(n)
  # ptca <- rbinom(n, 1, prob = 0.25)
  # angina <- rbinom(n,1,prob = 0.4)
  # weight <- rgamma(n, shape = 28, scale = 3)
  # heparin <- rbinom(n, 1, prob = 0.15)
  
  full.data <- data.frame(id = id, tstart = tstart, tstop = tstop, time.cov = t.cov, dbts = rep(dbts, each = 5))
  
  
  ################################## construct the observed data  #########################################
  u.trt <- runif(n)
  u.cen <- runif(n)
  epsilon <- rnorm(n, sd = 0.7)
  long.data <- stack.data <- base.data <- c()
  
  
  for(i in 1:n){
    time.data <- full.data[full.data$id == i,]
    x <- full.data[full.data$id == i,-c(1:3)]
    
    #trt.length <- trt.time.bell(u.trt[i], x, beta.trt, new.L, sp = sp) ### the bell shape dist
    trt.length <- trt.time.even(u.trt[i], x, beta.trt, new.L, baset1, baset2, baset3, baset4, baset5) ### the bell shape dist
    #trt.length <- cen.time(u.trt[i], x, base.trt, beta.trt, new.L)
    #trt.length <-  cen.time(u.trt[i], time.data[,-c(1,2)], 0.06, beta.trt, new.L)
    ad.event <- trt.time.even(u.cen[i], x, gamma.cen, new.L, basec1, basec2, basec3, basec4, basec5) ### the bell shape dist
    #ad.event <- cen.time.bell(u.cen[i], x, gamma.cen, new.L, sp = sp2) ### the bell shape dist
    if(trt.length > upp & ad.event > upp){
      status <- 0 # mean the subject is censored by end of study
      delta <- 1
      dur <- upp
    } else {
      status <- 1
      delta <- trt.length < ad.event
      dur <- ifelse(delta == 1, trt.length, ad.event)
      dur <- ifelse(dur >= upp, upp, dur) 
    }
    
    if(dur == 0) next
    
    real.data <- time.data[dur > new.L,]
    if(!is.null(nrow(real.data))){
      tstart <- real.data[,2]
      tstop <- c(real.data[-1,2],dur)
      
      status.c <- c(rep(0, length(tstart)-1), status*(1-delta))
      status.t <- c(rep(0, length(tstart)-1), status*delta)
      status.u <- status.c + status.t
      
      #real.data <- cbind(real.data[,1], tstart, tstop, real.data[,-c(1,2)], status.c, status.u, status.t)
      real.data$tstop <- tstop
      real.data$status.c <- status.c
      real.data$status.u <- status.u
      real.data$status.t <- status.t
      real.data$delta <- c(rep(0, length(tstart)-1), delta)
      
      outcome.data <- as.matrix(real.data[nrow(real.data),])
      baseline.data <- outcome.data
      baseline.data[4:5] <- as.numeric(real.data[1,4:5])
    }
    else{
      tstart <- 0
      tstop <- dur
      
      status.c <- c(rep(0, length(tstart)-1), status*(1-delta))
      status.t <- c(rep(0, length(tstart)-1), status*(delta))
      status.u <- status.c + status.t
      
      real.data$tstop <- tstop
      real.data$status.c <- status.c
      real.data$status.u <- status.u
      real.data$status.t <- status.t
      real.data$delta <- c(rep(0, length(tstart)-1), delta)
      
      outcome.data <- as.matrix(real.data)
      baseline.data <- outcome.data
    }
    
    if(type == "gaussian"){
      outcome <- inter + delta*sqrt(outcome.data[3]) + alpha %*% as.matrix(outcome.data[4:5]) - delta +
        log(outcome.data[3]) * (1-delta) + epsilon[i]
      # outcome <- inter + delta*sqrt(dur) + alpha[1] * 2 * exp(outcome.data[4]) +
      #   alpha[2] * 2 * (outcome.data[4] + outcome.data[5])^2 -
      #   delta + log(outcome.data[3]) * (1-delta) + epsilon[i]
    }
    else if(type == "binary"){
        lp <- delta*0.3*sqrt(outcome.data[3]) + alpha %*% as.matrix(outcome.data[4:5]) - delta +
        0.3*log(outcome.data[3]) * (1-delta)
      outcome <- rbinom(length(lp),1,prob = exp(lp)/(1 + exp(lp)))
    }

    
    outcome.data <- c(outcome.data, (1-delta)*outcome.data[3], delta*outcome.data[3], outcome)
    baseline.data <- c(baseline.data, (1-delta)*baseline.data[3], delta*baseline.data[3], outcome)
    
    base.data <- rbind(base.data, baseline.data)
    long.data <- rbind(long.data, outcome.data)
    stack.data <- rbind(stack.data, real.data)
  }
  colnames(long.data) <- colnames(base.data) <- c("id", "tstart", "tstop","time.cov", "dbts",
                           "status.c", "status.u", "status.t", "delta","tstatus.c","tstatus.t","outcome")
  colnames(stack.data) <- c("id", "tstart", "tstop","time.cov", "dbts",
                            "status.c", "status.u", "status.t","delta")
  full.data <- as.data.frame(full.data) ############# the data set used for cox model
  base.data <- as.data.frame(base.data)
  long.data <- as.data.frame(long.data)
  stack.data <- as.data.frame(stack.data)
  
  return(list(full.data = full.data, long.data = long.data, stack.data = stack.data, base.data = base.data))
}


####################################################################################################################
####################################################################################################################
####################################function to test the prob estimation ###########################################
####################################################################################################################
####################################################################################################################

true.cen <- function(base.cen, time.data, t, gamma){
  obs.data <- time.data[time.data$tstart < t, ]
  obs.data$tstop[length(obs.data$tstop)] <- t
  
  base <- base.cen * (obs.data$tstop - obs.data$tstart)
  lp <- as.matrix(obs.data[,4:5]) %*% gamma
  
  cen.sur <- exp(-sum(exp(lp) * base))
  
  return(cen.sur)
}

true.meanreg <- function(time.data, t, alpha){
  # dur <- data$tstop
  # delta <- data$status.t
  # cova <- as.matrix(data[,4:9])
  # output <- inter - delta + delta*sqrt(dur) + cova %*% alpha  +
  #   (1-delta)*log(dur)
  
  outcome.data <- time.data[tail(which(t > time.data[,2]), n = 1),]
  
  dur <- t
  #delta <- outcome.data$delta
  cova <- as.matrix(outcome.data[,4:5])
  output <- sqrt(dur) + cova %*% alpha
  
  return(output)
}


true.u <- function(t,i){
       id <- long.data$id[i]
       time.data <- as.matrix(full.data[full.data$id == id,])
       return(test.potential(t, time.data))
}

test.potential <- function(t, time.data){
  repl <- 2000
  delta.c <- c()
  delta.t <- c()
  delta.u <- c()
  
  u.trt <- runif(repl)
  u.cen <- runif(repl)
  
  x <- time.data[,-c(1:3)]
  
  for(i in 1:repl){
    trt.length <- trt.time.bell(u.trt[i], x, beta.trt, new.L, sp = sp) ### the bell shape dist
    #trt.length <- cen.time(u.trt[i], x, base.trt, beta.trt, new.L)
    #trt.length <-  cen.time(u.trt[i], time.data[,-c(1,2)], 0.06, beta.trt, new.L)
    ad.event <- cen.time(u.cen[i], x, base.cen, gamma.cen, new.L)
    #ad.event <- cen.time.bell(u.cen[i], x, gamma.cen, new.L, sp = sp2) ### the bell shape dist
    
    delta.c <- c(delta.c, ad.event < t)
    delta.t <- c(delta.t, trt.length < t)
    
    dur <- min(trt.length, ad.event)
    
    delta.u <- c(delta.u, dur < t) 
  }
  
  return(list(cen.sur = 1 - mean(delta.c), trt.sur = 1 - mean(delta.t), u.sur = 1 - mean(delta.u)))
}