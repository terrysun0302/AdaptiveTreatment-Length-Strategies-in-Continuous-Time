######################  created by Hao Sun at 08/14/2019 for Supporting Information, Biometrics

#========================================================================================================#
#==================================  functions for the proposed estimators ==============================#
#========================================================================================================#

### The overall function for estimation
#******************************************************************************#
#******************************************************************************#
#                                                                              #
# Inputs                                                                       #
#                                                                              #
#  test.time            a vector of time points for estimation                 #
#                                                                              #
#                                                                              #
#  long.data            the data with n rows, where n is the number of         #
#                       observations.                                          #
#                                                                              #
#                                                                              #
#  stack.data          expand each observation at different time points.       #
#                      see the simulation study for example                    #
#                                                                              #
#  type                the type of outcome models used in the estimation.      #
#                      "gam", generalized additive mode; "plm", partially      #
#                      linear model;"SL", super learner; "lm", linear model;   #
#                      "rf", random forest.                                    #
#******************************************************************************#

mr.fun.disc <- function(test.time, long.data, stack.data, type){
  ########## simulation studies
  
  mean.response.gam <- mean.response.plm <- mean.response.SL <- mean.response.rf <- mean.response.correct <-
    mean.response.gam.tr <- mean.response.plm.tr <- mean.response.SL.tr <-  mean.response.rf.tr <- mean.response.correct.tr <- 
    part1.gam <- part1.plm <- part1.SL <- part1.rf <- part1.correct <- part1.gam.tr <- part1.plm.tr <- part1.SL.tr <- part1.rf.tr <- part1.correct.tr <- 
    part2.gam <- part2.plm <- part2.SL <- part2.rf <- part2.correct <- part2.gam.tr <- part2.plm.tr <- part2.SL.tr <- part2.rf.tr <- part2.correct.tr <- 
    p1.weight <- p1.weight.tr <- c()
  
  num <- length(test.time)
  
  for(t.time in test.time){
    
    part1.mod <- mr.part1.disc(t.time, long.data, stack.data, type)
    part2.mod <- mr.part2.disc(t.time, long.data, stack.data, type)
    
    mean.response.gam <- c(mean.response.gam, part1.mod$p1.gam + part2.mod$p2.gam)
    mean.response.plm <- c(mean.response.plm, part1.mod$p1.plm + part2.mod$p2.plm)
    mean.response.SL <- c(mean.response.SL, part1.mod$p1.SL + part2.mod$p2.SL)
    mean.response.rf <- c(mean.response.rf, part1.mod$p1.rf + part2.mod$p2.rf)
    mean.response.correct <- c(mean.response.correct, part1.mod$p1.correct + part2.mod$p2.correct)
    
    # mean.response.gam.tr <- c(mean.response.gam.tr, part1.mod$p1.gam.tr + part2.mod$p2.gam.tr)
    # mean.response.plm.tr <- c(mean.response.plm.tr, part1.mod$p1.plm.tr + part2.mod$p2.plm.tr)
    # mean.response.SL.tr <- c(mean.response.SL.tr, part1.mod$p1.SL.tr + part2.mod$p2.SL.tr)
    # mean.response.rf.tr <- c(mean.response.rf.tr, part1.mod$p1.rf.tr + part2.mod$p2.rf.tr)
    # mean.response.correct.tr <- c(mean.response.correct.tr, part1.mod$p1.correct.tr + part2.mod$p2.correct.tr)
    
    part1.gam <- c(part1.gam, part1.mod$p1.gam)
    part1.plm <- c(part1.plm, part1.mod$p1.plm)
    part1.correct <- c(part1.correct, part1.mod$p1.correct)
    part1.SL <- c(part1.SL, part1.mod$p1.SL)
    part1.rf <- c(part1.rf, part1.mod$p1.rf)
    p1.weight <- c(p1.weight, part1.mod$p1.weight)
    
    # part1.gam.tr <- c(part1.gam.tr, part1.mod$p1.gam.tr)
    # part1.plm.tr <- c(part1.plm.tr, part1.mod$p1.plm.tr)
    # part1.correct.tr <- c(part1.correct.tr, part1.mod$p1.correct.tr)
    # part1.SL.tr <- c(part1.SL.tr, part1.mod$p1.SL.tr)
    # part1.rf.tr <- c(part1.rf.tr, part1.mod$p1.rf.tr)
    # p1.weight.tr <- c(p1.weight.tr, part1.mod$p1.weight.tr)
    
    part2.gam <- c(part2.gam,  part2.mod$p2.gam)
    part2.plm <- c(part2.plm,  part2.mod$p2.plm)
    part2.correct <- c(part2.correct, part2.mod$p2.correct)
    part2.SL <- c(part2.SL,  part2.mod$p2.SL)
    part2.rf <- c(part2.rf,  part2.mod$p2.rf)
    
    # part2.gam.tr <- c(part2.gam.tr,  part2.mod$p2.gam.tr)
    # part2.plm.tr <- c(part2.plm.tr,  part2.mod$p2.plm.tr)
    # part2.correct.tr <- c(part2.correct.tr, part2.mod$p2.correct.tr)
    # part2.SL.tr <- c(part2.SL.tr,  part2.mod$p2.SL.tr)
    # part2.rf.tr <- c(part2.rf.tr,  part2.mod$p2.rf.tr)
  }
  
  gam <- cbind(test.time, part1.gam, part2.gam,mean.response.gam)
  plm <- cbind(test.time, part1.plm, part2.plm,mean.response.plm)
  SL <- cbind(test.time, part1.SL, part2.SL, mean.response.SL)
  rf <- cbind(test.time, part1.rf, part2.rf, mean.response.rf)
  correct <- cbind(test.time, part1.correct, part2.correct, mean.response.correct)
  
  # gam.tr <- cbind(test.time, part1.gam.tr, part2.gam.tr,mean.response.gam.tr)
  # plm.tr <- cbind(test.time, part1.plm.tr, part2.plm.tr,mean.response.plm.tr)
  # SL.tr <- cbind(test.time, part1.SL.tr, part2.SL.tr, mean.response.SL.tr)
  # rf.tr <- cbind(test.time, part1.rf.tr, part2.rf.tr, mean.response.rf.tr)
  # correct.tr <- cbind(test.time, part1.correct.tr, part2.correct.tr, mean.response.correct.tr)
  
  return(list(gam = gam, plm = plm, SL = SL, rf = rf, correct = correct, p1.weight = p1.weight))
} 

### The first part of the overall estimation. 
#******************************************************************************#
#******************************************************************************#
#                                                                              #
# Inputs                                                                       #
#                                                                              #
#  t.time               a vector of time points for estimation                 #
#                                                                              #
#                                                                              #
#  long.data            the data with n rows, where n is the number of         #
#                       observations.                                          #
#                                                                              #
#                                                                              #
#  stack.data          expand each observation at different time points.       #
#                      see the simulation study for example                    #
#                                                                              #
#  type                the type of outcome models used in the estimation.      #
#                      "gam", generalized additive mode; "plm", partially      #
#                      linear model;"SL", super learner; "lm", linear model;   #
#                      "rf", random forest.                                    #
#******************************************************************************#
mr.part1.disc <- function(t.time, long.data, stack.data,type){
  n <- dim(long.data)[1]
  num.interval <- sum(new.L < t.time)-1
  #original.id <- long.data$id
  #long.data$id <- rep(1:length(original.id))
  ######## effective id, I(U >= t.time)
  long.id <- long.data$id[long.data$tstop > new.L[num.interval+1]] # the effective subjects are changed
  
  w <- long.data$w ############ when w = 1, it is the sample size
  w.long <- w[long.data$id %in% long.id] ############# the eligible w
  
  ############################ the data used for outcome regression
  effectivelong.data <- lapply(long.id, function(x) stack.data[stack.data$id == x,][num.interval+1,])
  effectivelong.data <- do.call(rbind,effectivelong.data)
  effectivelong.data$status.t <- 1
  effectivelong.data$status.c <- 0
  effectivelong.data$tstop <- t.time
  effectivelong.data$tstatus.t <- t.time
  effectivelong.data$tstatus.c <- 0
  #lp <- predict.Gam(gam.binary, effectivelong.data)
  #lp <- predict.Gam(gam.binary.t, effectivelong.data)
  #OR <- exp(lp)/(1 + exp(lp))
  #OR.t.joint <- predict.Gam(gam.linear, effectivelong.data)
  OR.gam <- OR.plm <- OR.SL <- OR.lm <- OR.rf <- c()
  
  #=============================================================================================#
  #======== in practice, the index here should be modified based on the design matrix ==========#
  #=============================================================================================#
  if("gam" %in% type) {OR.gam <- predict(gam.linear.t, effectivelong.data[,3:9], type = "response")}
  if("plm" %in% type) {OR.plm <- predict(plm.t, newdata = effectivelong.data[,3:9], type = "response")}
  if("SL" %in% type) {OR.SL <- predict(SL.t, effectivelong.data[,3:9], onlySL = T, type = "response")$pred}
  if("rf" %in% type) {OR.rf <- predict(rf.t, new_data = effectivelong.data[,3:9], type = "response")}
  if("lm" %in% type) {OR.lm <- predict(lm.t, effectivelong.data, type = "response")} ### correctly specified model
  
  ############################ the data used for cox model
  effectivestack.data <- stack.data[(stack.data$id %in% long.id) & (stack.data$tstart < t.time) ,]
  
  base.c.haz.interval <-  base.u.haz.interval <- base.t.haz.interval <- c()
  if(num.interval > 0){
    for(tt in 1:num.interval){
      base.c.haz.interval <- c(base.c.haz.interval, tail(base.c[base.c$time <= new.L[tt+1],1], n = 1))
      #base.u.haz.interval <- c(base.u.haz.interval, tail(base.u[base.u$time <= new.L[tt+1],1], n = 1))
      base.t.haz.interval <- c(base.t.haz.interval, tail(base.t[base.t$time <= new.L[tt+1],1], n = 1))
    }
    ################ here, we consider $P(C > t.time)$
    base.c.haz.interval <- c(base.c.haz.interval, tail(base.c[base.c$time <= t.time,1], n = 1))
    ################ here because we consider $P(u >= t.time)$
    #base.u.haz.interval <- c(base.u.haz.interval, tail(base.u[base.u$time < t.time,1], n = 1))
    base.t.haz.interval <- c(base.t.haz.interval, tail(base.t[base.t$time < t.time,1], n = 1))
    base.c.haz.interval.diff <- c(base.c.haz.interval[1], diff(base.c.haz.interval))
    #base.u.haz.interval.diff <- c(base.u.haz.interval[1], diff(base.u.haz.interval))
    base.t.haz.interval.diff <- c(base.t.haz.interval[1], diff(base.t.haz.interval))
  } else{
    base.c.haz.interval <- base.c.haz.interval.diff <- tail(base.c[base.c$time <= t.time,1], n = 1)
    #base.u.haz.interval <- base.u.haz.interval.diff <- tail(base.u[base.u$time < t.time,1], n = 1)
    base.t.haz.interval <- base.t.haz.interval.diff <- tail(base.t[base.t$time < t.time,1], n = 1)
  }
  
  ###############################  this part is used in the continuous set up ###################
  
  #=============================================================================================#
  #======== in practice, the index here should be modified based on the design matrix ==========#
  #=============================================================================================#
  
  effectivestack.data.X <- effectivestack.data[,4:9]
  effectivestack.data.c.hazard <- exp(as.matrix(effectivestack.data.X)  %*%  gamma.hat) * base.c.haz.interval.diff
  effectivestack.data.c.cumuhazard <- unname(tapply(effectivestack.data.c.hazard, 
                                                    (seq_along(effectivestack.data.c.hazard)-1)%/%(num.interval + 1), sum))
  #effectivestack.data.u.hazard <- exp(as.matrix(effectivestack.data.X)  %*%  xi.hat) * base.u.haz.interval.diff
  #effectivestack.data.u.cumuhazard <- unname(tapply(effectivestack.data.u.hazard, 
  #                                                 (seq_along(effectivestack.data.u.hazard)-1)%/%(num.interval + 1), sum))
  effectivestack.data.t.hazard <- exp(as.matrix(effectivestack.data.X)  %*%  eta.hat) * base.t.haz.interval.diff
  effectivestack.data.t.cumuhazard <- unname(tapply(effectivestack.data.t.hazard, 
                                                    (seq_along(effectivestack.data.t.hazard)-1)%/%(num.interval + 1), sum))
  
  ############################## this part is used as the weight function in k-stage ###########################
  weight.c.hazard <- effectivestack.data.c.hazard
  weight.c.hazard[(num.interval+1) * (1:length(long.id))] <- 0
  weight.c.cumuhazard <- unname(tapply(weight.c.hazard, 
                                       (seq_along(weight.c.hazard)-1)%/%(num.interval + 1), sum))
  
  weight.t.hazard <- effectivestack.data.t.hazard
  weight.t.hazard[(num.interval+1) * (1:length(long.id))] <- 0
  weight.t.cumuhazard <- unname(tapply(weight.t.hazard, 
                                       (seq_along(weight.t.hazard)-1)%/%(num.interval + 1), sum))
  
  ############## the censoring prob in estimand ##################
  Cen.sur <- as.vector(exp(-effectivestack.data.c.cumuhazard))
  #U.sur <- as.vector(exp(-effectivestack.data.u.cumuhazard))
  #T.sur <- as.vector(exp(-effectivestack.data.t.cumuhazard))
  
  #trueweight <- weightmat[long.id, num.interval+1]
  weight <- as.vector(exp(-weight.c.cumuhazard)) * as.vector(exp(-weight.t.cumuhazard))
  
  ################# the first part of the mean outcome: use the estimated weight function
  p1.gam <- sum(w.long * OR.gam * Cen.sur * 1/weight)/sum(w)
  p1.plm <- sum(w.long * OR.plm * Cen.sur * 1/weight)/sum(w)
  p1.SL <- sum(w.long * OR.SL * Cen.sur * 1/weight)/sum(w)
  p1.rf <- sum(w.long * OR.rf * Cen.sur * 1/weight)/sum(w)
  p1.correct <- sum(w.long * OR.lm * Cen.sur * 1/weight)/sum(w)
  p1.weight <- sum(1/weight)/n # the mean weight function
  
  # ################ check by using true weight
  # p1.gam.tr <- sum(w.long * OR.gam * Cen.sur * 1/trueweight)/sum(w)
  # p1.plm.tr <- sum(w.long * OR.plm * Cen.sur * 1/trueweight)/sum(w)
  # p1.SL.tr <- sum(w.long * OR.SL * Cen.sur * 1/trueweight)/sum(w)
  # p1.rf.tr <- sum(w.long * OR.rf * Cen.sur * 1/trueweight)/sum(w)
  # p1.correct.tr <- sum(w.long * OR.lm * Cen.sur * 1/trueweight)/sum(w)
  # p1.weight.tr <- sum(1/trueweight)/N # the mean weight function
  
  return(list("p1.gam" = p1.gam, 
              "p1.plm" = p1.plm, 
              "p1.SL" = p1.SL,
              "p1.rf" = p1.rf,
              "p1.correct" = p1.correct,
              "p1.weight" = p1.weight))
  
}

### The second part of the overall estimation. 
#******************************************************************************#
#******************************************************************************#
#                                                                              #
# Inputs                                                                       #
#                                                                              #
#  t.time               a vector of time points for estimation                 #
#                                                                              #
#                                                                              #
#  long.data            the data with n rows, where n is the number of         #
#                       observations.                                          #
#                                                                              #
#                                                                              #
#  stack.data          expand each observation at different time points.       #
#                      see the simulation study for example                    #
#                                                                              #
#  type                the type of outcome models used in the estimation.      #
#                      "gam", generalized additive mode; "plm", partially      #
#                      linear model;"SL", super learner; "lm", linear model;   #
#                      "rf", random forest.                                    #
#******************************************************************************#
############ a cleaner version of part2 mean function by using lapply : but is not faster #######################
########### update 4/18/2018: faster 2 times now; to be improved ################################################

mr.part2.disc <- function(t.time, long.data, stack.data, type){
  n <- dim(long.data)[1]
  w <- long.data$w ############ when w = 1, it is the sample size
  
  ################# part 2
  ############ the potential change points of hazard functions: the events time points and covariates change points
  ############ we only consider these points since dG(u) ! = 0 only at these points
  po.change <- unique(pmin(sort(unique(c(long.data$tstop[long.data$status.c == 1], new.L[-1]))), t.time))
  
  ### the base hazard for time covariates interval
  base_po_ind <- apply(as.matrix(po.change),1,function(x) {which(base.c[,2] > x)[1] - 1})
  
  c.base_po <- base.c[base_po_ind,1]
  #u.base_po <- base.u[base_po_ind,1]
  t.base_po <- base.t[base_po_ind,1]
  
  c.base_po_incre <- c(c.base_po[1], diff(c.base_po))
  #u.base_po_incre <- c(u.base_po[1], diff(u.base_po))
  t.base_po_incre <- c(t.base_po[1], diff(t.base_po))
  
  p2.joint <- p2.sepa <- p2.correct <- c()
  po.index <- rep(0, length(po.change))
  for(l in 1:length(new.L)){
    po.index <- po.index + (new.L[l] < po.change)
  }
  
  num.interval <- unlist(lapply(1:n, function(x) nrow(stack.data[stack.data$id == long.data$id[x],])-1))
  tr.length <- pmin(t.time, stop.L[num.interval+1]) # here the trt length can be extrapolated
  
  jumps <- lapply(1:n, function(x) po.change <= tr.length[x])
  ### the covariate value at each jumps
  
  po.idx <- lapply(1:n, function(x){po.index[jumps[[x]] == 1]})
  cov.po <- lapply(1:n, function(x) {Zi <- stack.data[stack.data$id == long.data$id[x],4:9];
  Zi[po.idx[[x]],]})
  
  temp.data <- lapply(1:n, function(x)
  {data.frame(cbind(tstop = po.change[jumps[[x]]], as.matrix(cov.po[[x]]),
                    status.c = rep(1,length(po.change[jumps[[x]]])),
                    status.t = rep(0,length(po.change[jumps[[x]]])),tstatus.c = po.change[jumps[[x]]],
                    tstatus.t = rep(0,length(po.change[jumps[[x]]]))))})
  
  temp.data <- do.call("rbind", temp.data)
  OR.po.gam <- OR.po.plm <- OR.po.SL <- OR.po.correct <- OR.po.rf <- c()
  
  #=============================================================================================#
  #======== in practice, the index here should be modified based on the design matrix ==========#
  #=============================================================================================#
  
  if("gam" %in% type)  OR.po.gam <- predict(gam.linear.c, temp.data, type = "response")
  if("plm" %in% type)  OR.po.plm <- predict(plm.c, newdata = temp.data, type = "response")
  if("SL" %in% type) OR.po.SL <- predict(SL.c, temp.data[,1:7], onlySL = T, type = "response")$pred
  if("rf" %in% type) OR.po.rf <- predict(rf.c, new_data = temp.data[,1:7], type = "response")
  if("lm" %in% type) OR.po.correct <- predict(lm.c, temp.data[,1:7], type = "response")
  
  incre_c.haz_po <- lapply(1:n, function(x) {exp(as.matrix(cov.po[[x]]) %*% gamma.hat) * c.base_po_incre[jumps[[x]]]})
  #incre_u.haz_po <- lapply(1:n, function(x) {exp(as.matrix(cov.po[[x]]) %*% xi.hat) * u.base_po_incre[jumps[[x]]]})
  incre_t.haz_po <- lapply(1:n, function(x) {exp(as.matrix(cov.po[[x]]) %*% eta.hat) * t.base_po_incre[jumps[[x]]]})
  c.surv_po <- lapply(1:n, function(x) exp(-cumsum(incre_c.haz_po[[x]])))
  #u.surv_po <- lapply(1:n, function(x) exp(-cumsum(incre_u.haz_po[[x]])))
  t.surv_po <- lapply(1:n, function(x) exp(-cumsum(incre_t.haz_po[[x]])))
  incre_c.surv_po <- lapply(1:n, function(x) c(c.surv_po[[x]][1]-1, diff(c.surv_po[[x]])))
  
  ################## find the denominator part of the weight function ####################
  ################# the index of endpoints (changepoints)
  w.idx <- lapply(1:n, function(x) {unlist(lapply(1:(num.interval[x]+1), function(y) {which(po.change >= new.L[y])[1]}))})
  w.c.sur <- lapply(1:n, function(x) {a <- c.surv_po[[x]][w.idx[[x]]]; a[1] <- 1;a})
  w.t.sur <- lapply(1:n, function(x) {a <- t.surv_po[[x]][w.idx[[x]]]; a[1] <- 1;a})
  w.sur <- lapply(1:n, function(x){w.c.sur[[x]] * w.t.sur[[x]]})
  weight <- lapply(1:n, function(x) {w.sur[[x]][po.idx[[x]]]/w[x]})
  
  p2.gam <- - sum(OR.po.gam*unlist(incre_c.surv_po)*1/unlist(weight))/sum(w)
  p2.plm <- - sum(OR.po.plm*unlist(incre_c.surv_po)*1/unlist(weight))/sum(w)
  p2.SL <- - sum(OR.po.SL*unlist(incre_c.surv_po)*1/unlist(weight))/sum(w)
  p2.rf <- - sum(OR.po.rf*unlist(incre_c.surv_po)*1/unlist(weight))/sum(w)
  p2.correct <- - sum(OR.po.correct*unlist(incre_c.surv_po)*1/unlist(weight))/sum(w)
  
  
  # trueweight <- lapply(1:n, function(x) weightmat[x,po.idx[[x]]])
  # p2.gam.tr <- - sum(OR.po.gam*unlist(incre_c.surv_po)*1/unlist(trueweight))/sum(w)
  # p2.plm.tr <- - sum(OR.po.plm*unlist(incre_c.surv_po)*1/unlist(trueweight))/sum(w)
  # p2.SL.tr <- - sum(OR.po.SL*unlist(incre_c.surv_po)*1/unlist(trueweight))/sum(w)
  # p2.rf.tr <- - sum(OR.po.rf*unlist(incre_c.surv_po)*1/unlist(trueweight))/sum(w)
  # p2.correct.tr <- - sum(OR.po.correct*unlist(incre_c.surv_po)*1/unlist(trueweight))/sum(w)
  
  return(list("p2.gam" = p2.gam,
              "p2.plm" = p2.plm,
              "p2.SL" = p2.SL,
              "p2.rf" = p2.rf,
              "p2.correct" = p2.correct))
}


#========================================================================================================#
#==================================  functions for the Lu and Johnson(2015)==============================#
#========================================================================================================#

####################################################  the mean response ############################################
mr.part1.ind <- function(t.time, base.data){
  N <- nrow(base.data)
  w <- base.data$w
  
  num <- nrow(base.data)
  long.id <- base.data$id[base.data$tstop >= t.time]
  
  effectivelong.data <- base.data
  effectivelong.data$status.t <- 1
  effectivelong.data$status.c <- 0
  effectivelong.data$tstop <- t.time
  effectivelong.data$tstatus.t <- t.time
  effectivelong.data$tstatus.c <- 0
  effectivelong.data$dur <- t.time
  #lp <- predict.Gam(gam.binary, effectivelong.data)
  #lp <- predict.Gam(gam.binary.t, effectivelong.data)
  #OR <- exp(lp)/(1 + exp(lp))
  OR.t <- predict(lm.t, effectivelong.data, type = "response")
  OR.sepa <- predict(gam.linear.t, effectivelong.data, type = "response")
  OR.SL <- predict(SL.t, effectivelong.data, onlySL = T, type = "response")$pred
  
  base.c.haz <- base.c[tail(which(base.c[,2] <= t.time), n = 1),1]
  
  effectivestack.data.X <- effectivelong.data[,4:9]
  effectivestack.data.c.cumuhazard <- exp(as.matrix(effectivestack.data.X)  %*%  gamma.hat) * base.c.haz
  
  Cen.sur <- as.vector(exp(-effectivestack.data.c.cumuhazard))
  
  ################# the first part of the mean outcome
  p1.sepa <- sum(w*OR.sepa * Cen.sur)/sum(w)
  p1.cor <- sum(w*OR.t * Cen.sur)/sum(w)
  p1.SL <- sum(w*OR.SL * Cen.sur)/sum(w)
  
  return(list(p1.sepa = p1.sepa, p1.correct = p1.cor, p1.SL = p1.SL))
}

mr.part2.ind <- function(t.time, long.data){
  N <- dim(long.data)[1]
  w <- long.data$w ############ when w = 1, it is the sample size
  
  ################# part 2
  ############ the potential change points of hazard functions: the events time points and covariates change points
  ############ we only consider these points since dG(u) ! = 0 only at these points
  po.change <- unique(pmin(sort(unique(c(long.data$tstop[long.data$status.c == 1], new.L[-1]))), t.time))
  
  ### the base hazard for time covariates interval
  base_po_ind <- apply(as.matrix(po.change),1,function(x) {which(base.c[,2] > x)[1] - 1})
  
  c.base_po <- base.c[base_po_ind,1]
  t.base_po <- base.t[base_po_ind,1]
  
  c.base_po_incre <- c(c.base_po[1], diff(c.base_po))
  t.base_po_incre <- c(t.base_po[1], diff(t.base_po))
  
  p2.sepa <- p2.correct <- c()
  
  #num.interval <- unlist(lapply(1:n, function(x) nrow(stack.data[stack.data$id == long.data$id[x],])-1))
  #num.interval <- unlist(lapply(1:n, function(x) nrow(stack.data[stack.data$id == x,])-1))
  
  cov.po <- lapply(1:n, function(x) {Zi <- base.data[x,4:9];
  Zi[rep(1, length(po.change)),]})
  
  temp.data <- lapply(1:n, function(x)
  {data.frame(cbind(tstop = po.change, cov.po[[x]],
                    status.c = 1,
                    status.t = 0, tstatus.c = po.change,
                    tstatus.t = 0))})
  temp.data <- do.call("rbind", temp.data)
  OR.po.sepa <- predict.Gam(gam.linear.c, temp.data,type = "response")
  OR.po.correct <- predict.Gam(lm.c, temp.data, type = "response")
  OR.po.SL <- predict(SL.c, temp.data[,1:7], onlySL = T, type = "response")$pred
  
  incre_c.haz_po <- lapply(1:n, function(x) {as.numeric(exp(as.numeric(base.data[x,4:9]) %*% gamma.hat)) * c.base_po_incre})
  #incre_u.haz_po <- lapply(1:n, function(x) {exp(as.matrix(cov.po[[x]]) %*% xi.hat) * u.base_po_incre[jumps[[x]]]})
  incre_t.haz_po <- lapply(1:n, function(x) {as.numeric(exp(as.numeric(base.data[x,4:9]) %*% eta.hat)) * t.base_po_incre})
  c.surv_po <- lapply(1:n, function(x) exp(-cumsum(incre_c.haz_po[[x]])))
  #u.surv_po <- lapply(1:n, function(x) exp(-cumsum(incre_u.haz_po[[x]])))
  t.surv_po <- lapply(1:n, function(x) exp(-cumsum(incre_t.haz_po[[x]])))
  incre_c.surv_po <- lapply(1:n, function(x) c(c.surv_po[[x]][1]-1, diff(c.surv_po[[x]])))
  
  weight <- rep(w, times = length(po.change))
  p2.sepa <- - sum(OR.po.sepa*unlist(incre_c.surv_po))/sum(w)
  p2.correct <- - sum(OR.po.correct*unlist(incre_c.surv_po))/sum(w)
  p2.SL <- - sum(OR.po.SL*unlist(incre_c.surv_po))/sum(w)
  
  return(list("p2.sepa" = p2.sepa,
              "p2.correct" = p2.correct,"p2.SL" = p2.SL))
}

mr.fun.ind <- function(test.time, long.data){
  ########## simulation studies
  
  #mean.response.SL <- mean.response.gam <- mean.response.correct <- 
  #  part1.gam <- part1.SL <- part1.correct <- part2.gam <- part2.SL <- part2.correct <- p1.weight <- c()
  
  mean.response.joint <- mean.response.sepa <- mean.response.correct <- mean.response.SL <- 
    part1.sepa <- part1.correct <- part1.SL <- part2.sepa <- part2.correct <- part2.SL <- p1.weight <- c()
  
  num <- length(test.time)
  
  for(t.time in test.time){
    
    part1.mod <- mr.part1.ind(t.time, long.data)
    part2.mod <- mr.part2.ind(t.time, long.data)
    
    #mean.response.joint <- c(mean.response.joint, part1.mod$p1.joint + part2.mod$p2.joint)
    mean.response.sepa <- c(mean.response.sepa, part1.mod$p1.sepa + part2.mod$p2.sepa)
    mean.response.correct <- c(mean.response.correct, part1.mod$p1.correct + part2.mod$p2.correct)
    mean.response.SL <- c(mean.response.SL, part1.mod$p1.SL + part2.mod$p2.SL)
    
    #part1.joint <- c(part1.joint, part1.mod$p1.joint)
    part1.sepa <- c(part1.sepa, part1.mod$p1.sepa)
    part1.correct <- c(part1.correct, part1.mod$p1.correct)
    part1.SL <- c(part1.SL, part1.mod$p1.SL)
    
    #part2.joint <- c(part2.joint, part2.mod$p2.joint)
    part2.sepa <- c(part2.sepa,  part2.mod$p2.sepa)
    part2.correct <- c(part2.correct, part2.mod$p2.correct)
    part2.SL <- c(part2.SL, part2.mod$p2.SL)
    
    #p1.weight <- c(p1.weight, part1.mod$p1.weight)
  }
  
  #joint <- cbind(test.time, part1.joint, part2.joint,mean.response.joint)
  sepa <- cbind(test.time, part1.sepa, part2.sepa,mean.response.sepa)
  correct <- cbind(test.time, part1.correct, part2.correct, mean.response.correct)
  SL <- cbind(test.time, part1.SL, part2.SL, mean.response.SL)
  
  return(list(sepa = sepa,correct = correct, SL=SL))
} 