#================================================================================================================#
#======= this is an example of using the proposed method with a pseudo data set =================================#
#======= you can e-mail brent_johnson@urmc.rochester.edu for the infusion trial data set ========================#
#================================================================================================================#



#=============================================================#
#========== Step 1: laod data, functions and packages=========#
#=============================================================#

path = setwd("~/Desktop/Supporting Information/")
source(paste0(path,"/R example/MeanFuncDiscInfusion.r",sep=""))
load(paste0(path,"/R data/long.data",".rda",sep=""))
load(paste0(path,"/R data/stack.data",".rda",sep=""))
load(paste0(path,"/R data/base.data",".rda",sep=""))

require(survival)
require(gam)
require(randomForest)
require("SuperLearner")
require("earth")
require("gam")
require("ranger")
require(KernSmooth)
#library(np)
require("hal9001")




#==============================================================#
#=============== Step 2: fit different outcome models,=========#
#============== which will be inserted into the estimation ====#
#==============================================================#

### the number of reasampling in se estimation
K <- 150
### the landmarks
new.L <- c(0,6,12,18,24)
stop.L <- c(6,12,18,24,30)
### the time points of interest
test.time <- c(16,18,20,22,24)

### some arries to save the resampling results
NPB.td.gam <- NPB.td.correct <- NPB.td.SL <- array(numeric(0), dim = c(length(test.time),3,K))
NPB.ti.gam <- NPB.ti.correct <- NPB.ti.SL <- array(numeric(0), dim = c(length(test.time),3,K))


#=========================================#
#=== step 2.1: fit the Q, G functions ====#
#=========================================#
######### fit the assignment process model
c.cox <- coxph(Surv(tstart, tstop, status.c) ~ time.cov+dbts+ptca+angina+wght+heparin,data=stack.data)
gamma.hat <- coefficients(c.cox)
### the cumulative baseline hazard function
base.c <- basehaz(c.cox, centered = F)

########## time dependent cox model for U = min(T, C)
t.cox <- coxph(Surv(tstart, tstop, status.t) ~ time.cov+dbts+ptca+angina+wght+heparin,data=stack.data)
eta.hat <- coefficients(t.cox)
### the cumulative baseline hazard function
base.t <- basehaz(t.cox, centered = F)



#=========================================#
#=== step 2.2: fit the m1, m0 functions ==#
#=========================================#

#########  the outcome regression #########
### gam
gam.linear.t <- gam(outcome ~ s(tstop) + time.cov + dbts + ptca + angina + wght + heparin,
                    data = long.data[long.data$status.c == 0, ], family=binomial(link='logit'))
gam.linear.c <- gam(outcome ~ s(tstop) + time.cov + dbts + ptca + angina + wght + heparin,
                    data = long.data[long.data$status.c == 1,], family=binomial(link='logit'))
### parametric 
lm.t <- glm(outcome ~ log(tstop) + time.cov + dbts + ptca + angina + wght + heparin,
            data = long.data[long.data$status.t == 1, ], family=binomial(link='logit'))
lm.c <- glm(outcome ~ log(tstop) + time.cov + dbts + ptca + angina + wght + heparin,
            data = long.data[long.data$status.c == 1,], family=binomial(link='logit'))

### super learner
# "SL.earth","SL.glm.interaction","SL.mean","SL.ranger","SL.randomForest"
sl.lib=c("SL.gam","SL.glm", "SL.earth","SL.glm.interaction","SL.mean","SL.randomForest")
xt <- long.data[long.data$status.c == 0, 3:9]
xc <- long.data[long.data$status.c == 1, 3:9]
outcomet <- long.data[long.data$status.c == 0,]$outcome
outcomec <- long.data[long.data$status.c == 1,]$outcome
SL.t <- SuperLearner(Y=outcomet, X=xt, SL.library=sl.lib, family=binomial(link='logit'))
SL.c <- SuperLearner(Y=outcomec, X=xc, SL.library=sl.lib, family=binomial(link='logit'))


#==============================================================#
#=============== Step 3: the proposed estimation ==============#
#==============================================================#
# the proposed estimator
mr.mod.td <- mr.fun.disc(test.time, long.data, stack.data, type = c("lm", "gam", "SL"))

### then we do the bootstrap-based variance estimation
for(k in 1:K){
  set.seed(k)
  print(k)
  
  idx <- sample(n, replace = T)
  re.long <- long.data[idx,]
  re.long$id <- 1:n
  re.stack <- lapply(1:n, function(x){a <- stack.data[stack.data$id == idx[x],]; a$id <- x; a})
  re.stack <- Reduce(rbind, re.stack)
  re.base <- base.data[idx,]
  
  ########## fit the support functions 
  c.cox <- coxph(Surv(tstart, tstop, status.c) ~ time.cov+dbts + ptca + angina + wght + heparin, data=re.stack)
  gamma.hat <- coefficients(c.cox)
  ### the cumulative baseline hazard function
  base.c <- basehaz(c.cox, centered = F)
  
  t.cox <- coxph(Surv(tstart, tstop, status.t) ~ time.cov+dbts + ptca + angina + wght + heparin, data=re.stack)
  eta.hat <- coefficients(t.cox)
  ### the cumulative baseline hazard function
  base.t <- basehaz(t.cox, centered = F)
  
  
  ################################  the outcome regression #################################################
  ### semiparametric
  gam.linear.t <- gam(outcome ~ s(tstop) + time.cov + dbts + ptca + angina + wght + heparin,
                      data = re.long[re.long$status.c == 0, ], family=binomial(link='logit'))
  gam.linear.c <- gam(outcome ~ s(tstop) + time.cov + dbts + ptca + angina + wght + heparin,
                      data = re.long[re.long$status.c == 1,], family=binomial(link='logit'))
  ### parametric model
  lm.t <- glm(outcome ~ log(tstop) + time.cov + dbts + ptca + angina + wght + heparin,
             data = re.long[re.long$status.c == 0, ], family=binomial(link='logit'))
  lm.c <- glm(outcome ~ log(tstop) + time.cov + dbts + ptca + angina + wght + heparin,
             data = re.long[re.long$status.c == 1,], family=binomial(link='logit'))
  
  ### super learner
  # "SL.earth","SL.glm.interaction","SL.mean","SL.ranger","SL.randomForest"
  sl.lib=c("SL.gam","SL.glm", "SL.earth","SL.glm.interaction","SL.mean","SL.randomForest")
  xt <- re.long[re.long$status.c == 0, 3:9]
  xc <- re.long[re.long$status.c == 1, 3:9]
  outcomet <- re.long[re.long$status.c == 0,]$outcome
  outcomec <- re.long[re.long$status.c == 1,]$outcome
  SL.t <- SuperLearner(Y=outcomet, X=xt, SL.library=sl.lib, family=binomial(link='logit'))
  SL.c <- SuperLearner(Y=outcomec, X=xc, SL.library=sl.lib, family=binomial(link='logit'))
  
  
  N <- nrow(re.long)
  
  ######### estimate the mean
  re.td <- mr.fun.disc(test.time, re.long, re.stack, type = c("lm", "gam", "SL"))
  
  ######### save the result to calculate the resampling variance estimation
  NPB.td.gam[,,k] <- re.td$gam[,-1]
  NPB.td.correct[,,k] <- re.td$correct[,-1]
  NPB.td.SL[,,k] <- re.td$SL[,-1]
}




#==============================================================#
#============== Step 4: compare with Lucy's method ============#
#==============================================================#
### the lucy's estimator

c.cox <- coxph(Surv(tstop, status.c) ~ time.cov + dbts + ptca + angina + wght + heparin,data= base.data)
gamma.hat <- coefficients(c.cox)
### the cumulative baseline hazard function
base.c <- basehaz(c.cox, centered = F)

gam.linear.t <- gam(outcome ~ s(tstop) + time.cov + dbts + ptca + angina + wght + heparin,
                    data = base.data[base.data$status.c == 0, ], family=binomial(link='logit'))
gam.linear.c <- gam(outcome ~ s(tstop) + time.cov + dbts + ptca + angina + wght + heparin,
                    data = base.data[base.data$status.c == 1,], family=binomial(link='logit'))
### parametric 
lm.t <- glm(outcome ~ log(tstop) + time.cov + dbts + ptca + angina + wght + heparin,
            data = base.data[base.data$status.t == 1, ], family=binomial(link='logit'))
lm.c <- glm(outcome ~ log(tstop) + time.cov  + dbts + ptca + angina + wght + heparin,
            data = base.data[base.data$status.c == 1,], family=binomial(link='logit'))

### super learner
# "SL.earth","SL.glm.interaction","SL.mean","SL.ranger","SL.randomForest"
sl.lib=c("SL.gam","SL.glm", "SL.earth","SL.glm.interaction","SL.mean","SL.randomForest")
xt <- base.data[base.data$status.c == 0, 3:9]
xc <- base.data[base.data$status.c == 1, 3:9]
outcomet <- base.data[base.data$status.c == 0,]$outcome
outcomec <- base.data[base.data$status.c == 1,]$outcome
SL.t <- SuperLearner(Y=outcomet, X=xt, SL.library=sl.lib, family=binomial(link='logit'))
SL.c <- SuperLearner(Y=outcomec, X=xc, SL.library=sl.lib, family=binomial(link='logit'))

mr.mod.ti <- mr.fun.ind(test.time, base.data)

### the bootstrap based variance estimation
for(k in 1:K){
  set.seed(k)
  print(k)
  
  idx <- sample(n, replace = T)
  re.base <- base.data[idx,]
  
  #################################################################################
  ##################################  NPB part  ###################################
  #################################################################################
  ########## permute the data and refit all support functions
  ########## fit the support functions 
  #############################################################################
  c.cox <- coxph(Surv(tstop, status.c) ~ time.cov + dbts + ptca + angina + wght + heparin,data= re.base)
  gamma.hat <- coefficients(c.cox)
  ### the cumulative baseline hazard function
  base.c <- basehaz(c.cox, centered = F)
  
  #################################  the outcome regression ###################
  ### the semiparametric model
  gam.linear.t <- gam(outcome ~ s(tstop) + time.cov + dbts + ptca + angina + wght + heparin,
                      data = re.base[re.base$status.c == 0, ], family=binomial(link='logit'))
  gam.linear.c <- gam(outcome ~ s(tstop) + time.cov + dbts + ptca + angina + wght + heparin,
                      data = re.base[re.base$status.c == 1,], family=binomial(link='logit'))
  
  #### the parametric model
  lm.t <- glm(outcome ~ log(tstop) + time.cov + dbts + ptca + angina + wght + heparin,
              data = re.base[re.base$status.c == 0, ], family=binomial(link='logit'))
  lm.c <- glm(outcome ~ log(tstop) + time.cov + dbts + ptca + angina + wght + heparin,
              data = re.base[re.base$status.c == 1,], family=binomial(link='logit'))
  
  ### super learner
  # "SL.earth","SL.glm.interaction","SL.mean","SL.ranger","SL.randomForest"
  sl.lib=c("SL.gam","SL.glm", "SL.earth","SL.glm.interaction","SL.mean","SL.randomForest")
  xt <- re.base[re.base$status.c == 0, 3:9]
  xc <- re.base[re.base$status.c == 1, 3:9]
  outcomet <- re.base[re.base$status.c == 0,]$outcome
  outcomec <- re.base[re.base$status.c == 1,]$outcome
  SL.t <- SuperLearner(Y=outcomet, X=xt, SL.library=sl.lib, family=binomial(link='logit'))
  SL.c <- SuperLearner(Y=outcomec, X=xc, SL.library=sl.lib, family=binomial(link='logit'))
  N <- nrow(re.base)  
  
  ######### estimate the mean
  re.ti <- mr.fun.ind(test.time, re.base)
  
  ######### save the result to calculate the resampling variance estimation
  NPB.ti.gam[,,k] <- re.ti$sepa[,-1]
  NPB.ti.correct[,,k] <- re.ti$correct[,-1]
  NPB.ti.SL[,,k] <- re.ti$SL[,-1]
}


#==============================================================#
#========= Step 5: summarize and save the results =============#
#==============================================================#
################################################## summarize the results ########################################
mean.gam.td <- mr.mod.td$gam[,-1]
mean.correct.td <- mr.mod.td$correct[,-1]
mean.SL.td <- mr.mod.td$SL[,-1]
mean.gam.ti <- mr.mod.ti$sepa[,-1]
mean.correct.ti <- mr.mod.ti$correct[,-1]
mean.SL.ti <- mr.mod.ti$SL[,-1]

############################# save the results when time.cov is mdoeled in the OR #################################
mean.tv <- cbind(mean.SL.td[,3], mean.gam.td[,3], mean.correct.td[,3], mean.SL.ti[,3], mean.gam.ti[,3], mean.correct.ti[,3])
colnames(mean.tv) <- c("SL(Proposed)","PLM(proposed)", "LM(Proposed)","SL(LJ)", "PLM(LJ)", "LM(LJ)")
rownames(mean.tv) <- c("16", "18", "20", "22", "24")
sd.tv <- cbind(apply(NPB.td.SL, c(1,2),sd)[,3], apply(NPB.td.gam, c(1,2),sd)[,3], apply(NPB.td.correct, c(1,2),sd)[,3], 
               apply(NPB.ti.SL, c(1,2),sd)[,3], apply(NPB.ti.gam, c(1,2),sd)[,3], apply(NPB.ti.correct, c(1,2),sd)[,3])


############################# save the results when time.cov is not mdoeled in the OR #############################
mean.fix <- cbind(mean.SL.td[,3], mean.gam.td[,3], mean.correct.td[,3], 
                  mean.SL.ti[,3], mean.gam.ti[,3], mean.correct.ti[,3])
colnames(mean.fix) <- colnames(mean.tv) <- c("SL(Proposed)","PLM(proposed)", "LM(Proposed)","SL(LJ)", "PLM(LJ)", "LM(LJ)")
rownames(mean.fix) <- c("16", "18", "20", "22", "24")
sd.fix <- cbind(apply(NPB.td.SL, c(1,2),sd)[,3], apply(NPB.td.gam, c(1,2),sd)[,3], apply(NPB.td.correct, c(1,2),sd)[,3], 
                apply(NPB.ti.SL, c(1,2),sd)[,3], apply(NPB.ti.gam, c(1,2),sd)[,3], apply(NPB.ti.correct, c(1,2),sd)[,3])


realdata <- list(mean.tv = mean.tv, sd.tv = sd.tv, mean.fix = mean.fix, sd.fix = sd.fix)
save(realdata, file = "realdata.rda")

