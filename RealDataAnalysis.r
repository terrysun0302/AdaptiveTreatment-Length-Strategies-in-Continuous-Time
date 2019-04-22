#================================================================================================================#
#======= this is an example of using the proposed method with a private data set ================================#
#======= you can e-mail brent_johnson@urmc.rochester.edu for the infusion trial data set ========================#
#================================================================================================================#


#=============================================================#
#========== Step 1: data clean for the application===========#
#=============================================================#
### read the data
#source("C:/Users/hsun16/Desktop/HSun/Research Project/Causal inference/Treatment length/data/InfDat.fanli.r")
source("C:/Users/hsun16/Desktop/HSun/Research Project/Causal inference/Treatment length/data/ReadEspritData.r")

### the file for proposed estimator
### we simulate the data based on the real data #################################
############### the function when cosnidering time dependent confounders
source("C:/Users/hsun16/Desktop/HSun/Research Project/Causal inference/Treatment length/data/MeanFuncDiscInfusion.r")

### construct the data set for our functions
library(survival)
library(gam)
############################### since some people don't have time-dependent covariate ############################
new.esprit.dat <- esprit.dat[which((esprit.dat$id %in% unique(Inf.CK.dat$id)) == 1),]

any.event <- ifelse(new.esprit.dat$dth+new.esprit.dat$mi+new.esprit.dat$pci+new.esprit.dat$cabg>0,1,0)

new.Y <- any.event
new.L <- c(0, 6, 12, 18, 24)

total.num <- unique(length(new.esprit.dat$id))
new.id <- new.esprit.dat$id
N <- length(new.id)

stack.data <- c() ### used for proposed estimator
long.data <- c() ### used for proposed estimator
base.data <- c() ### used for lucy's estimator

for(i in 1:total.num){
  id <- new.id[i]
  dat <- new.esprit.dat[new.esprit.dat$id == id,]
  dur <- dat$dur
  cens <- dat$cens
  ### baseline covariate ###
  dbts <- dat$dbts
  ptca <- dat$ptca
  angina <- dat$angina
  wght <- dat$wght
  heparin <- dat$heparin
  ### time dependent covariates ###
  ### some people miss the baseline one ###
  ck.obs <- Inf.CK.dat[Inf.CK.dat$id == id, ]
  ck.obs <- ck.obs[order(ck.obs$time),]
  ck.impute <- matrix(numeric(0), ncol = ncol(ck.obs), nrow = length(new.L))
  miss.part <- new.L %in% ck.obs$time
  
  ck.impute[which(miss.part==1),] <- as.matrix(ck.obs)
  
  
  if(miss.part[1] == FALSE){ # we miss the baseline value of time dependent covariate
    ck.impute[1,] <- c(id = id, time = 0, value = ck.obs$value[1])
  }
  for(m in 2:length(new.L)){
    if(miss.part[m] == FALSE){
      ck.impute[m,] <-  ck.impute[m-1,] 
      ck.impute[m,2] <- new.L[m]  
    }
  }
  
  ar <- new.L <= dur
  ar.time <- new.L[ar]
  time.cov <- ck.impute[ar,3]
  
  tstart <- ar.time
  tstop <- c(ar.time[-1], dur)
  time.data <- cbind(id = id, tstart = tstart, tstop=tstop, time.cov=time.cov, dbts = dbts,
                     ptca=ptca, angina = angina, wght = wght , heparin = heparin)
  if(tail(tstart, n=1) == tail(tstop, n=1)){
    time.data <- head(time.data, -1)
  }
  
  ################### here, we treat censoring as the event #################################################
  status.c <- c(rep(0, nrow(time.data)-1), 1-cens) 
  status.t <- c(rep(0, nrow(time.data)-1), cens) 
  status.u <- c(rep(0, nrow(time.data)-1), 1) 
  
  time.data <- cbind(time.data, status.c, status.t, status.u)
  
  stack.data <- rbind(stack.data, time.data)
  long.data <- rbind(long.data,  time.data[nrow(time.data),])
}
long.data <- cbind(long.data, outcome = new.Y)
stack.data <- as.data.frame(stack.data) ############# the data set used for cox model
long.data <- as.data.frame(long.data)
base.time.cov <- stack.data$time.cov[stack.data$tstart==0]
base.data <- long.data
base.data$time.cov <- base.time.cov

base.data$w <- 1
long.data$w <- 1
stack.data$w <- 1

base.data$id <- c(1:N)
long.data$id <- c(1:N)
stack.data$id <- rep(1:N, times = table(stack.data$id))

y <- long.data$outcome
X <- long.data[, 4:10]
time <- long.data$tstop


#==============================================================#
#=============== Step 2: laod packages ========================#
#==============================================================#
########## time dependent cox model for censoring process
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
#=============== Step 3: fit different outcome models,=========#
#============== which will be inserted into the estimation ====#
#==============================================================#
K <- 150
new.L <- c(0,6,12,18,24)
stop.L <- c(6,12,18,24,30)
test.time <- c(16,18,20,22,24)

NPB.td.gam <- NPB.td.correct <- NPB.td.SL <- array(numeric(0), dim = c(length(test.time),3,K))
NPB.ti.gam <- NPB.ti.correct <- NPB.ti.SL <- array(numeric(0), dim = c(length(test.time),3,K))
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

################################  the outcome regression #################################################
### semiparametric
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
#=============== Step 4: the proposed method ==================#
#==============================================================#
# the proposed estimator
mr.mod.td <- mr.fun.disc(test.time, long.data, stack.data, type = c("lm", "gam", "SL"))

# then we do the bootstrap-based variance estimation
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
#============== Step 5: compare with Lucy's method ============#
#==============================================================#
# the lucy's estimator

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

# the bootstrap based variance estimation
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
#========= Step 6: summarize and save the results =============#
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
