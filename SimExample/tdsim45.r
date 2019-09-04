###########created by Hao Sun at 4/20/2018

#########################################################################################################################
######################### simulation study with time dependent confounders ##############################################
#########################################################################################################################


#############################################  simulate the data based on the real data #################################
############### the function when cosnidering time dependent confounders
source("/gpfs/fs1/home/hsun16/my simulation/Treatment length/sim2018/tdcase/DatGen.r")
source("/gpfs/fs1/home/hsun16/my simulation/Treatment length/sim2018/tdcase/TrtLenFun.r")
source("/gpfs/fs1/home/hsun16/my simulation/Treatment length/sim2018/tdcase/MeanFuncDisc.r")
############### the function when cosnidering only baseline covariates
source("/gpfs/fs1/home/hsun16/my simulation/Treatment length/sim2018/ticase/MeanFuncInd.r")

library(survival)
library(gam)
library(np)
##############################################################################################################################
################ simulation set ups #########################################################################################
############################################################################################################################## 
test.time <- c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29)
new.L <- c(0,6,12,18,24)
stop.L <- c(6,12,18,24,30)
K <- 150 ### the resamping numbers
n <- 500 ### the sample size

#######################################################################################################################
############################## the baseline hazard function for td case with 45% censoring ############################
############################## hard to find it ... ####################################################################
############################## old one ################################################################################
#baset1 <- 0.014; baset2 <- 0.015; baset3 <- 0.02; baset4 <- 0.035; baset5 <- 0.12
#basec1 <- 0.01; basec2 <- 0.01; basec3 <- 0.015; basec4 <- 0.028; basec5 <- 0.08
############################## new one ################################################################################
baset1 <- 0.0155; baset2 <- 0.02; baset3 <- 0.028; baset4 <- 0.044; baset5 <- 0.14
basec1 <- 0.0135; basec2 <- 0.016; basec3 <- 0.024; basec4 <- 0.038; basec5 <- 0.11


########################################################################################################################
############################### the baseline hazard function for td case with 20% censoring ############################
############################### hard to find it ... ####################################################################
############################### old one ################################################################################
#baset1 <- 0.02; baset2 <- 0.028; baset3 <- 0.037; baset4 <- 0.06; baset5 <- 0.23
#basec1 <- 0.0035; basec2 <- 0.0035; basec3 <- 0.0035; basec4 <- 0.004; basec5 <- 0.016
############################### new one ################################################################################
#baset1 <- 0.022; baset2 <- 0.03; baset3 <- 0.040; baset4 <- 0.062; baset5 <- 0.23
#basec1 <- 0.006; basec2 <- 0.0075; basec3 <- 0.01; basec4 <- 0.01; basec5 <- 0.04


################################ output folder ################################################################################
today <- paste0("td",n,"ae45wellnewbinary")
#today <- paste0("td",n,"ae45wellnewgaussian")
output_folder <- paste0("output_",today)
if (!file.exists(output_folder)){dir.create(output_folder)}

## SLURM_ARRAY_TASK_ID on bluehive
SLURM_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(100+SLURM_ID)


library(survival)
library(gam)
library(randomForest)
require("SuperLearner")
require("earth")
require("gam")
require("ranger")
require(KernSmooth)

NPB.td.gam <- NPB.td.plm <- NPB.td.SL <- NPB.td.correct <- array(numeric(0), dim = c(length(test.time),3,K))
NPB.td.tr.gam <- NPB.td.tr.plm <- NPB.td.tr.SL <- NPB.td.tr.correct <- array(numeric(0), dim = c(length(test.time),3,K))
NPB.ti.gam <- NPB.ti.SL <- NPB.ti.correct <- array(numeric(0), dim = c(length(test.time),3,K))
NPB.ti.tr.gam <- NPB.ti.tr.correct <- NPB.ti.tr.SL  <- array(numeric(0), dim = c(length(test.time),3,K))

################################################################################
########################## generate the data ###################################
################################################################################
  dat.mod <- dat.gen(n, base.cen, sp, type = "binary")
  
  ########## data structure needed for proposed method
  full.data <- dat.mod$full.data; full.data$w <- 1
  long.data <- dat.mod$long.data; long.data$w <- 1
  stack.data <- dat.mod$stack.data; stack.data$w <- 1
  
  ########## data structure needed for lucy's method
  base.data <- dat.mod$base.data;base.data$w <- 1
  
  ########### generate the true weight ########################################
  weightmat <- true.weight(new.L, full.data)

  #############################################################################
  ########## the estimation by using proposed method ##########################
  #############################################################################
  ########## fit the support functions 
  c.cox <- coxph(Surv(tstart, tstop, status.c) ~ time.cov+dbts,data=stack.data)
  gamma.hat <- coefficients(c.cox)
  ### the cumulative baseline hazard function
  base.c <- basehaz(c.cox, centered = F)
  
  t.cox <- coxph(Surv(tstart, tstop, status.t) ~ time.cov+dbts,data=stack.data)
  eta.hat <- coefficients(t.cox)
  ### the cumulative baseline hazard function
  base.t <- basehaz(t.cox, centered = F)
  
  
  ################################  the outcome regression #####################
  ### spline regression
  gam.linear.t <- gam(outcome ~ s(tstop) + time.cov + dbts,
                     data = long.data[long.data$status.c == 0, ], family=binomial(link='logit'))
  gam.linear.c <- gam(outcome ~ s(tstop) +  time.cov + dbts,
                      data = long.data[long.data$status.c == 1,], family=binomial(link='logit'))
  #gam.linear.t <- gam(outcome ~ s(tstop) + time.cov + dbts,
  #                   data = long.data[long.data$status.c == 0, ])
  #gam.linear.c <- gam(outcome ~ s(tstop) +  time.cov + dbts,
  #                    data = long.data[long.data$status.c == 1,])


  ### partial linear model
  #bw.t <- npplregbw(formula= outcome ~ time.cov + dbts|tstop, data = long.data[long.data$status.c == 0, ], regtype = "ll")
  #plm.t <- npplreg(bws=bw.t)
  #bw.c <- npplregbw(formula= outcome ~ time.cov + dbts|tstop, data = long.data[long.data$status.c == 1, ], regtype = "ll")
  #plm.c <- npplreg(bws=bw.c)

  ### nonparametric super learner
  #sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glm.interaction","SL.mean","SL.ranger","SL.randomForest")
  #xt <- long.data[long.data$status.c == 0, 3:5]
  #xc <- long.data[long.data$status.c == 1, 3:5]
  #outcomet <- long.data[long.data$status.c == 0,]$outcome
  #outcomec <- long.data[long.data$status.c == 1,]$outcome
  #SL.t <- SuperLearner(Y=outcomet, X=xt, SL.library=sl.lib)
  #SL.c <- SuperLearner(Y=outcomec, X=xc, SL.library=sl.lib)

  ### randomforest
  #rf.t <- randomForest(xt,outcomet)
  #rf.c <- randomForest(xc,outcomec)  

  ### parametric
  lm.t <- glm(outcome ~ sqrt(tstop) + time.cov + dbts,
             data = long.data[long.data$status.c == 0, ], family=binomial(link='logit'))
  lm.c <- glm(outcome ~ log(tstop) +  time.cov + dbts,
             data = long.data[long.data$status.c == 1,], family=binomial(link='logit'))
  #lm.t <- glm(outcome ~ sqrt(tstop) + time.cov + dbts,
  #           data = long.data[long.data$status.c == 0, ])
  #lm.c <- glm(outcome ~ log(tstop) +  time.cov + dbts,
  #           data = long.data[long.data$status.c == 1,])

  N <- nrow(long.data)
  
  #weightmat <- true.weight(new.L, full.data)
  ######### estimate the mean
  mr.mod.td <- mr.fun.disc(test.time, long.data, stack.data, type = c("lm", "gam"))
  
  #################################################################################
  ########################### variance estimation #################################
  #################################################################################
 for(k in 1:K){
        set.seed(k+SLURM_ID)
	print(k)
	 
	idx <- sample(n, replace = T)
	re.long <- long.data[idx,]
	re.long$id <- 1:n
	re.full <- lapply(1:n, function(x){a <- full.data[full.data$id == idx[x],]; a$id <- x; a})
	re.full <- Reduce(rbind, re.full)
	re.stack <- lapply(1:n, function(x){a <- stack.data[stack.data$id == idx[x],]; a$id <- x; a})
	re.stack <- Reduce(rbind, re.stack)
	re.base <- base.data[idx,]
	#re.weight <- weightmat[idx,]


	#################################################################################
        ##################################  NPB part  ###################################
 	#################################################################################
	########## permute the data and refit all support functions
	########## fit the support functions 
  	c.cox <- coxph(Surv(tstart, tstop, status.c) ~ time.cov+dbts,data=re.stack)
  	gamma.hat <- coefficients(c.cox)
  	### the cumulative baseline hazard function
  	base.c <- basehaz(c.cox, centered = F)
  
     	t.cox <- coxph(Surv(tstart, tstop, status.t) ~ time.cov+dbts,data=re.stack)
  	eta.hat <- coefficients(t.cox)
  	### the cumulative baseline hazard function
  	base.t <- basehaz(t.cox, centered = F)
  
  
  	################################  the outcome regression #################################################
  	### semiparametric
  	gam.linear.t <- gam(outcome ~ s(tstop) + time.cov + dbts,
                      	data = re.long[re.long$status.c == 0, ], family=binomial(link='logit'))
  	gam.linear.c <- gam(outcome ~ s(tstop) +  time.cov + dbts,
                    	 data = re.long[re.long$status.c == 1,], family=binomial(link='logit'))
        #gam.linear.t <- gam(outcome ~ s(tstop) + time.cov + dbts,
        #              	data = re.long[re.long$status.c == 0, ])
  	#gam.linear.c <- gam(outcome ~ s(tstop) +  time.cov + dbts,
        #            	 data = re.long[re.long$status.c == 1,])


        ##bw.t <- npplregbw(formula= outcome ~ time.cov + dbts|tstop, data = re.long[re.long$status.c == 0, ], regtype = "ll")
        #plm.t <- npplreg(bws=as.matrix(unlist(bw.t$bandwidth)), formula= outcome ~ time.cov + dbts|tstop, data = re.long[re.long$status.c == 0, ], regtype = "ll")
        ##bw.c <- npplregbw(formula= outcome ~ time.cov + dbts|tstop, data = re.long[re.long$status.c == 1, ], regtype = "ll")
        #plm.c <- npplreg(bws=as.matrix(unlist(bw.c$bandwidth)), formula= outcome ~ time.cov + dbts|tstop, data = re.long[re.long$status.c == 1, ], regtype = "ll")

  	### nonparametric
  	#sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glm.interaction","SL.mean","SL.ranger")
  	#xt <- re.long[re.long$status.c == 0, 3:5]
  	#xc <- re.long[re.long$status.c == 1, 3:5]
  	#outcomet <- re.long[re.long$status.c == 0,]$outcome
  	#outcomec <- re.long[re.long$status.c == 1,]$outcome
  	#SL.t <- SuperLearner(Y=outcomet, X=xt, SL.library=sl.lib)
  	#SL.c <- SuperLearner(Y=outcomec, X=xc, SL.library=sl.lib)
  	### parametric
  	lm.t <- glm(outcome ~ sqrt(tstop) + time.cov + dbts,
             data = re.long[re.long$status.c == 0, ], family=binomial(link='logit'))
  	lm.c <- glm(outcome ~ log(tstop) +  time.cov + dbts,
             data = re.long[re.long$status.c == 1,], family=binomial(link='logit'))
        #lm.t <- glm(outcome ~ sqrt(tstop) + time.cov + dbts,
        #     data = re.long[re.long$status.c == 0, ])
  	#lm.c <- glm(outcome ~ log(tstop) +  time.cov + dbts,
        #     data = re.long[re.long$status.c == 1,])

  	N <- nrow(re.long)
  
  	######### estimate the mean
  	re.td <- mr.fun.disc(test.time, re.long, re.stack, type = c("lm", "gam"))

	######### save the result to calculate the resampling variance estimation
	NPB.td.gam[,,k] <- re.td$gam[,-1]
	NPB.td.plm[,,k] <- re.td$plm[,-1]
	NPB.td.correct[,,k] <- re.td$correct[,-1]
	#NPB.td.tr.gam[,,k] <- re.td$gam.tr[,-1]
	#NPB.td.tr.plm[,,k] <- re.td$plm.tr[,-1]
	#NPB.td.tr.correct[,,k] <- re.td$correct.tr[,-1]

	

	#################################################################################
        ################################## MB   part  ###################################
 	#################################################################################
	#set.seed(k+SLURM_ID)
	#w <- rexp(n, rate = 1)
	#re.long <- long.data
	#re.long$w <- w
	#re.stack <- stack.data
	#re.stack$w <- rep(w, table(re.stack$id))

	#c.cox <- coxph(Surv(tstart, tstop, status.c) ~ time.cov+dbts, weight = w, data=re.stack)
	#gamma.hat <- coefficients(c.cox)
	#base.c <- basehaz(c.cox, centered = F)

	########## time dependent cox model for U = min(T, C)
	#t.cox <- coxph(Surv(tstart, tstop, status.t) ~ time.cov+dbts,data=re.stack, weight = w)
	#eta.hat <- coefficients(t.cox)
	#base.t <- basehaz(t.cox, centered = F)


	#################################  the outcome regression #################################################
	#gam.linear.t <- gam(outcome ~ s(tstop) + time.cov + dbts,
	#                    data = re.long[re.long$status.c == 0, ], weight = w)
	#gam.linear.c <- gam(outcome ~ s(tstop) +  time.cov + dbts,
     	#            data = re.long[re.long$status.c == 1,], weight = w)

	#lm.t <- lm(outcome ~ sqrt(tstop) + time.cov + dbts,
        #   data = re.long[re.long$status.c == 0, ], weights = w)
	#lm.c <- lm(outcome ~ log(tstop) +  time.cov + dbts,
        #   data = re.long[re.long$status.c == 1,], weights = w)
	#N <- nrow(re.long)
  
  	######### estimate the mean
  	#re.td <- mr.fun.disc(test.time, re.long, re.stack)
	######### save the result to calculate the resampling variance estimation
	#MB.td.sepa[,,k] <- re.td$sepa[,-1]
	#MB.td.correct[,,k] <- re.td$correct[,-1]
 }

  #############################################################################
  ############################# end ################ ##########################
  #############################################################################
  
  
  #############################################################################
  ########## the estimation by using lucy's method ############################
  #############################################################################
  c.cox <- coxph(Surv(tstop, status.c) ~ time.cov+dbts,data= base.data)
  gamma.hat <- coefficients(c.cox)
  ### the cumulative baseline hazard function
  base.c <- basehaz(c.cox, centered = F)
  
  #################################  the outcome regression ###################
  gam.linear.t <- gam(outcome ~ s(tstop) + time.cov + dbts,
                      data = base.data[base.data$status.c == 0, ], family=binomial(link='logit'))
  gam.linear.c <- gam(outcome ~ s(tstop) +  time.cov + dbts,
                      data = base.data[base.data$status.c == 1,], family=binomial(link='logit'))
  #gam.linear.t <- gam(outcome ~ s(tstop) + time.cov + dbts,
  #                    data = base.data[base.data$status.c == 0, ])
  #gam.linear.c <- gam(outcome ~ s(tstop) +  time.cov + dbts,
  #                    data = base.data[base.data$status.c == 1,])

  
  #xt <- base.data[base.data$status.c == 0, 3:5]
  #xc <- base.data[base.data$status.c == 1, 3:5]
  #outcomet <- base.data[base.data$status.c == 0,]$outcome
  #outcomec <- base.data[base.data$status.c == 1,]$outcome
  #SL.t <- SuperLearner(Y=outcomet, X=xt, SL.library=sl.lib)
  #SL.c <- SuperLearner(Y=outcomec, X=xc, SL.library=sl.lib)
  
  ####### parametric model 
  lm.t <- glm(outcome ~ sqrt(tstop) + time.cov + dbts,
             data = base.data[base.data$status.c == 0, ], family=binomial(link='logit'))
  lm.c <- glm(outcome ~ log(tstop) +  time.cov + dbts,
             data = base.data[base.data$status.c == 1,], family=binomial(link='logit'))
  #lm.t <- glm(outcome ~ sqrt(tstop) + time.cov + dbts,
  #           data = base.data[base.data$status.c == 0, ])
  #lm.c <- glm(outcome ~ log(tstop) +  time.cov + dbts,
  #           data = base.data[base.data$status.c == 1,])

  N <- nrow(base.data)  

  ######### estimate the mean
  mr.mod.ti <- mr.fun.ind(test.time, base.data)
  
  #################################################################################
  ########################### variance estimation #################################
  #################################################################################
 

 for(k in 1:K){
   set.seed(k+SLURM_ID)
   print(k)
   
   idx <- sample(n, replace = T)
   re.base <- base.data[idx,]
   
   #################################################################################
   ##################################  NPB part  ###################################
   #################################################################################
   ########## permute the data and refit all support functions
   ########## fit the support functions 
   #############################################################################
   c.cox <- coxph(Surv(tstop, status.c) ~ time.cov+dbts,data= re.base)
   gamma.hat <- coefficients(c.cox)
   ### the cumulative baseline hazard function
   base.c <- basehaz(c.cox, centered = F)
   
   #################################  the outcome regression ###################
   gam.linear.t <- gam(outcome ~ s(tstop) + time.cov + dbts,
                       data = re.base[re.base$status.c == 0, ], family=binomial(link='logit'))
   gam.linear.c <- gam(outcome ~ s(tstop) +  time.cov + dbts,
                       data = re.base[re.base$status.c == 1,], family=binomial(link='logit'))
   #gam.linear.t <- gam(outcome ~ s(tstop) + time.cov + dbts,
   #                    data = re.base[re.base$status.c == 0, ])
   #gam.linear.c <- gam(outcome ~ s(tstop) +  time.cov + dbts,
   #                   data = re.base[re.base$status.c == 1,])

   
   #xt <- re.base[re.base$status.c == 0, 3:5]
   #xc <- re.base[re.base$status.c == 1, 3:5]
   #outcomet <- re.base[re.base$status.c == 0,]$outcome
   #outcomec <- re.base[re.base$status.c == 1,]$outcome
   #SL.t <- SuperLearner(Y=outcomet, X=xt, SL.library=sl.lib)
   #SL.c <- SuperLearner(Y=outcomec, X=xc, SL.library=sl.lib)
   lm.t <- glm(outcome ~ sqrt(tstop) + time.cov + dbts,
              data = re.base[re.base$status.c == 0, ], family=binomial(link='logit'))
   lm.c <- glm(outcome ~ log(tstop) +  time.cov + dbts,
              data = re.base[re.base$status.c == 1,], family=binomial(link='logit'))
   #lm.t <- glm(outcome ~ sqrt(tstop) + time.cov + dbts,
   #           data = re.base[re.base$status.c == 0, ])
   #lm.c <- glm(outcome ~ log(tstop) +  time.cov + dbts,
   #           data = re.base[re.base$status.c == 1,])

   N <- nrow(re.base)  
   
   ######### estimate the mean
   re.ti <- mr.fun.ind(test.time, re.base)
   
   ######### save the result to calculate the resampling variance estimation
   NPB.ti.gam[,,k] <- re.ti$gam[,-1]
   NPB.ti.correct[,,k] <- re.ti$correct[,-1]
   #NPB.ti.SL[,,k] <- re.ti$SL[,-1]
   
   
   
   #################################################################################
#   ##################################  MB part  ###################################
#   #################################################################################
#   #set.seed(k+SLURM_ID)
#   #w <- rexp(n, rate = 1)
#   #re.base <- long.data
#   #re.base$w <- w
#   
#   #c.cox <- coxph(Surv(tstop, status.c) ~ time.cov+dbts, weight = w, data=re.base)
#   #gamma.hat <- coefficients(c.cox)
#   #base.c <- basehaz(c.cox, centered = F)
#   
#   #################################  the outcome regression #################################################
#   #gam.linear.t <- gam(outcome ~ s(tstop) + time.cov + dbts,
#   #                    data = re.base[re.base$status.c == 0, ], weight = w)
#   #gam.linear.c <- gam(outcome ~ s(tstop) +  time.cov + dbts,
#   #	            data = re.base[re.base$status.c == 1,], weight = w)
#   
#   #lm.t <- lm(outcome ~ sqrt(tstop) + time.cov + dbts,
#   #   data = re.base[re.base$status.c == 0, ], weights = w)
#   #lm.c <- lm(outcome ~ log(tstop) +  time.cov + dbts,
#   #   data = re.base[re.base$status.c == 1,], weights = w)
#   #N <- nrow(re.base)
#   
#   ######### estimate the mean
#   #re.ti <- mr.fun.ind(test.time, re.base)
#   ######### save the result to calculate the resampling variance estimation
#   #MB.ti.sepa[,,k] <- re.ti$sepa[,-1]
#   #MB.ti.correct[,,k] <- re.ti$correct[,-1]
}

  #############################################################################
  ############################# end ################ ##########################
  #############################################################################
  
  ########## record the result ################################################
  mean.gam.td <- mr.mod.td$gam[,-1]
  mean.plm.td <- mr.mod.td$plm[,-1]
  mean.correct.td <- mr.mod.td$correct[,-1]
  mean.SL.td <- mr.mod.td$SL[,-1]
  p1.weight <- mr.mod.td$p1.weight

  #mean.gam.td.tr <- mr.mod.td$gam.tr[,-1]
  #mean.plm.td.tr <- mr.mod.td$plm.tr[,-1]
  #mean.correct.td.tr <- mr.mod.td$correct.tr[,-1]
  #mean.SL.td.tr <- mr.mod.td$SL.tr[,-1]
  #p1.weight.tr <- mr.mod.td$p1.weight.tr


  mean.gam.ti <- mr.mod.ti$gam[,-1]
  mean.correct.ti <- mr.mod.ti$correct[,-1]
  #mean.SL.ti <- mr.mod.ti$SL[,-1]



result.td <- list("gam.est" = mean.gam.td, "plm.est" = mean.plm.td, "correct.est" = mean.correct.td, "SL.est" = mean.SL.td,"p1.weight" = p1.weight,
	 "NPB.gam" = NPB.td.gam, "NPB.plm" = NPB.td.plm, "NPB.correct" = NPB.td.correct, "NPB.SL" = NPB.td.SL)
result.ti <- list("gam.est" = mean.gam.ti, "correct.est" = mean.correct.ti,"NPB.gam" = NPB.ti.gam, "NPB.correct" = NPB.ti.correct)
result = list(result.td = result.td, result.ti = result.ti)
save(result, file = paste0(output_folder,"/mean",SLURM_ID,".rda"))
