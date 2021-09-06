library(depmixS4)

nsim <- 1000
nrep <- 100
nt <- 50

# ## simulation 1
set.seed(1234)
randomSeeds <- sample(seq(1,nsim*1000),nsim)
# 
prior <- c(8,1,1)
prior <- prior/sum(prior)
transition <- 5*diag(3) + 1
transition <- transition/rowSums(transition)
means <- c(-1,0,1)
sds <- c(1,1,1)
pmiss <- c(.05,.25,.5)

truepars1 <- c(prior,as.numeric(t(transition)),as.numeric(rbind(means,sds)))
truepars2 <- c(prior,as.numeric(t(transition)),as.numeric(rbind(means,sds,1-pmiss,pmiss)))

out1 <- rep(0.0, length=nsim)
out2 <- rep(0.0, length=nsim)
for(sim in 1:nsim) {
  set.seed(randomSeeds[sim])
  truestate <- matrix(nrow=nt,ncol=nrep)
  for(i in 1:nrep) {
    truestate[1,i] <- sample(1:3,size=1,prob=prior)
    for(t in 2:nt) {
      truestate[t,i] <- sample(1:3,size=1,prob=transition[truestate[t-1,i],])
    }
  }
  dat <- data.frame(trueState=as.numeric(truestate),trial=1:nt)
  dat$trueResponse <- rnorm(nrow(dat),mean=means[dat$trueState],sd=sds[dat$trueState])
  dat$missing <- rbinom(nrow(dat),size=1,prob=pmiss[dat$trueState])
  dat$response <- dat$trueResponse
  dat$response[dat$missing==1] <- NA
  
  mod <- depmix(list(response~1),family=list(gaussian()),data=dat,nstates=3,ntimes=rep(nt,nrep))
  mod <- setpars(mod,truepars1)
  out1[sim] <- sum(posterior(mod, type="global") == dat$trueState)/nrow(dat)
  
  mod <- depmix(list(response~1,missing~1),family=list(gaussian(),multinomial("identity")),data=dat,nstates=3,ntimes=rep(nt,nrep))
  mod <- setpars(mod,truepars2)
  out2[sim] <- sum(posterior(mod, type="global") == dat$trueState)/nrow(dat)
}
sim1_MAR_acc <- out1
sim1_MNAR_acc <- out2

# ## simulation 2
set.seed(1234)
randomSeeds <- sample(seq(1,nsim*1000),nsim)
# 
prior <- c(8,1,1)
prior <- prior/sum(prior)
transition <- 5*diag(3) + 1
transition <- transition/rowSums(transition)
means <- c(-1,0,1)
sds <- c(1,1,1)
pmiss <- c(.25,.25,.25)

truepars1 <- c(prior,as.numeric(t(transition)),as.numeric(rbind(means,sds)))
truepars2 <- c(prior,as.numeric(t(transition)),as.numeric(rbind(means,sds,1-pmiss,pmiss)))

out1 <- rep(0.0, length=nsim)
out2 <- rep(0.0, length=nsim)
for(sim in 1:nsim) {
  set.seed(randomSeeds[sim])
  truestate <- matrix(nrow=nt,ncol=nrep)
  for(i in 1:nrep) {
    truestate[1,i] <- sample(1:3,size=1,prob=prior)
    for(t in 2:nt) {
      truestate[t,i] <- sample(1:3,size=1,prob=transition[truestate[t-1,i],])
    }
  }
  dat <- data.frame(trueState=as.numeric(truestate),trial=1:nt)
  dat$trueResponse <- rnorm(nrow(dat),mean=means[dat$trueState],sd=sds[dat$trueState])
  dat$missing <- rbinom(nrow(dat),size=1,prob=pmiss[dat$trueState])
  dat$response <- dat$trueResponse
  dat$response[dat$missing==1] <- NA
  
  mod <- depmix(list(response~1),family=list(gaussian()),data=dat,nstates=3,ntimes=rep(nt,nrep))
  mod <- setpars(mod,truepars1)
  out1[sim] <- sum(posterior(mod, type="global") == dat$trueState)/nrow(dat)
  
  mod <- depmix(list(response~1,missing~1),family=list(gaussian(),multinomial("identity")),data=dat,nstates=3,ntimes=rep(nt,nrep))
  mod <- setpars(mod,truepars2)
  out2[sim] <- sum(posterior(mod, type="global") == dat$trueState)/nrow(dat)
}
sim2_MAR_acc <- out1
sim2_MNAR_acc <- out2

# ## simulation 3
set.seed(1234)
randomSeeds <- sample(seq(1,nsim*1000),nsim)
# 
prior <- c(8,1,1)
prior <- prior/sum(prior)
transition <- 5*diag(3) + 1
transition <- transition/rowSums(transition)
means <- c(-1,0,1)
sds <- c(3,3,3)
pmiss <- c(.05,.25,.5)

truepars1 <- c(prior,as.numeric(t(transition)),as.numeric(rbind(means,sds)))
truepars2 <- c(prior,as.numeric(t(transition)),as.numeric(rbind(means,sds,1-pmiss,pmiss)))

out1 <- rep(0.0, length=nsim)
out2 <- rep(0.0, length=nsim)
for(sim in 1:nsim) {
  set.seed(randomSeeds[sim])
  truestate <- matrix(nrow=nt,ncol=nrep)
  for(i in 1:nrep) {
    truestate[1,i] <- sample(1:3,size=1,prob=prior)
    for(t in 2:nt) {
      truestate[t,i] <- sample(1:3,size=1,prob=transition[truestate[t-1,i],])
    }
  }
  dat <- data.frame(trueState=as.numeric(truestate),trial=1:nt)
  dat$trueResponse <- rnorm(nrow(dat),mean=means[dat$trueState],sd=sds[dat$trueState])
  dat$missing <- rbinom(nrow(dat),size=1,prob=pmiss[dat$trueState])
  dat$response <- dat$trueResponse
  dat$response[dat$missing==1] <- NA
  
  mod <- depmix(list(response~1),family=list(gaussian()),data=dat,nstates=3,ntimes=rep(nt,nrep))
  mod <- setpars(mod,truepars1)
  out1[sim] <- sum(posterior(mod, type="global") == dat$trueState)/nrow(dat)
  
  mod <- depmix(list(response~1,missing~1),family=list(gaussian(),multinomial("identity")),data=dat,nstates=3,ntimes=rep(nt,nrep))
  mod <- setpars(mod,truepars2)
  out2[sim] <- sum(posterior(mod, type="global") == dat$trueState)/nrow(dat)
}
sim3_MAR_acc <- out1
sim3_MNAR_acc <- out2

# ## simulation 4
set.seed(1234)
randomSeeds <- sample(seq(1,nsim*1000),nsim)
# 
prior <- c(8,1,1)
prior <- prior/sum(prior)
transition <- 5*diag(3) + 1
transition <- transition/rowSums(transition)
means <- c(-1,0,1)
sds <- c(3,3,3)
pmiss <- c(.25,.25,.25)

truepars1 <- c(prior,as.numeric(t(transition)),as.numeric(rbind(means,sds)))
truepars2 <- c(prior,as.numeric(t(transition)),as.numeric(rbind(means,sds,1-pmiss,pmiss)))

out1 <- rep(0.0, length=nsim)
out2 <- rep(0.0, length=nsim)
for(sim in 1:nsim) {
  set.seed(randomSeeds[sim])
  truestate <- matrix(nrow=nt,ncol=nrep)
  for(i in 1:nrep) {
    truestate[1,i] <- sample(1:3,size=1,prob=prior)
    for(t in 2:nt) {
      truestate[t,i] <- sample(1:3,size=1,prob=transition[truestate[t-1,i],])
    }
  }
  dat <- data.frame(trueState=as.numeric(truestate),trial=1:nt)
  dat$trueResponse <- rnorm(nrow(dat),mean=means[dat$trueState],sd=sds[dat$trueState])
  dat$missing <- rbinom(nrow(dat),size=1,prob=pmiss[dat$trueState])
  dat$response <- dat$trueResponse
  dat$response[dat$missing==1] <- NA
  
  mod <- depmix(list(response~1),family=list(gaussian()),data=dat,nstates=3,ntimes=rep(nt,nrep))
  mod <- setpars(mod,truepars1)
  out1[sim] <- sum(posterior(mod, type="global") == dat$trueState)/nrow(dat)
  
  mod <- depmix(list(response~1,missing~1),family=list(gaussian(),multinomial("identity")),data=dat,nstates=3,ntimes=rep(nt,nrep))
  mod <- setpars(mod,truepars2)
  out2[sim] <- sum(posterior(mod, type="global") == dat$trueState)/nrow(dat)
}
sim4_MAR_acc <- out1
sim4_MNAR_acc <- out2


## simulation 5
set.seed(1234)
randomSeeds <- sample(seq(1,nsim*1000),nsim)
out <- rep(list(vector("list",3)),nsim)

prior <- c(8,1,1)
prior <- prior/sum(prior)
transition <- 5*diag(3) + 1
transition <- transition/rowSums(transition)
means <- c(-1,0,1)
sds <- c(1,1,1)
pmiss <- rep(mean(1/(1+exp(-(-5+.125*1:50)))),3)

truepars1 <- c(prior,as.numeric(t(transition)),as.numeric(rbind(means,sds)))
truepars2 <- c(prior,as.numeric(t(transition)),as.numeric(rbind(means,sds,1-pmiss,pmiss)))
truepars3 <- c(prior,as.numeric(t(transition)),as.numeric(rbind(means,sds,c(-5,-5,-5),c(.125,.125,.125))))

out1 <- rep(0.0, length=nsim)
out2 <- rep(0.0, length=nsim)
out3 <- rep(0.0, length=nsim)
for(sim in 1:nsim) {
  set.seed(randomSeeds[sim])
  truestate <- matrix(nrow=nt,ncol=nrep) 
  for(i in 1:nrep) {
    truestate[1,i] <- sample(1:3,size=1,prob=prior)
    for(t in 2:nt) {
      truestate[t,i] <- sample(1:3,size=1,prob=transition[truestate[t-1,i],])
    }
  }
  dat <- data.frame(trueState=as.numeric(truestate),trial=1:nt)
  dat$trueResponse <- rnorm(nrow(dat),mean=means[dat$trueState],sd=sds[dat$trueState])
  dat$missing <- rbinom(nrow(dat),size=1,prob=1/(1+exp(-(-5 + .125*1:50))))
  dat$response <- dat$trueResponse
  dat$response[dat$missing==1] <- NA
  
  set.seed(randomSeeds[sim])
  mod <- depmix(list(response~1),family=list(gaussian()),data=dat,nstates=3,ntimes=rep(nt,nrep))
  mod <- setpars(mod,truepars1)
  out1[sim] <- sum(posterior(mod, type="global") == dat$trueState)/nrow(dat)
  
  mod <- depmix(list(response~1,missing~1),family=list(gaussian(),multinomial("identity")),data=dat,nstates=3,ntimes=rep(nt,nrep))
  mod <- setpars(mod,truepars2)
  out2[sim] <- sum(posterior(mod, type="global") == dat$trueState)/nrow(dat)
  
  mod <- depmix(list(response~1,missing~trial),family=list(gaussian(),binomial()),data=dat,nstates=3,ntimes=rep(nt,nrep))
  mod <- setpars(mod,truepars3)
  out3[sim] <- sum(posterior(mod, type="global") == dat$trueState)/nrow(dat)
  
}
sim5_MAR_acc <- out1
sim5_MNAR_acc <- out2
sim5_MNAR_time_acc <- out3

save(sim1_MAR_acc, sim1_MNAR_acc, sim2_MAR_acc, sim2_MNAR_acc, sim3_MAR_acc, sim3_MNAR_acc, sim4_MAR_acc, sim4_MNAR_acc, sim5_MAR_acc, sim5_MNAR_acc, sim5_MAR_time_acc, file="classification_accuracy.RData")