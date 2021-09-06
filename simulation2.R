library(depmixS4)

nsim <- 1000
nrep <- 100
nt <- 50

## simulation 2
set.seed(202109012)
randomSeeds <- sample(seq(1,nsim*1000),nsim)
out <- rep(list(vector("list",3)),nsim)

prior <- c(8,1,1)
prior <- prior/sum(prior)
transition <- 5*diag(3) + 1
transition <- transition/rowSums(transition)
means <- c(-1,0,1)
sds <- c(1,1,1)
pmiss <- c(.25,.25,.25)

truepars1 <- c(prior,as.numeric(t(transition)),as.numeric(rbind(means,sds)))
truepars2 <- c(prior,as.numeric(t(transition)),as.numeric(rbind(means,sds,1-pmiss,pmiss)))
truepars3 <- c(prior,as.numeric(t(transition)),as.numeric(rbind(means,sds,log(pmiss/(1-pmiss)),c(0,0,0))))

generate_data <- function() {
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
  return(dat)
}

for(sim in 1:nsim) {
  set.seed(randomSeeds[sim])
  dat <- generate_data()
  
  set.seed(randomSeeds[sim])
  mod <- depmix(list(response~1),family=list(gaussian()),data=dat,nstates=3,ntimes=rep(nt,nrep))
  mod <- setpars(mod,truepars1)
  ok <- FALSE
  nTry <- 1
  while(!ok && nTry < 21) {
    fmod <- try(fit(mod,emcontrol=em.control(maxit = 4000, random.start = nTry > 1),verbose=FALSE))
    if(!inherits(fmod,"try-error")) {
      if(fmod@message == "Log likelihood converged to within tol. (relative change)") {
        cat("mod1 simulation",sim,"iteration",nTry,"completed.\n")
        ok <- TRUE
        out[[sim]][[1]] <- list(pars=getpars(fmod),logLik=logLik(fmod),viterbi=posterior(fmod, type="global"),trueState=dat$trueState)
      } else {
        cat("mod1 simulation",sim,"iteration",nTry,"failed.\n")
        nTry <- nTry + 1
      }
    } else {
      cat("mod1 simulation",sim,"iteration",nTry,"failed.\n")
      nTry <- nTry + 1
    }
  }
  
  
  set.seed(randomSeeds[sim])
  mod <- depmix(list(response~1,missing~1),family=list(gaussian(),multinomial("identity")),data=dat,nstates=3,ntimes=rep(nt,nrep))
  mod <- setpars(mod,truepars2)
  ok <- FALSE
  nTry <- 1
  while(!ok && nTry < 21) {
    fmod <- try(fit(mod,emcontrol=em.control(maxit = 4000, random.start = nTry > 1),verbose=FALSE))
    if(!inherits(fmod,"try-error")) {
      if(fmod@message == "Log likelihood converged to within tol. (relative change)") {
        cat("mod2 simulation",sim,"iteration",nTry,"completed.\n")
        ok <- TRUE
        out[[sim]][[2]] <- list(pars=getpars(fmod),logLik=logLik(fmod),viterbi=posterior(fmod, type="global"),trueState=dat$trueState)
      } else {
        cat("mod2 simulation",sim,"iteration",nTry,"failed.\n")
        nTry <- nTry + 1
      }
    } else {
      cat("mod2 simulation",sim,"iteration",nTry,"failed.\n")
      nTry <- nTry + 1
    }
  }
}

save(out,file="simulation2.Rdata")


## rerun failed simulations
load("simulation2.Rdata")
failed <- c()
for(i in 1:nsim) {
  if(is.null(out[[i]][[1]])) {
    cat("simulation",i,"of model 1 failed\n")
    failed <- c(failed,i)
  }
  if(is.null(out[[i]][[2]])) {
    cat("simulation",i,"of model 2 failed\n")
    failed <- c(failed,i)
  }
}
failed <- unique(failed)

set.seed(202109012)
randomSeeds <- sample(seq(nsim*1000,nsim*1000 + nsim*10000),nsim)

for(sim in failed) {
  set.seed(randomSeeds[sim])
  dat <- generate_data()
  
  set.seed(randomSeeds[sim])
  mod <- depmix(list(response~1),family=list(gaussian()),data=dat,nstates=3,ntimes=rep(nt,nrep))
  mod <- setpars(mod,truepars1)
  ok <- FALSE
  nTry <- 1
  while(!ok && nTry < 21) {
    fmod <- try(fit(mod,emcontrol=em.control(maxit = 5000, random.start = nTry > 1),verbose=FALSE))
    if(!inherits(fmod,"try-error")) {
      if(fmod@message == "Log likelihood converged to within tol. (relative change)") {
        cat("mod1 simulation",sim,"iteration",nTry,"completed.\n")
        ok <- TRUE
        out[[sim]][[1]] <- list(pars=getpars(fmod),logLik=logLik(fmod),viterbi=posterior(fmod, type="global"),trueState=dat$trueState)
      } else {
        cat("mod1 simulation",sim,"iteration",nTry,"failed.\n")
        nTry <- nTry + 1
      }
    } else {
      cat("mod1 simulation",sim,"iteration",nTry,"failed.\n")
      nTry <- nTry + 1
    }
  }
  
  set.seed(randomSeeds[sim])
  mod <- depmix(list(response~1,missing~1),family=list(gaussian(),multinomial("identity")),data=dat,nstates=3,ntimes=rep(nt,nrep))
  mod <- setpars(mod,truepars2)
  ok <- FALSE
  nTry <- 1
  while(!ok && nTry < 21) {
    fmod <- try(fit(mod,emcontrol=em.control(maxit = 5000, random.start = nTry > 1),verbose=FALSE))
    if(!inherits(fmod,"try-error")) {
      if(fmod@message == "Log likelihood converged to within tol. (relative change)") {
        cat("mod2 simulation",sim,"iteration",nTry,"completed.\n")
        ok <- TRUE
        out[[sim]][[2]] <- list(pars=getpars(fmod),logLik=logLik(fmod),viterbi=posterior(fmod, type="global"),trueState=dat$trueState)
      } else {
        cat("mod2 simulation",sim,"iteration",nTry,"failed.\n")
        nTry <- nTry + 1
      }
    } else {
      cat("mod2 simulation",sim,"iteration",nTry,"failed.\n")
      nTry <- nTry + 1
    }
  }
}
save(out,file="simulation2.Rdata")