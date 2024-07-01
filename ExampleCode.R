#################################################################
### Simulate data from a joint model for the marker evolution ###
### the visiting process (gap-time), and two competing risks  ###
#################################################################

# Set the working directory
setwd("F:/works/Methodological/VisitingProcessJMnew/pdfs/Biostatistics/Rev1")

# Set the seed
set.seed(2)

# Maximum study duration
tmax <- 10

# Fixed-effect parameters for the linear mixed model modeling the marker evolution
betaMark = c(17.20,4.83,-2.8)

if (!require("mvtnorm")) install.packages("mvtnorm")

# Dtrue covariance matrix of random effects
Dtrue = matrix(c(43.20,-7.56,-2.69,-7.56,7.68,-7.38,-2.69,-7.38,31.90),nr = 3, nc = 3)

# Within-individual standard deviation
sde = 2.81

# Number of individuals
N <- 1000

# Marker random effects
b <- rmvnorm(N,sigma = Dtrue)

# Parameters of the visiting-process model (gap time parameterization)
betaVis <- c(-1.15,0.02,-0.02,-1.50,0.02,0.20)

# Lists with the visit times and the observed marker values at the visit times
yList <- list()
timesList <- list()

# Parameters of the baseline hazard of the visiting model
basVisPar1 <- log(1.2)
basVisPar2 <- 0.85
basVisPar3 <- 10

# To center previous gap times
vstar <- 0.15

# log(time+cc)
cc <- 1

for (i in 1:N)
{
  # Simulate the marker value at baseline (t=0)
  mi0 <- c(cbind(1,log(0+cc),(0/10)^3) %*% (betaMark + b[i,]))
  yList[[i]] <- rnorm(1,mean = mi0,sd = sde)
  timesList[[i]] <- 0 
  
  # Gap times
  diffTime <- 0
  
  # Simulation until reaching the maximum study duration
  while( tail(timesList[[i]],1) < tmax)
  {
    # Start calendar time
    tij <- tail(timesList[[i]],1)
    
    # Marker value at the start of the interval
    yij <- tail(yList[[i]],1)
    
    # Occassion
    ord <- length(timesList[[i]])
    
    # Hazard function of time between visits
    hazGaptime <- function(s)
    {
      # s: time between visits
      # t: calendar time since origin
      t <- tij + s
      xit <- cbind(1,log(t+cc),(t/10)^3 )
      zit <- xit
      mit <- c(xit %*% betaMark + zit %*% b[i,])
      mitDer <- c(cbind(0,1/(t+cc),3*(t^2)/10^3 ) %*% (betaMark+b[i,]))
      
      # Baseline hazard
      basH = basVisPar3*dlnorm(s,meanlog = basVisPar1,sdlog = basVisPar2 )
      
      basH*exp(betaVis[1] + betaVis[2]*yij + betaVis[3]*tij + betaVis[4]*(diffTime-vstar)*(ord>1)
               + betaVis[5]*mi0 + betaVis[6]*mitDer)
    }
    
    # Survival function of the time between visits
    SurvGaptime <- function(time)
    {
      # Cumulative hazard function
      chaz <- integrate(hazGaptime,lower = 0,upper = time)$val
      
      # Survival function
      exp(-chaz)
    }
    
    # Inverse CDF method (to be performed numerically)
    ff <- function(time)
    {
      log(SurvGaptime(time)) - log(u)
    }
    ff <- Vectorize(ff)
    
    u <- runif(1)
    # Simulation of survival time using inverse CD4 theorem
    fit = try( uniroot(ff,lower = .Machine$double.eps,upper = tmax,tol = .Machine$double.eps),silent = T)
    
    if ( class(fit) %in% "try-error")
    {
      diffTime <- tmax
    } else{
      diffTime <- fit$root
    }
    
    # Time points at which marker measurements are to be collected
    timesList[[i]] <- c(timesList[[i]],tij+diffTime)
    
    # Simulate marker measurement
    xit <- cbind(1, log(tail(timesList[[i]],1)+cc) ,(tail(timesList[[i]],1)/10)^3 )
    zit <- xit
    mit <- c(xit %*% betaMark + zit %*% b[i,])
    yList[[i]] <- c(yList[[i]],rnorm(1,mean = mit,sd = sde))
  }
  # Remove data after study termination 
  timesList[[i]] <- timesList[[i]][-length(timesList[[i]])]
  yList[[i]] <- yList[[i]][-length(yList[[i]])]
  
  print(i)
  
}

# Longitudinal data.frame
dataFull <- data.frame( y = unlist(yList),times = unlist(timesList),
                        id = rep(1:N,sapply(yList,length)))

dataFull$lag1_y <- unlist(tapply(dataFull$y,dataFull$id,FUN = function(x){
  c(0,x[-length(x)])
}))
dataFull$lag1_gap <- unlist(tapply(dataFull$times,dataFull$id,function(x) {c(0,diff(x)-vstar)}))

###########################################################
### Simulation of survival time  (competing risks)      ###
### Two competing risks with the cause-specific hazards ###
### depending on previous marker values and             ###
### `true` current marker value and slope               ###
###########################################################

# Add a binary group
dataFull$group <- rbinom(N,size = 1,prob = 0.5)[dataFull$id]

# Parameters of the cause-specific hazards for cause 1
betasDrop1 <- c(-1.70,-0.20,0.5)
assocDrop1 <- c(-0.05,-0.10)

basH1Drop_1 <- log(9)
basH2Drop_1 <- 0.15

# Parameters of the cause-specific hazards for cause 2
betasDrop2 <- c(-4.66,-0.02,1.4,0.5)
assocDrop2 <- c(-0.05,-0.20)

basH1Drop_2 <- log(18)
basH2Drop_2 <- 0.55

# Number of observation by individual
ni <- table(dataFull$id)

# Order of measurements
dataFull$ord <- unlist(tapply(dataFull$id,dataFull$id,function(x){1:length(x)}))

# Survival dataset (one observation per individual)
dataid <- dataFull[dataFull$ord==1,]

# Observed survival time and failure indicator
dataid$time <- NA
dataid$status <- NA

# Function evaluating the cause-specific hazard for cause 1
HazDrop1 = function(time)
{
  # Exponentiated linear predictor at all time points
  exbetas1 <- c(exp(xMat1 %*% betasDrop1))
  
  # Time interval of time
  int <- findInterval(time,timePoints)
  
  # Baseline hazard
  basDrop1 <- exp(basH1Drop_1 - basH2Drop_1*time)
  
  # True marker value and rate of change at time
  mi0 <- c(cbind(1,log(0+cc),(0/10)^3 ) %*% (betaMark + b[i,]))
  mitDer <- c(cbind(0,1/(time+cc),3*(time^2)/10^3 ) %*% (betaMark+b[i,]))
  
  (exbetas1[int]*basDrop1)*exp(assocDrop1[1]*mi0 + assocDrop1[2]*mitDer)
  
}

# Function evaluating the cause-specific hazard for cause 2
HazDrop2 = function(time)
{
  # Exponentiated linear predictor at all time points
  exbetas2 <- c(exp(xMat2 %*% betasDrop2))
  
  # Time interval of time
  int <- findInterval(time,timePoints)
  
  # Baseline hazard
  basDrop2 <- exp(basH1Drop_2 - basH2Drop_2*time)
  
  # True marker value and rate of change at time
  mi0 <- c(cbind(1,log(0+cc),(0/10)^3 ) %*% (betaMark + b[i,]))
  mitDer <- c(cbind(0,1/(time+cc),3*(time^2)/10^3 ) %*% (betaMark+b[i,]))
  
  (exbetas2[int]*basDrop2)*exp(assocDrop2[1]*mi0 + assocDrop2[2]*mitDer)
  
}

# Function evaluating the all-cause survival function
SurvFun = function(time)
{
  # Number of time points on individual i
  nii <- length(timePoints)
  
  # Exponentiated linear predictor at all time points
  exbetas1 <- c(exp(xMat1 %*% betasDrop1))
  exbetas2 <- c(exp(xMat2 %*% betasDrop2))
  
  # Time interval of time
  int <- findInterval(time,timePoints)
  
  # Calculation of integrals of the hazard function for all relevant integrals
  integrals <- rep(NA,int)
  for (j in 1:int)
  {
    ff = function(x)
    {
      basDrop1 <- exp(basH1Drop_1 - basH2Drop_1*x)
      basDrop2 <- exp(basH1Drop_2 - basH2Drop_2*x)
      mi0 <- c(cbind(1,log(0+cc),(0/10)^3) %*% (betaMark + b[i,]))
      mitDer <- c(cbind(0,1/(x+cc),3*(x^2)/10^3) %*% (betaMark+b[i,]))
      
      temp <- (basDrop1*exbetas1[j])*exp(assocDrop1[1]*mi0 + assocDrop1[2]*mitDer)
      temp <- temp + (basDrop2*exbetas2[j])*exp(assocDrop2[1]*mi0 + assocDrop2[2]*mitDer)
      temp
    }
    
    # Upper time point
    if (j<int)
    {
      up <- timePoints[j+1]
    } else{
      up <- time
    }
    
    integrals[j] <- integrate(ff,lower = timePoints[j],upper = up  )$val
  }
  integrals <- integrals
  
  # Cumulative hazard function at time
  chaz <- sum(integrals)
  
  # Survival function
  exp(-chaz)
  
}

# Draws from the uniform distribution on (0,1) (for simulation using inversion)
u = runif(N)
table(u<.Machine$double.eps)
u[u<.Machine$double.eps] = .Machine$double.eps

ff = function(time)
{
  log(SurvFun(time)) - log(u[i])
}
ff <- Vectorize(ff)

for (i in 1:N)
{
  # Time points for the ith individual
  timePoints = dataFull$times[dataFull$id==i]
  
  # Design matrix for the ith individual
  xMat1 = cbind(1,dataFull$y[dataFull$id==i],dataFull$group[dataFull$id==i])
  xMat2 = cbind(1,dataFull$y[dataFull$id==i],dataFull$lag1_gap[dataFull$id==i],dataFull$group[dataFull$id==i])
  
  # Simulation of survival time using inverse CD4 theorem
  fit = try( uniroot(ff,lower = .Machine$double.eps,upper = tmax,tol = .Machine$double.eps),silent = T)
  
  if ( class(fit) %in% "try-error")
  {
    # Error denotes right censoring
    dataid$status[i] <- 0
    dataid$time[i] <- tmax
  } else{
    # Survival time observed (before tmax)
    if (fit$root==0)
    {
      fit$root <- .Machine$double.eps
    }
    dataid$time[i] <- fit$root
    
    ##################################
    ### Simulate the failure cause ###
    ##################################
    pcause1 <- HazDrop1(fit$root)
    pcause1 <- pcause1/(pcause1+HazDrop2(fit$root))
    
    # Simulate the cause
    cause1 = rbinom(n = 1,size = 1,prob = pcause1)
    dataid$status[i] = 1*cause1 + 2*(1-cause1)
    
  }
  print(i)
}

# Create missingness indicator
dataFull$mis = F
dataFull$N = rep(ni,ni)

# Now we need to delete the missing measurements 
# from the "longitudinal" data frame
etimes = rep(dataid$time,ni)
dataFull$mis = dataFull$times > etimes

# Observed data
dim(dataFull)
dataLong = dataFull[dataFull$mis == F,]
dim(dataLong)

dataLong$time = dataid$time[dataLong$id]
dataLong$status = dataid$status[dataLong$id]

ni = table(dataLong$id)
dataLong$N = rep(ni,ni)

############################
# Create failure indicator #
############################
dataLong$deltaCause1 <- 1*(dataLong$status==1)
dataLong$deltaCause1[duplicated(dataLong$id,fromLast = T)] <- 0
dataLong$deltaCause2 <- 1*(dataLong$status==2)
dataLong$deltaCause2[duplicated(dataLong$id,fromLast = T)] <- 0
ids = which(ni!=1)

######################
# Create entry times #
######################

tstart = dataLong$times

#####################
# Create exit times #
#####################

tstop = split(tstart,dataLong$id)
time = split(dataLong$time,dataLong$id)

for (i in ids)
{
  tstop[[i]] = c(tstop[[i]][2:ni[i]], time[[i]][ni[i]])
}

for (i in (1:N)[-ids])
{
  tstop[[i]] = time[[i]]
}

tstop = unlist(tstop,use.names = F)

dataLong$tstart <- tstart
dataLong$tstop <- tstop
rm(time,ids)

##################################################
### Gap time model for the time between visits ###
##################################################
dataLong$gapTimesVis <- dataLong$tstop - dataLong$tstart
dataLong$deltaVis <- 1*(dataLong$ord<dataLong$N)

######################
### Fit the models ###
######################

# Import the function
source("PropJMvisDropCompMLE_fit.R")

if (!require("nlme")) install.packages("nlme")
if (!require("survival")) install.packages("survival")
if (!require("matrixcalc")) install.packages("matrixcalc")
if (!require("splines2")) install.packages("splines2")
if (!require("MASS")) install.packages("MASS")

################################################################
### Fit the proposed model using a gap-time visiting process ###
################################################################

# Fit the linear mixed model
fitlme = try(lme(y ~ I(log(times+1)) + I((times/10)^3),
                 random = ~ I(log(times+1)) + I((times/10)^3)|id,data = dataLong,
                 control = list(apVar = T,opt = "optim" , returnObject = F, 
                                maxIter = 100, msMaxIter = 100, niterEM = 150)),
             silent = T)
summary(fitlme)

# Cox model for the gap times between visits
fitCoxVisGap <- coxph(Surv(gapTimesVis,deltaVis) ~ y + times + lag1_gap + cluster(id), 
                   data = dataLong, control = coxph.control(timefix = FALSE))
summary(fitCoxVisGap)

# Cox model for cause 1
fitCoxDrop1 <- coxph(Surv(tstart,tstop,deltaCause1) ~ y + group + cluster(id), data = dataLong,
                     control = coxph.control(timefix = FALSE))
summary(fitCoxDrop1)

# Cox model for cause 2
fitCoxDrop2 <- coxph(Surv(tstart,tstop,deltaCause2) ~ y + lag1_gap + group + cluster(id), data = dataLong,
                     control = coxph.control(timefix = FALSE))
summary(fitCoxDrop2)

##############################################
### Fit the proposed model using gap times ###
##############################################
fitPropGap <- PropJMvisDropCompMLE(fitlme,fitCoxDrop1,fitCoxDrop2,fitCoxVisGap,
                                   nknotsDrop1 = 1,nknotsDrop2 = 1,nknotsVis = 3,nGH = 4)
round(fitPropGap$sumLong,3)
round(fitPropGap$sumVis,3)
round(fitPropGap$sumDrop1,3)
round(fitPropGap$sumDrop2,3)


##################################################
### Fit the proposed model using calendar time ###
##################################################
fitCoxVisVcal <- coxph(Surv(tstart,tstop,deltaVis) ~ y + lag1_gap + cluster(id), 
                        data = dataLong,control = coxph.control(timefix = FALSE))
summary(fitCoxVisVcal)

# Fit the proposed model
fitPropCal <- PropJMvisDropCompMLE(fitlme, fitCoxDrop1,fitCoxDrop2, fitCoxVisVcal,
                                   nknotsDrop1 = 1,nknotsDrop2 = 1,nknotsVis = 3,nGH = 4)

round(fitPropCal$sumLong,3)
round(fitPropCal$sumVis,3)
round(fitPropCal$sumDrop1,3)
round(fitPropCal$sumDrop2,3)
