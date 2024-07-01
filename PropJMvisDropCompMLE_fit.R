###################################################################
### Function fitting the proposed approach for                  ###
### modeling longitudinal data, visit process and survival data ###
###################################################################

PropJMvisDropCompMLE <- function(fitlme,fitCoxDrop1,fitCoxDrop2,fitCoxVis,maxit = 50,epsilon = 1e-6,
                                 nknotsDrop1 = 3,nknotsDrop2 = 3,nknotsVis = 3,nGH = 10,alpha = 0.05,timeVar = "times",
                                 Gauleg_points = 30,useGauleg = F,epsDerXZ = 1e-7,ttBas = 0)
{
  # Needs a lme object for the longitudinal data
  # coxph objects for dropout data and time between visits
  
  if (!(class(fitlme)=="lme" & class(fitCoxDrop1)=="coxph" & class(fitCoxDrop2)=="coxph" & class(fitCoxVis)=="coxph"))
  {
    break
  }
  
  # Approximates the gradient of any function using forward difference
  fd <-
    function (x, f,n = length(x),..., eps = 1e-05) {
      #n <- length(x)
      res <- numeric(n)
      ex <- pmax(abs(x), 1)
      f0 <- f(x, ...)
      for (i in 1:n) {
        x1 <- x
        x1[i] <- x[i] + eps * ex[i]
        diff.f <- c(f(x1, ...) - f0)
        diff.x <- x1[i] - x[i]
        res[i] <- diff.f / diff.x
      }
      res
    }
  
  # Approximates the gradient of any function using forward difference
  fd.vec <-
    function (x, f, ..., eps = 1e-05) {
      n <- length(x)
      res <- matrix(NA,n,n)
      ex <- pmax(abs(x), 1)
      f0 <- f(x,n, ...)
      for (i in 1:n) {
        x1 <- x
        x1[i] <- x[i] + eps * ex[i]
        diff.f <- c(f(x1,n=i, ...) - f0[1:i])
        diff.x <- x1[i] - x[i]
        res[i,1:i] <- diff.f / diff.x
        res[1:i,i] <- diff.f / diff.x
      }
      res
    }
  
  fd.vec <-
    function (x, f, ..., eps = 1e-05) {
      n <- length(x)
      res <- matrix(0, n, n)
      ex <- eps * (abs(x) + eps)
      f0 <- f(x, ...)
      for (i in seq_len(n)) {
        x1 <- x
        x1[i] <- x[i] + ex[i]
        diff.f <- c(f(x1, ...) - f0)
        diff.x <- x1[i] - x[i]
        res[, i] <- diff.f / diff.x
      }
      0.5 * (res + t(res))
    }
  
  # Approximates the gradient of any function using central difference
  cd <-
    function (x, f, ..., eps = 1e-07) {
      n <- length(x)
      res <- numeric(n)
      ex <- pmax(abs(x), 1)
      for (i in 1:n) {
        x1 <- x2 <- x
        x1[i] <- x[i] + eps * ex[i]
        x2[i] <- x[i] - eps * ex[i]
        diff.f <- c(f(x1, ...) - f(x2, ...))
        diff.x <- x1[i] - x2[i]
        res[i] <- diff.f / diff.x
      }
      res
    }
  
  # Approximates the hessian of any function using central difference
  cd.vec <-
    function (x, f, ..., eps = 1e-04) {
      n <- length(x)
      res <- matrix(0, n, n)
      ex <- pmax(abs(x), 1)
      for (i in 1:n) {
        x1 <- x2 <- x
        x1[i] <- x[i] + eps * ex[i]
        x2[i] <- x[i] - eps * ex[i]
        diff.f <- c(f(x1, ...) - f(x2, ...))
        diff.x <- x1[i] - x2[i]
        res[, i] <- diff.f / diff.x
      }
      0.5 * (res + t(res))
    }
  
  # Function to evaluate the GH points
  gauher <- function(n) {# Gauss-Hermite: returns x,w so that
    #\int_-\infty^\infty exp(-x^2) f(x) dx \doteq \sum w_i f(x_i)
    EPS <- 3.e-14
    PIM4 <- .7511255444649425
    MAXIT <- 10
    m <- trunc((n+1)/2)
    x <- w <- rep(-1,n)
    for (i in 1:m) {
      if (i==1) {
        z <- sqrt(2*n+1)-1.85575*(2*n+1)^(-.16667)
      } else if(i==2) {
        z <- z-1.14*n^.426/z
      } else if (i==3) {
        z <- 1.86*z-.86*x[1]
      } else if (i==4) {
        z <- 1.91*z-.91*x[2]
      } else {
        z <- 2.*z-x[i-2]
      }
      for (its in 1:MAXIT) {
        p1 <- PIM4
        p2 <- 0
        for (j in 1:n) {
          p3 <- p2
          p2 <- p1
          p1 <- z*sqrt(2/j)*p2-sqrt((j-1)/j)*p3
        }
        pp <- sqrt(2*n)*p2
        z1 <- z
        z <- z1-p1/pp
        if(abs(z-z1) <= EPS) break
      }
      x[i] <- z
      x[n+1-i] <- -z
      w[i] <- 2/(pp*pp)
      w[n+1-i] <- w[i]
    }
    list(x=x,w=w)
  }
  
  mnormL = function(n,mu,L)
  {
    # n: number of draws
    # mu: mean vector
    # L: cholesky of precision matrix
    
    p = ncol(L)
    
    pivot = FALSE
    tol = -1
    
    # Upper triangular matrix
    z = t.default(sobol(n,dim = p,normal = T))
    
    sample = .Internal(backsolve(L, z, k = p, upper.tri = TRUE, transpose = FALSE)) + mu
    
    return(t.default(sample))
  }
  
  gauleg <- function(n,a=-1,b=1) {# Gauss-Legendre: returns x,w so that
    #\int_a^b f(x) dx \doteq \sum w_i f(x_i)
    EPS <- 3.e-14
    m <- trunc((n+1)/2)
    xm <- 0.5*(b+a)
    xl <- 0.5*(b-a)
    x <- w <- rep(-1,n)
    for (i in 1:m) {
      z <- cos(pi*(i-.25)/(n+.5))
      z1 <- z+1
      while (abs(z-z1) > EPS) {
        p1 <- 1
        p2 <- 0
        for (j in 1:n) {# recursively evaluates pn(x)
          p3 <- p2
          p2 <- p1
          p1 <- ((2*j-1)*z*p2-(j-1)*p3)/j
        }
        pp <- n*(z*p1-p2)/(z*z-1)
        z1 <- z
        z <- z1-p1/pp #Newton iteration
      }
      x[i] <- xm-xl*z
      x[n+1-i] <- xm+xl*z
      w[i] <- 2*xl/((1-z*z)*pp*pp)
      w[n+1-i] <- w[i]
    }
    list(x=x,w=w)
  }
  
  # 7-15 Gauss-Kronrod quadradure points and weights
  gaussKronrod <-
    function (k = 15) {
      sk <- c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, -0.405845151377397166906606412076961, 0,
              0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
              -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
              0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 0.991455371120812639206854697526329)
      wk15 <- c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 0.190350578064785409913256402421014,
                0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
                0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
                0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 0.022935322010529224963732008058970)
      wk7 <- c(0.129484966168869693270611432679082, 0.279705391489276667901467771423780, 0.381830050505118944950369775488975, 
               0.417959183673469387755102040816327, 0.381830050505118944950369775488975, 0.279705391489276667901467771423780, 0.129484966168869693270611432679082)
      if (k == 7) 
        list(sk = sk[1:7], wk = wk7)
      else
        list(sk = sk, wk = wk15)
    }
  
  #################################
  ### Prepare longitudinal data ###
  #################################
  
  # Vector with the observed marker values
  y <- as.vector(getResponse(fitlme))
  
  # Number of marker measurements
  Nrow <- length(y)
  
  # Get the data from the LME object
  dataLme <- fitlme$data
  
  # Design matrix of fixed effects
  xLong <- model.matrix(formula(fitlme),dataLme)
  
  # Dimension of fixed effects vector
  pLong <- ncol(xLong)
  
  xLong <- as.matrix(xLong[1:Nrow,1:pLong])
  colnames(xLong) = NULL
  rownames(xLong) = NULL
  
  # Design matrix of random effects
  zLong <- model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLme)
  
  # Number of random effects
  qLong <- ncol(zLong)
  
  zLong <- as.matrix(zLong[1:Nrow,1:qLong])
  colnames(zLong) = NULL
  rownames(zLong) = NULL
  
  # Grouping variable
  idLong <- as.numeric(as.character(fitlme$groups[,1]))

  # Number of subjects
  N <- length(unique(idLong))

  # Number of measurements by idLong
  ni <- table(idLong);names(ni) <- NULL
  
  # Indices of id's
  ids = list()
  for (i in 1:N){ids[[i]] = which(idLong==i)}
  
  # Split design matrices by idLong
  zLongList = lapply(split(data.frame(zLong),idLong),function(x) {row.names(x)=NULL;x = as.matrix(x)})
  zzLongList = lapply(zLongList,FUN=function(x) crossprod(x) ) # cross products
  xLongList = lapply(split(data.frame(xLong),idLong),function(x) {row.names(x) = NULL;x = as.matrix(x)})
  yList = split(y,idLong)
  for (i in 1:N){colnames(xLongList[[i]])<-NULL;colnames(zLongList[[i]])<-NULL}
  
  # Define vech (vector-half operator)
  vech <- function(x){t( t( x[!upper.tri(x)] )) }
  pivot <- FALSE
  tol <- -1
  
  # The duplication matrix transforming vech(A) to vec(A)
  Gqq = duplication.matrix(qLong)
  Kqq = commutation.matrix(r = qLong,c = qLong)
  
  # Prediction of the random effects based on the lme model
  brePred <- as.matrix(ranef(fitlme))
  rownames(brePred) <- NULL; colnames(brePred) <- NULL
  
  # Log-likelihood
  logl <- rep(NA,N)
  logVisDrop <- rep(NA,N)
  
  # Number of longitudinal model parameters
  nParLong <- pLong+1+qLong*(qLong+1)/2
  
  # Weights and quadrature points using the GH rule
  gh = gauher(nGH)

  ll1 = list()
  ll2 = list()
  for (j in 1:qLong)
  {
    ll1[[j]] = gh$x
    ll2[[j]] = gh$w
  }

  bgrid = as.matrix(expand.grid(ll1))
  tbgrid = t(bgrid)
  weights =  expand.grid(ll2)
  weights = apply(weights,1,prod)
  
  ############################################
  ### Prepare data for the dropout process ###
  ############################################
  nVis <- nrow(model.matrix(fitCoxVis))
  # Number of obs per subject in the dropout model
  temp = as.matrix(model.frame(fitCoxVis))
  idVis = as.numeric(temp[,ncol(temp)]);names(idVis)=NULL
  niVis = as.numeric(table(idVis))
  
  # Indices of id's
  idsVis = list()
  for (i in 1:N){idsVis[[i]] = which(idVis==i)}
  
  # Entry and exit time, and dropout indicator
  dVis <- get(as.character(fitCoxVis$call$data),envir = .GlobalEnv)
  if (ncol(fitCoxVis$y)==3)
  {
    tstartVis = fitCoxVis$y[,1]
    tstopVis = fitCoxVis$y[,2]
    deltaVis = fitCoxVis$y[,3]
    deltaVisSum <- tapply(deltaVis,idVis,sum)
    names(deltaVisSum) <- NULL
    names(tstartVis)=NULL;names(tstopVis)=NULL;names(deltaVis)=NULL
    
    # Calendar visit time (start of the interval)
    timesVis <- tstartVis
    calendarVisModel <- T
  } else if (ncol(fitCoxVis$y)==2)
  {
    tstartVis = rep(0,nVis)
    tstopVis = fitCoxVis$y[,1]
    deltaVis = fitCoxVis$y[,2]
    deltaVisSum <- tapply(deltaVis,idVis,sum)
    names(deltaVisSum) <- NULL
    names(tstartVis)=NULL;names(tstopVis)=NULL;names(deltaVis)=NULL
    
    # Calendar visit time (start of the interval)
    timesVis <- dVis[,timeVar]
    calendarVisModel <- F
  }
  
  # Design matrix for the dropout model (using all rows)
  xVis <- model.matrix(fitCoxVis)
  
  # Percentiles of the internal knots for B-splines
  centilesVis <- seq(0,1,length.out = nknotsVis + 2)[-c(1,nknotsVis + 2)]
  IntKnotsVis <- quantile(tstopVis[deltaVis==1],probs = centilesVis)
  bKnotsVis <- range(tstopVis)
  
  # B-splines basis matrices at visit times (tstop)
  BstVis <- bSpline(tstopVis,knots = IntKnotsVis,Boundary.knots = bKnotsVis,intercept = T)
  
  # Design matrices for the fixed and random effects at calendar visit times (tstop)
  #dataLmeid <- dataLme[rep(1,length(timesVis)),]
  dVis <- get(as.character(fitCoxVis$call$data),envir = .GlobalEnv)
  dataLmeid <- dVis
  dataLmeid[,timeVar] = ttBas
  
  xTimeVis = model.matrix(formula(fitlme),data = dataLmeid)
  colnames(xTimeVis) = NULL
  rownames(xTimeVis) = NULL
  
  zTimeVis = model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeid)
  colnames(zTimeVis) = NULL
  rownames(zTimeVis) = NULL
  
  # Prediction of the true marker values at calendar visit times (tstop)
  mtPredVis <- c(xTimeVis %*% fixed.effects(fitlme) + rowSums(zTimeVis*brePred[idVis,]))
  
  dataLmeid[,timeVar] = tstopVis*calendarVisModel + (timesVis + tstopVis)*(1-calendarVisModel)
  # Derivatives of the design matrices for the fixed and random effects at the observed survival times
  dataLmeidFOW <- dataLmeid
  dataLmeidFOW[,timeVar] <- dataLmeid[,timeVar] + epsDerXZ
  
  dataLmeidBAC <- dataLmeid
  dataLmeidBAC[,timeVar] <- dataLmeid[,timeVar] - epsDerXZ
  
  # Check for negative times
  aa <- which(dataLmeidBAC[,timeVar]<0)
  dataLmeidBAC[aa,timeVar] <- dataLmeid[,timeVar][aa]
  
  xTimeFOWVis <- model.matrix(formula(fitlme),data = dataLmeidFOW)
  zTimeFOWVis <- model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidFOW)
  xTimeBACVis <- model.matrix(formula(fitlme),data = dataLmeidBAC)
  zTimeBACVis <- model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidBAC)
  
  xTimeDerVis <- (xTimeFOWVis-xTimeBACVis)/(dataLmeidFOW[,timeVar]-dataLmeidBAC[,timeVar])
  zTimeDerVis <- (zTimeFOWVis-zTimeBACVis)/(dataLmeidFOW[,timeVar]-dataLmeidBAC[,timeVar])
  colnames(xTimeDerVis) <- NULL;rownames(xTimeDerVis) <- NULL
  colnames(zTimeDerVis) <- NULL;rownames(zTimeDerVis) <- NULL
  
  rm(list = c("dataLmeidFOW","dataLmeidBAC","xTimeFOWVis","zTimeFOWVis","xTimeBACVis","zTimeBACVis"))
  mtDerVis <- c(xTimeDerVis %*% fixed.effects(fitlme) + rowSums(zTimeDerVis*brePred[idVis,]))
  
  # Gauss-Legendre rules (required for integration of the hazard function)
  if (useGauleg == F)
  {
    gk = gaussKronrod()
    gk <- data.frame(sk = gk$sk,wk = gk$wk)
    gk <- gk[order(gk$sk),]
    pps = gk$sk
    wws = gk$wk
  } else {
    gk = gauleg(Gauleg_points)
    pps = gk$x
    wws = gk$w
  }
  
  # Additional time points for integration
  p1Vis <- rep( (tstartVis + tstopVis)/2,each = length(pps))
  p2Vis <- rep( (tstopVis - tstartVis)/2,each = length(pps))
  ppsVisLong = rep(pps,nrow(xVis))
  wwsVisLong = rep(wws,nrow(xVis))
  indVis = rep(1:nVis,each = length(pps))
  timesVisForInt <- p2Vis*ppsVisLong + p1Vis
  
  # Design matrix for the survival model (replicated rows for integration)
  xVisForInt <- xVis[indVis,];if(ncol(xVis)==1){xVisForInt <- as.matrix(xVisForInt)}
  
  # B-spline basis matrix (replicated rows for integration)
  BstVisForInt <- bSpline(timesVisForInt,knots = IntKnotsVis,Boundary.knots = bKnotsVis,intercept = T)
  
  # Design matrices of the fixed and random effects at additional time points for integration
  #dataLmeidForInt <- dataLmeid[rep(1,length(timesVisForInt)),]
  dataLmeidForInt <- dVis[indVis,]
  dataLmeidForInt[,timeVar] <- ttBas
  idVisInd <- rep(idVis,each = length(pps))
  
  # Indices of id's for integration the hazard function of the dropout process
  idsVisForInt = list()
  for (i in 1:N){idsVisForInt[[i]] = which(idVisInd==i)}
  
  xTimeVisForInt <- model.matrix(formula(fitlme),data = dataLmeidForInt)
  rownames(xTimeVisForInt) <- NULL
  colnames(xTimeVisForInt) = NULL
  
  zTimeVisForInt = model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidForInt)
  rownames(zTimeVisForInt) <- NULL
  colnames(zTimeVisForInt) <- NULL
  
  # True marker values at the additional time points
  mtPredVisForInt <- c(xTimeVisForInt %*% fixed.effects(fitlme) + 
                         rowSums(zTimeVisForInt*(brePred[idVis,][indVis,])))
  
  dataLmeidForInt[,timeVar] <- timesVisForInt + rep(timesVis,each = length(pps))*(1-calendarVisModel)
  # Derivative of true marker values at the additional time points
  dataLmeidForIntFOW <- dataLmeidForInt
  dataLmeidForIntFOW[,timeVar] <- dataLmeidForInt[,timeVar] + epsDerXZ
  
  dataLmeidForIntBAC <- dataLmeidForInt
  dataLmeidForIntBAC[,timeVar] <- dataLmeidForInt[,timeVar] - epsDerXZ
  
  # Check for negative times
  aa <- which(dataLmeidForIntBAC[,timeVar] < 0)
  dataLmeidForIntBAC[aa,timeVar] <- timesVisForInt[aa]
  
  xTimeVisForIntFOW <- model.matrix(formula(fitlme),data = dataLmeidForIntFOW)
  zTimeVisForIntFOW <- model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidForIntFOW)
  xTimeVisForIntBAC <- model.matrix(formula(fitlme),data = dataLmeidForIntBAC)
  zTimeVisForIntBAC <- model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidForIntBAC)
  
  xTimeDerVisForInt <- (xTimeVisForIntFOW - xTimeVisForIntBAC)/(dataLmeidForIntFOW[,timeVar]-
                                                                  dataLmeidForIntBAC[,timeVar])
  zTimeDerVisForInt <- (zTimeVisForIntFOW - zTimeVisForIntBAC)/(dataLmeidForIntFOW[,timeVar]-
                                                                  dataLmeidForIntBAC[,timeVar])
  colnames(xTimeDerVisForInt) <- NULL;rownames(xTimeDerVisForInt) <- NULL;
  colnames(zTimeDerVisForInt) <- NULL;rownames(zTimeDerVisForInt) <- NULL;
  
  rm(list = c("dataLmeidForIntFOW","dataLmeidForIntBAC","xTimeVisForIntFOW","zTimeVisForIntFOW",
              "xTimeVisForIntBAC","zTimeVisForIntBAC"))
  mtDerVisForInt <- c(xTimeDerVisForInt %*% fixed.effects(fitlme) + 
                        rowSums(zTimeDerVisForInt*(brePred[idVis,][indVis,])))
  
  
  tstartVisForInt <- timesVisForInt
  tstopVisForInt <- unlist(tapply(tstartVisForInt,idVisInd,function(x) c(x[2:length(x)],-1)  ))
  tstopVisForInt[!duplicated(indVis,fromLast = T)] <- tstopVis
  deltaVisForInt <- rep(0,length(timesVisForInt))
  deltaVisForInt[!duplicated(indVis,fromLast = T)] <- deltaVis
  
  # Starting values for the parameters of the dropout process
  # Covariates + predictions of the random effects based on the lme
  fitVisStart <- coxph( Surv(tstartVisForInt,tstopVisForInt,deltaVisForInt) ~ xVisForInt + 
                          mtPredVisForInt + mtDerVisForInt,
                        control = coxph.control(timefix = FALSE))
  summary(fitVisStart)
  
  # Obtain starting values for the B-spline parameters
  # by approximating the baseline cumulative based on the Cox model
  predChaz = basehaz(fitVisStart,centered = F)
  predChaz = predChaz[predChaz$hazard!=0,]
  predChaz$diff <- c(NA,diff(predChaz$hazard))
  predChaz <- predChaz[which(predChaz$diff!=0),]
  predChaz$logdiff <- log(predChaz$diff)
  predChaz <- predChaz[is.finite(predChaz$logdiff),]
  
  # B-spline basis matrix at the times where the hazard has been evaluated
  splVis <- bSpline(predChaz$time,knots = IntKnotsVis,Boundary.knots = bKnotsVis,intercept = T)
  
  # Starting values for the B-spline parameters
  fit <- lm(predChaz$logdiff ~ -1 + splVis)
  
  lbVis <- length(coef(fitVisStart))
  lassocVis <- lbVis - ncol(xVis)
  lpsiVis <- ncol(BstVis)
  lVis <- lbVis + lpsiVis
  paramVis <- c(coef(fitVisStart),coef(fit))
  names(paramVis) <- NULL
  paramVis_bup <- paramVis
  
  # Remove unnecessary objects
  rm(list = c("predChaz","fit","splVis","tstartVisForInt","tstopVisForInt",
              "deltaVisForInt","dataLmeidForInt"))
  
  fVis = function(paramVis)
  {
    # Parameters of the dropout model (regarding the observed marker values and other covariates)
    betasVis <- paramVis[1:ncol(xVis)]
    
    # Association parameters (true marker value and true marker slope)
    assocVis <- paramVis[(ncol(xVis)+1):lbVis]
    
    # Parameters for the B-spline basis matrix
    psiVis <- paramVis[-(1:lbVis)]
    
    # Log hazard at the visit times (tstop)
    loghazard <- c(xVis %*% betasVis + BstVis %*% psiVis + assocVis[1]*mtPredVis+
                     assocVis[2]*mtDerVis)
    
    # Log hazard at additional time points (for integration)
    loghazardForInt <- c((xVisForInt %*% betasVis) + (BstVisForInt %*% psiVis)) + (
      assocVis[1]*mtPredVisForInt + assocVis[2]*mtDerVisForInt)
    
    # Integration of the hazard
    hazInt <- p2Vis*(exp(loghazardForInt)*wwsVisLong)
    
    # Log likelihood
    logl <- sum(loghazard*deltaVis) - sum(hazInt)
    
    -logl
  }
  
  # Gradient of the log-likelihood for the dropout process
  fdGradVis <- function(x,n=length(x)){fd(x,fVis,n)}
  cdGradVis <- function(x){cd(x,fVis)}
  
  fitContrVis <- try(optim(paramVis,fn = fVis,gr = cdGradVis,method = "BFGS",
                           control = list(trace = 0,maxit = 300)),silent = TRUE)
  if (!"try-error" %in% class(fitContrVis))
  {
    paramVis <- fitContrVis$par
  } else{paramVis <- paramVis_bup}
  
  ######################################################
  ### Prepare data for the dropout process (cause 1) ###
  ######################################################
  nDrop1 <- nrow(model.matrix(fitCoxDrop1))
  # Number of obs per subject in the dropout model
  temp = as.matrix(model.frame(fitCoxDrop1))
  idDrop1 = as.numeric(temp[,ncol(temp)]);names(idDrop1)=NULL
  niDrop1 = as.numeric(table(idDrop1))
  
  # Entry and exit time, and dropout indicator
  if (ncol(fitCoxDrop1$y)==3)
  {
    tstartDrop1 = fitCoxDrop1$y[,1]
    tstopDrop1 = fitCoxDrop1$y[,2]
    deltaDrop1 = fitCoxDrop1$y[,3]
    deltaDrop1Sum <- tapply(deltaDrop1,idDrop1,sum)
    names(deltaDrop1Sum) <- NULL
    names(tstartDrop1)=NULL;names(tstopDrop1)=NULL;names(deltaDrop1)=NULL
  } else if (ncol(fitCoxDrop1$y)==2)
  {
    tstartDrop1 = rep(0,nDrop1)
    tstopDrop1 = fitCoxDrop1$y[,1]
    deltaDrop1 = fitCoxDrop1$y[,2]
    deltaDrop1Sum <- tapply(deltaDrop1,idDrop1,sum)
    names(deltaDrop1Sum) <- NULL
    names(tstartDrop1)=NULL;names(tstopDrop1)=NULL;names(deltaDrop1)=NULL
  }
  
  # Observed survival times 
  timesDrop1 <- tstopDrop1[!duplicated(idDrop1,fromLast = T)]
  
  # Design matrix for the dropout model (using all rows)
  xDrop1 <- model.matrix(fitCoxDrop1)
  
  # Design matrix for the dropout model at the last observation
  xDrop1Final <- xDrop1[!duplicated(idDrop1,fromLast = T),];if(ncol(xDrop1)==1){xDrop1Final <- as.matrix(xDrop1Final)}
  
  # Percentiles of the internal knots for B-splines
  centilesDrop1 <- seq(0,1,length.out = nknotsDrop1 + 2)[-c(1,nknotsDrop1 + 2)]
  IntKnotsDrop1 <- quantile(tstopDrop1[deltaDrop1==1],probs = centilesDrop1)
  bKnotsDrop1 <- range(tstopDrop1)
  
  # B-splines basis matrices at survival times
  BstDrop1Final <- bSpline(timesDrop1,knots = IntKnotsDrop1,Boundary.knots = bKnotsDrop1,intercept = T)

  # Design matrices for the fixed and random effects at the observed survival times
  dataLmeid <- dataLme[!duplicated(idLong),]
  dataLmeid[,timeVar] = ttBas
  
  xTimeDrop1 = model.matrix(formula(fitlme),data = dataLmeid)
  colnames(xTimeDrop1) = NULL
  rownames(xTimeDrop1) = NULL
  
  zTimeDrop1 = model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeid)
  colnames(zTimeDrop1) = NULL
  rownames(zTimeDrop1) = NULL
  
  # Prediction of the true marker values at the observed survival times
  mtPredDrop1Final <- c(xTimeDrop1 %*% fixed.effects(fitlme) + rowSums(zTimeDrop1*brePred))
  
  dataLmeid[,timeVar] = timesDrop1
  # Derivatives of the design matrices for the fixed and random effects at the observed survival times
  dataLmeidFOW <- dataLmeid
  dataLmeidFOW[,timeVar] <- timesDrop1 + epsDerXZ
  
  dataLmeidBAC <- dataLmeid
  dataLmeidBAC[,timeVar] <- timesDrop1 - epsDerXZ
  
  # Check for negative times
  aa <- which(dataLmeidBAC[,timeVar]<0)
  dataLmeidBAC[aa,timeVar] <- timesDrop1[aa]
  
  xTimeFOWDrop1 <- model.matrix(formula(fitlme),data = dataLmeidFOW)
  zTimeFOWDrop1 <- model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidFOW)
  xTimeBACDrop1 <- model.matrix(formula(fitlme),data = dataLmeidBAC)
  zTimeBACDrop1 <- model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidBAC)
  
  xTimeDerDrop1 <- (xTimeFOWDrop1-xTimeBACDrop1)/(dataLmeidFOW[,timeVar]-dataLmeidBAC[,timeVar])
  zTimeDerDrop1 <- (zTimeFOWDrop1-zTimeBACDrop1)/(dataLmeidFOW[,timeVar]-dataLmeidBAC[,timeVar])
  colnames(xTimeDerDrop1) <- NULL;rownames(xTimeDerDrop1) <- NULL
  colnames(zTimeDerDrop1) <- NULL;rownames(zTimeDerDrop1) <- NULL
  
  rm(list = c("dataLmeidFOW","dataLmeidBAC","xTimeFOWDrop1","zTimeFOWDrop1","xTimeBACDrop1","zTimeBACDrop1"))
  mtDerDrop1Final <- c(xTimeDerDrop1 %*% fixed.effects(fitlme) + rowSums(zTimeDerDrop1*brePred))
  
  # Gauss-Legendre rules (required for integration of the hazard function)
  if (useGauleg == F)
  {
    gk = gaussKronrod()
    gk <- data.frame(sk = gk$sk,wk = gk$wk)
    gk <- gk[order(gk$sk),]
    pps = gk$sk
    wws = gk$wk
  } else {
    gk = gauleg(Gauleg_points)
    pps = gk$x
    wws = gk$w
  }
  
  # Additional time points for integration
  p1Drop1 <- rep( (tstartDrop1 + tstopDrop1)/2,each = length(pps))
  p2Drop1 <- rep( (tstopDrop1 - tstartDrop1)/2,each = length(pps))
  ppsDrop1Long = rep(pps,nrow(xDrop1))
  wwsDrop1Long = rep(wws,nrow(xDrop1))
  indDrop1 = rep(1:nDrop1,each = length(pps))
  timesDrop1ForInt <- p2Drop1*ppsDrop1Long + p1Drop1

  # Design matrix for the survival model (replicated rows for integration)
  xDrop1ForInt <- xDrop1[indDrop1,];if(ncol(xDrop1)==1){xDrop1ForInt <- as.matrix(xDrop1ForInt)}
  
  # B-spline basis matrix (replicated rows for integration)
  BstDrop1ForInt <- bSpline(timesDrop1ForInt,knots = IntKnotsDrop1,Boundary.knots = bKnotsDrop1,intercept = T)
  
  # Design matrices of the fixed and random effects at additional time points for integration
  #dataLmeidForInt <- dataLmeid[rep(1,length(timesDrop1ForInt)),]
  dSurv1 <- get(as.character(fitCoxDrop1$call$data),envir = .GlobalEnv)
  dataLmeidForInt <- dSurv1[indDrop1,]
  dataLmeidForInt[,timeVar] <- ttBas
  idDrop1Ind <- rep(idDrop1,each = length(pps))
  
  # Indices of id's for integration the hazard function of the dropout process
  idsDrop1ForInt = list()
  for (i in 1:N){idsDrop1ForInt[[i]] = which(idDrop1Ind==i)}
  
  xTimeDrop1ForInt <- model.matrix(formula(fitlme),data = dataLmeidForInt)
  rownames(xTimeDrop1ForInt) <- NULL
  colnames(xTimeDrop1ForInt) = NULL
  
  zTimeDrop1ForInt = model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidForInt)
  rownames(zTimeDrop1ForInt) <- NULL
  colnames(zTimeDrop1ForInt) <- NULL
  
  # True marker values at the additional time points
  mtPredDrop1ForInt <- c(xTimeDrop1ForInt %*% fixed.effects(fitlme) + 
                        rowSums(zTimeDrop1ForInt*(brePred[idDrop1,][indDrop1,])))
  
  dataLmeidForInt[,timeVar] <- timesDrop1ForInt
  # Derivative of true marker values at the additional time points
  dataLmeidForIntFOW <- dataLmeidForInt
  dataLmeidForIntFOW[,timeVar] <- timesDrop1ForInt + epsDerXZ
  
  dataLmeidForIntBAC <- dataLmeidForInt
  dataLmeidForIntBAC[,timeVar] <- timesDrop1ForInt - epsDerXZ
  
  # Check for negative times
  aa <- which(dataLmeidForIntBAC[,timeVar] < 0)
  dataLmeidForIntBAC[aa,timeVar] <- timesDrop1ForInt[aa]
  
  xTimeDrop1ForIntFOW <- model.matrix(formula(fitlme),data = dataLmeidForIntFOW)
  zTimeDrop1ForIntFOW <- model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidForIntFOW)
  xTimeDrop1ForIntBAC <- model.matrix(formula(fitlme),data = dataLmeidForIntBAC)
  zTimeDrop1ForIntBAC <- model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidForIntBAC)
  
  xTimeDerDrop1ForInt <- (xTimeDrop1ForIntFOW - xTimeDrop1ForIntBAC)/(dataLmeidForIntFOW[,timeVar]-
                                                                   dataLmeidForIntBAC[,timeVar])
  zTimeDerDrop1ForInt <- (zTimeDrop1ForIntFOW - zTimeDrop1ForIntBAC)/(dataLmeidForIntFOW[,timeVar]-
                                                                   dataLmeidForIntBAC[,timeVar])
  colnames(xTimeDerDrop1ForInt) <- NULL;rownames(xTimeDerDrop1ForInt) <- NULL;
  colnames(zTimeDerDrop1ForInt) <- NULL;rownames(zTimeDerDrop1ForInt) <- NULL;
  
  rm(list = c("dataLmeidForIntFOW","dataLmeidForIntBAC","xTimeDrop1ForIntFOW","zTimeDrop1ForIntFOW",
              "xTimeDrop1ForIntBAC","zTimeDrop1ForIntBAC"))
  mtDerDrop1ForInt <- c(xTimeDerDrop1ForInt %*% fixed.effects(fitlme) + 
                       rowSums(zTimeDerDrop1ForInt*(brePred[idDrop1,][indDrop1,])))
  
  
  tstartDrop1ForInt <- timesDrop1ForInt
  tstopDrop1ForInt <- unlist(tapply(tstartDrop1ForInt,idDrop1Ind,function(x) c(x[2:length(x)],NA)  ))
  tstopDrop1ForInt[!duplicated(idDrop1Ind,fromLast = T)] <- timesDrop1
  deltaDrop1ForInt <- rep(0,length(timesDrop1ForInt))
  deltaDrop1ForInt[!duplicated(idDrop1Ind,fromLast = T)] <- deltaDrop1Sum
  
  # Starting values for the parameters of the dropout process
  # Covariates + predictions of the random effects based on the lme
  fitDrop1Start <- coxph( Surv(tstartDrop1ForInt,tstopDrop1ForInt,deltaDrop1ForInt) ~ xDrop1ForInt + 
                        mtPredDrop1ForInt + mtDerDrop1ForInt,
                         control = coxph.control(timefix = FALSE))
  summary(fitDrop1Start)
  
  
  # Obtain starting values for the B-spline parameters
  # by approximating the baseline cumulative based on the Cox model
  predChaz = basehaz(fitDrop1Start,centered = F)
  predChaz = predChaz[predChaz$hazard!=0,]
  predChaz$diff <- c(NA,diff(predChaz$hazard))
  predChaz <- predChaz[which(predChaz$diff!=0),]
  predChaz$logdiff <- log(predChaz$diff)
  predChaz <- predChaz[is.finite(predChaz$logdiff),]
  
  # B-spline basis matrix at the times where the hazard has been evaluated
  splDrop1 <- bSpline(predChaz$time,knots = IntKnotsDrop1,Boundary.knots = bKnotsDrop1,intercept = T)
  
  # Starting values for the B-spline parameters
  fit <- lm(predChaz$logdiff ~ -1 + splDrop1)
  
  lbDrop1 <- length(coef(fitDrop1Start))
  lassocDrop1 <- lbDrop1 - ncol(xDrop1)
  lpsiDrop1 <- ncol(BstDrop1Final)
  lDrop1 <- lbDrop1 + lpsiDrop1
  paramDrop1 <- c(coef(fitDrop1Start),coef(fit))
  names(paramDrop1) <- NULL
  paramDrop1_bup <- paramDrop1
  
  # Remove unnecessary objects
  rm(list = c("predChaz","fit","splDrop1","tstartDrop1ForInt","tstopDrop1ForInt",
              "deltaDrop1ForInt","dataLmeidForInt"))
  
  fDrop1 = function(paramDrop1)
  {
    # Parameters of the dropout model (regarding the observed marker values and other covariates)
    betasDrop1 <- paramDrop1[1:ncol(xDrop1)]
    
    # Association parameters (true marker value and true marker slope)
    assocDrop1 <- paramDrop1[(ncol(xDrop1)+1):lbDrop1]
    
    # Parameters for the B-spline basis matrix
    psiDrop1 <- paramDrop1[-(1:lbDrop1)]
    
    # Log hazard at the observed survival time
    loghazard <- c(xDrop1Final %*% betasDrop1 + BstDrop1Final %*% psiDrop1 + assocDrop1[1]*mtPredDrop1Final+
                   assocDrop1[2]*mtDerDrop1Final)
    
    # Log hazard at additional time points (for integration)
    loghazardForInt <- c((xDrop1ForInt %*% betasDrop1) + (BstDrop1ForInt %*% psiDrop1)) + (
                        assocDrop1[1]*mtPredDrop1ForInt + assocDrop1[2]*mtDerDrop1ForInt)
    
    # Integration of the hazard
    hazInt <- p2Drop1*(exp(loghazardForInt)*wwsDrop1Long)
    
    # Log likelihood
    logl <- sum(loghazard*deltaDrop1Sum) - sum(hazInt)
    
    -logl
  }
  
  # Gradient of the log-likelihood for the dropout process
  fdGradDrop1 <- function(x,n=length(x)){fd(x,fDrop1,n)}
  cdGradDrop1 <- function(x){cd(x,fDrop1)}
  
  fitContrDrop1 <- try(optim(paramDrop1,fn = fDrop1,gr = cdGradDrop1,method = "BFGS",
                            control = list(trace = 0,maxit = 300)),silent = TRUE)
  if (!"try-error" %in% class(fitContrDrop1))
  {
    paramDrop1 <- fitContrDrop1$par
  } else{paramDrop1 <- paramDrop1_bup}
  
  ######################################################
  ### Prepare data for the dropout process (cause 2) ###
  ######################################################
  nDrop2 <- nrow(model.matrix(fitCoxDrop2))
  # Number of obs per subject in the dropout model
  temp = as.matrix(model.frame(fitCoxDrop2))
  idDrop2 = as.numeric(temp[,ncol(temp)]);names(idDrop2)=NULL
  niDrop2 = as.numeric(table(idDrop2))
  
  # Entry and exit time, and dropout indicator
  if (ncol(fitCoxDrop2$y)==3)
  {
    tstartDrop2 = fitCoxDrop2$y[,1]
    tstopDrop2 = fitCoxDrop2$y[,2]
    deltaDrop2 = fitCoxDrop2$y[,3]
    deltaDrop2Sum <- tapply(deltaDrop2,idDrop2,sum)
    names(deltaDrop2Sum) <- NULL
    names(tstartDrop2)=NULL;names(tstopDrop2)=NULL;names(deltaDrop2)=NULL
  } else if (ncol(fitCoxDrop2$y)==2)
  {
    tstartDrop2 = rep(0,nDrop2)
    tstopDrop2 = fitCoxDrop2$y[,1]
    deltaDrop2 = fitCoxDrop2$y[,2]
    deltaDrop2Sum <- tapply(deltaDrop2,idDrop2,sum)
    names(deltaDrop2Sum) <- NULL
    names(tstartDrop2)=NULL;names(tstopDrop2)=NULL;names(deltaDrop2)=NULL
  }
  
  # Observed survival times 
  timesDrop2 <- tstopDrop2[!duplicated(idDrop2,fromLast = T)]
  
  # Design matrix for the dropout model (using all rows)
  xDrop2 <- model.matrix(fitCoxDrop2)
  
  # Design matrix for the dropout model at the last observation
  xDrop2Final <- xDrop2[!duplicated(idDrop2,fromLast = T),];if(ncol(xDrop2)==1){xDrop2Final <- as.matrix(xDrop2Final)}
  
  # Percentiles of the internal knots for B-splines
  centilesDrop2 <- seq(0,1,length.out = nknotsDrop2 + 2)[-c(1,nknotsDrop2 + 2)]
  IntKnotsDrop2 <- quantile(tstopDrop2[deltaDrop2==1],probs = centilesDrop2)
  bKnotsDrop2 <- range(tstopDrop2)
  
  # B-splines basis matrices at survival times
  BstDrop2Final <- bSpline(timesDrop2,knots = IntKnotsDrop2,Boundary.knots = bKnotsDrop2,intercept = T)
  
  # Design matrices for the fixed and random effects at the observed survival times
  dataLmeid <- dataLme[!duplicated(idLong),]
  dataLmeid[,timeVar] = ttBas
  
  xTimeDrop2 = model.matrix(formula(fitlme),data = dataLmeid)
  colnames(xTimeDrop2) = NULL
  rownames(xTimeDrop2) = NULL
  
  zTimeDrop2 = model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeid)
  colnames(zTimeDrop2) = NULL
  rownames(zTimeDrop2) = NULL
  
  # Prediction of the true marker values at the observed survival times
  mtPredDrop2Final <- c(xTimeDrop2 %*% fixed.effects(fitlme) + rowSums(zTimeDrop2*brePred))
  
  dataLmeid[,timeVar] = timesDrop2
  # Derivatives of the design matrices for the fixed and random effects at the observed survival times
  dataLmeidFOW <- dataLmeid
  dataLmeidFOW[,timeVar] <- timesDrop2 + epsDerXZ
  
  dataLmeidBAC <- dataLmeid
  dataLmeidBAC[,timeVar] <- timesDrop2 - epsDerXZ
  
  # Check for negative times
  aa <- which(dataLmeidBAC[,timeVar]<0)
  dataLmeidBAC[aa,timeVar] <- timesDrop2[aa]
  
  xTimeFOWDrop2 <- model.matrix(formula(fitlme),data = dataLmeidFOW)
  zTimeFOWDrop2 <- model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidFOW)
  xTimeBACDrop2 <- model.matrix(formula(fitlme),data = dataLmeidBAC)
  zTimeBACDrop2 <- model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidBAC)
  
  xTimeDerDrop2 <- (xTimeFOWDrop2-xTimeBACDrop2)/(dataLmeidFOW[,timeVar]-dataLmeidBAC[,timeVar])
  zTimeDerDrop2 <- (zTimeFOWDrop2-zTimeBACDrop2)/(dataLmeidFOW[,timeVar]-dataLmeidBAC[,timeVar])
  colnames(xTimeDerDrop2) <- NULL;rownames(xTimeDerDrop2) <- NULL
  colnames(zTimeDerDrop2) <- NULL;rownames(zTimeDerDrop2) <- NULL
  
  rm(list = c("dataLmeidFOW","dataLmeidBAC","xTimeFOWDrop2","zTimeFOWDrop2","xTimeBACDrop2","zTimeBACDrop2"))
  mtDerDrop2Final <- c(xTimeDerDrop2 %*% fixed.effects(fitlme) + rowSums(zTimeDerDrop2*brePred))
  
  # Gauss-Legendre rules (required for integration of the hazard function)
  if (useGauleg == F)
  {
    gk = gaussKronrod()
    gk <- data.frame(sk = gk$sk,wk = gk$wk)
    gk <- gk[order(gk$sk),]
    pps = gk$sk
    wws = gk$wk
  } else {
    gk = gauleg(Gauleg_points)
    pps = gk$x
    wws = gk$w
  }
  
  # Additional time points for integration
  p1Drop2 <- rep( (tstartDrop2 + tstopDrop2)/2,each = length(pps))
  p2Drop2 <- rep( (tstopDrop2 - tstartDrop2)/2,each = length(pps))
  ppsDrop2Long = rep(pps,nrow(xDrop2))
  wwsDrop2Long = rep(wws,nrow(xDrop2))
  indDrop2 = rep(1:nDrop2,each = length(pps))
  timesDrop2ForInt <- p2Drop2*ppsDrop2Long + p1Drop2
  
  # Design matrix for the survival model (replicated rows for integration)
  xDrop2ForInt <- xDrop2[indDrop2,];if(ncol(xDrop2)==1){xDrop2ForInt <- as.matrix(xDrop2ForInt)}
  
  # B-spline basis matrix (replicated rows for integration)
  BstDrop2ForInt <- bSpline(timesDrop2ForInt,knots = IntKnotsDrop2,Boundary.knots = bKnotsDrop2,intercept = T)
  
  # Design matrices of the fixed and random effects at additional time points for integration
  #dataLmeidForInt <- dataLmeid[rep(1,length(timesDrop2ForInt)),]
  dSurv2 <- get(as.character(fitCoxDrop2$call$data),envir = .GlobalEnv)
  dataLmeidForInt <- dSurv2[indDrop2,]
  dataLmeidForInt[,timeVar] <- ttBas
  idDrop2Ind <- rep(idDrop2,each = length(pps))
  
  # Indices of id's for integration the hazard function of the dropout process
  idsDrop2ForInt = list()
  for (i in 1:N){idsDrop2ForInt[[i]] = which(idDrop2Ind==i)}
  
  xTimeDrop2ForInt <- model.matrix(formula(fitlme),data = dataLmeidForInt)
  rownames(xTimeDrop2ForInt) <- NULL
  colnames(xTimeDrop2ForInt) = NULL
  
  zTimeDrop2ForInt = model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidForInt)
  rownames(zTimeDrop2ForInt) <- NULL
  colnames(zTimeDrop2ForInt) <- NULL
  
  # True marker values at the additional time points
  mtPredDrop2ForInt <- c(xTimeDrop2ForInt %*% fixed.effects(fitlme) + 
                           rowSums(zTimeDrop2ForInt*(brePred[idDrop2,][indDrop2,])))
  
  dataLmeidForInt[,timeVar] <- timesDrop2ForInt
  # Derivative of true marker values at the additional time points
  dataLmeidForIntFOW <- dataLmeidForInt
  dataLmeidForIntFOW[,timeVar] <- timesDrop2ForInt + epsDerXZ
  
  dataLmeidForIntBAC <- dataLmeidForInt
  dataLmeidForIntBAC[,timeVar] <- timesDrop2ForInt - epsDerXZ
  
  # Check for negative times
  aa <- which(dataLmeidForIntBAC[,timeVar] < 0)
  dataLmeidForIntBAC[aa,timeVar] <- timesDrop2ForInt[aa]
  
  xTimeDrop2ForIntFOW <- model.matrix(formula(fitlme),data = dataLmeidForIntFOW)
  zTimeDrop2ForIntFOW <- model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidForIntFOW)
  xTimeDrop2ForIntBAC <- model.matrix(formula(fitlme),data = dataLmeidForIntBAC)
  zTimeDrop2ForIntBAC <- model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidForIntBAC)
  
  xTimeDerDrop2ForInt <- (xTimeDrop2ForIntFOW - xTimeDrop2ForIntBAC)/(dataLmeidForIntFOW[,timeVar]-
                                                                        dataLmeidForIntBAC[,timeVar])
  zTimeDerDrop2ForInt <- (zTimeDrop2ForIntFOW - zTimeDrop2ForIntBAC)/(dataLmeidForIntFOW[,timeVar]-
                                                                        dataLmeidForIntBAC[,timeVar])
  colnames(xTimeDerDrop2ForInt) <- NULL;rownames(xTimeDerDrop2ForInt) <- NULL;
  colnames(zTimeDerDrop2ForInt) <- NULL;rownames(zTimeDerDrop2ForInt) <- NULL;
  
  rm(list = c("dataLmeidForIntFOW","dataLmeidForIntBAC","xTimeDrop2ForIntFOW","zTimeDrop2ForIntFOW",
              "xTimeDrop2ForIntBAC","zTimeDrop2ForIntBAC"))
  mtDerDrop2ForInt <- c(xTimeDerDrop2ForInt %*% fixed.effects(fitlme) + 
                          rowSums(zTimeDerDrop2ForInt*(brePred[idDrop2,][indDrop2,])))
  
  
  tstartDrop2ForInt <- timesDrop2ForInt
  tstopDrop2ForInt <- unlist(tapply(tstartDrop2ForInt,idDrop2Ind,function(x) c(x[2:length(x)],NA)  ))
  tstopDrop2ForInt[!duplicated(idDrop2Ind,fromLast = T)] <- timesDrop2
  deltaDrop2ForInt <- rep(0,length(timesDrop2ForInt))
  deltaDrop2ForInt[!duplicated(idDrop2Ind,fromLast = T)] <- deltaDrop2Sum
  
  # Starting values for the parameters of the dropout process
  # Covariates + predictions of the random effects based on the lme
  fitDrop2Start <- coxph( Surv(tstartDrop2ForInt,tstopDrop2ForInt,deltaDrop2ForInt) ~ xDrop2ForInt + 
                            mtPredDrop2ForInt + mtDerDrop2ForInt,
                          control = coxph.control(timefix = FALSE))
  summary(fitDrop2Start)
  
  
  # Obtain starting values for the B-spline parameters
  # by approximating the baseline cumulative based on the Cox model
  predChaz = basehaz(fitDrop2Start,centered = F)
  predChaz = predChaz[predChaz$hazard!=0,]
  predChaz$diff <- c(NA,diff(predChaz$hazard))
  predChaz <- predChaz[which(predChaz$diff!=0),]
  predChaz$logdiff <- log(predChaz$diff)
  predChaz <- predChaz[is.finite(predChaz$logdiff),]
  
  # B-spline basis matrix at the times where the hazard has been evaluated
  splDrop2 <- bSpline(predChaz$time,knots = IntKnotsDrop2,Boundary.knots = bKnotsDrop2,intercept = T)
  
  # Starting values for the B-spline parameters
  fit <- lm(predChaz$logdiff ~ -1 + splDrop2)
  
  lbDrop2 <- length(coef(fitDrop2Start))
  lassocDrop2 <- lbDrop2 - ncol(xDrop2)
  lpsiDrop2 <- ncol(BstDrop2Final)
  lDrop2 <- lbDrop2 + lpsiDrop2
  paramDrop2 <- c(coef(fitDrop2Start),coef(fit))
  names(paramDrop2) <- NULL
  paramDrop2_bup <- paramDrop2
  
  # Remove unnecessary objects
  rm(list = c("predChaz","fit","splDrop2","tstartDrop2ForInt","tstopDrop2ForInt",
              "deltaDrop2ForInt","dataLmeidForInt"))
  
  fDrop2 = function(paramDrop2)
  {
    # Parameters of the dropout model (regarding the observed marker values and other covariates)
    betasDrop2 <- paramDrop2[1:ncol(xDrop2)]
    
    # Association parameters (true marker value and true marker slope)
    assocDrop2 <- paramDrop2[(ncol(xDrop2)+1):lbDrop2]
    
    # Parameters for the B-spline basis matrix
    psiDrop2 <- paramDrop2[-(1:lbDrop2)]
    
    # Log hazard at the observed survival time
    loghazard <- c(xDrop2Final %*% betasDrop2 + BstDrop2Final %*% psiDrop2 + assocDrop2[1]*mtPredDrop2Final+
                     assocDrop2[2]*mtDerDrop2Final)
    
    # Log hazard at additional time points (for integration)
    loghazardForInt <- c((xDrop2ForInt %*% betasDrop2) + (BstDrop2ForInt %*% psiDrop2)) + (
      assocDrop2[1]*mtPredDrop2ForInt + assocDrop2[2]*mtDerDrop2ForInt)
    
    # Integration of the hazard
    hazInt <- p2Drop2*(exp(loghazardForInt)*wwsDrop2Long)
    
    # Log likelihood
    logl <- sum(loghazard*deltaDrop2Sum) - sum(hazInt)
    
    -logl
  }
  
  # Gradient of the log-likelihood for the dropout process
  fdGradDrop2 <- function(x,n=length(x)){fd(x,fDrop2,n)}
  cdGradDrop2 <- function(x){cd(x,fDrop2)}
  
  fitContrDrop2 <- try(optim(paramDrop2,fn = fDrop2,gr = cdGradDrop2,method = "BFGS",
                             control = list(trace = 0,maxit = 300)),silent = TRUE)
  if (!"try-error" %in% class(fitContrDrop2))
  {
    paramDrop2 <- fitContrDrop2$par
  } else{paramDrop2 <- paramDrop2_bup}
  
  # Matrix logarithm of the covariance matrix of the random effects
  decom <- eigen(getVarCov(fitlme)) 
  logDre <- decom$vectors %*% diag(log(decom$values)) %*% t(decom$vectors)
  
  # Starting values: obtained from lme
  paramStart <- c(fixed.effects(fitlme),log(fitlme$sigma^2),vech(logDre),
                  paramVis,paramDrop1,paramDrop2)
  names(paramStart) <- NULL
  
  # Function approximating the log-likelihood
  logMarg <- function(param)
  {
    # Parameters of the marker model
    betaMark <- param[1:pLong] # Fixed-effects of the LMM
    sigma2 <- exp(param[(pLong+1)]) # within-individual variance
    vechLogDre <- param[(pLong + 2):nParLong] # vech(log(Dre)) 
    
    # Parameters of the visiting model
    paramVis <- param[(nParLong + 1):(nParLong + lVis)]
    betasVis <- paramVis[1:ncol(xVis)]
    assocVis <- paramVis[(ncol(xVis)+1):lbVis]
    psiVis <- paramVis[-(1:lbVis)]
    
    # Parameters of the dropout model (cause 1)
    paramDrop1 <- param[(nParLong + lVis+1):(nParLong + lVis + lDrop1)]
    betasDrop1 <- paramDrop1[1:ncol(xDrop1)]
    assocDrop1 <- paramDrop1[(ncol(xDrop1)+1):lbDrop1]
    psiDrop1 <- paramDrop1[-(1:lbDrop1)]
    
    # Parameters of the dropout model (cause 2)
    paramDrop2 <- param[-(1:(nParLong + lVis + lDrop1))]
    betasDrop2 <- paramDrop2[1:ncol(xDrop2)]
    assocDrop2 <- paramDrop2[(ncol(xDrop2)+1):lbDrop2]
    psiDrop2 <- paramDrop2[-(1:lbDrop2)]
    
    # Matrix logarithm of Dre
    logDre <- matrix(Gqq %*% vechLogDre,nrow = qLong,ncol = qLong,byrow = F)
    decom <- eigen(logDre)
    
    # Variance-covariance of the random effects
    Dre <- decom$vectors %*% diag(exp(decom$values)) %*% t(decom$vectors)
    DreInv <- solve(Dre)
    
    # Get cholesky factor of the variance-covariance matrix of the random effects
    LDre <- chol(Dre)
    LDreT <- t(LDre)
    
    ######################
    ### Visiting model ###
    ######################
    # Log-hazard at visit times (tstop)
    loghazardVis <- c(xVis %*% betasVis + BstVis %*% psiVis)
    
    # Log-hazard at additional time points for integration (ignoring the association parameters)
    loghazardVisForInt <- c((xVisForInt %*% betasVis) + (BstVisForInt %*% psiVis))
    
    # True marker value at the observed survival time (ignoring random effects)
    mtPredVis <- c(xTimeVis %*% betaMark)
    
    # True marker value at additional time points for integration (ignoring random effects)
    mtPredVisForInt <- c(xTimeVisForInt %*% betaMark)
    
    # Derivative of the true marker value at the observed survival time (ignoring random effects)
    mtDerVis <- c(xTimeDerVis %*% betaMark)
    
    # Derivative of the true marker value at additional time points for integration (ignoring random effects)
    mtDerVisForInt <- c(xTimeDerVisForInt %*% betaMark)
    
    ###############################
    ### Dropout model (cause 1) ###
    ###############################
    
    # Log-hazard at the observed survival time (ignoring the association parameters)
    loghazardDrop1 <- c(xDrop1Final %*% betasDrop1 + BstDrop1Final %*% psiDrop1)
    
    # Log-hazard at additional time points for integration (ignoring the association parameters)
    loghazardDrop1ForInt <- c((xDrop1ForInt %*% betasDrop1) + (BstDrop1ForInt %*% psiDrop1))
    
    # True marker value at the observed survival time (ignoring random effects)
    mtPredDrop1Final <- c(xTimeDrop1 %*% betaMark)
    
    # True marker value at additional time points for integration (ignoring random effects)
    mtPredDrop1ForInt <- c(xTimeDrop1ForInt %*% betaMark)
    
    # Derivative of the true marker value at the observed survival time (ignoring random effects)
    mtDerDrop1Final <- c(xTimeDerDrop1 %*% betaMark)
    
    # Derivative of the true marker value at additional time points for integration (ignoring random effects)
    mtDerDrop1ForInt <- c(xTimeDerDrop1ForInt %*% betaMark)
    
    ###############################
    ### Dropout model (cause 2) ###
    ###############################
    
    # Log-hazard at the observed survival time (ignoring the association parameters)
    loghazardDrop2 <- c(xDrop2Final %*% betasDrop2 + BstDrop2Final %*% psiDrop2)
    
    # Log-hazard at additional time points for integration (ignoring the association parameters)
    loghazardDrop2ForInt <- c((xDrop2ForInt %*% betasDrop2) + (BstDrop2ForInt %*% psiDrop2))
    
    # True marker value at the observed survival time (ignoring random effects)
    mtPredDrop2Final <- c(xTimeDrop2 %*% betaMark)
    
    # True marker value at additional time points for integration (ignoring random effects)
    mtPredDrop2ForInt <- c(xTimeDrop2ForInt %*% betaMark)
    
    # Derivative of the true marker value at the observed survival time (ignoring random effects)
    mtDerDrop2Final <- c(xTimeDerDrop2 %*% betaMark)
    
    # Derivative of the true marker value at additional time points for integration (ignoring random effects)
    mtDerDrop2ForInt <- c(xTimeDerDrop2ForInt %*% betaMark)
    
    for (i in 1:N)
    {
      ####################
      ### Marker model ###
      ####################
      zi = zLongList[[i]]
      zzi = zzLongList[[i]]
      
      # Residuals
      ri <- yList[[i]] - xLongList[[i]] %*% betaMark
      nii <- ni[i]
      
      # Vi: variance-covariance of Yi
      Vi <- diag(ni[i])*sigma2 + tcrossprod(zi %*% LDreT )
      Li = .Internal(La_chol(Vi, pivot, tol)) # Cholesky factor of Vi
      
      # Quadratic form, ri^{\top}Vi^{-1}ri
      tmp = .Internal(backsolve(Li, ri,nii,TRUE,TRUE))
      rss = sum(tmp^2)
      
      logl[i] <- - sum(log(diag(Li))) - 0.5*nii*log(2*pi) - 0.5*rss
      
      ################################################################
      ### Contribution of the random effects in the visiting model ###
      ################################################################
      # Posterior distribution of b_{i}|y_{i}
      
      # Inverse covariance matrix of b_{i}|y_{i}
      Linf = DreInv + zzi/sigma2
      
      # Cholesky of Linf 
      Bi = .Internal(La_chol(Linf, pivot, tol))
      temp = crossprod(zi,ri/sigma2)
      
      # Mean of b_{i}|y_{i}
      mub = c(.Internal( backsolve(Bi, .Internal(backsolve(Bi,temp,qLong,TRUE,TRUE)) ,qLong,TRUE,FALSE )  ))
      
      # Approximation of the integral using quasi Monte Carlo
      # br = mnormL(1500,mu = mub,L = Bi)
      # assocVisBr <- c(br %*% assocVis)
      # logVis[i] = log(mean(exp(logVis1[i] + (nii-1)*assocVisBr - logVis2[i]*exp(assocVisBr))))
      
      ####################################################################
      ### Approximation of the integral using Gauss-Hermite quadrature ###
      ####################################################################
      # Quadrature points
      rgrid = t(sqrt(2)*.Internal(backsolve(Bi,tbgrid,qLong,T,F)) + c(mub))
      trgrid <- t(rgrid)
      
      # Indices for measurements of individual i in the visiting dataset
      jjVis <- idsVis[[i]]
      
      # True marker values at survival time (for all GH points for the random effects)
      mtPredVisi <- mtPredVis[jjVis] + zTimeVis[jjVis,] %*% trgrid
      
      # Derivative of the true marker values at survival time (for all GH points for the random effects)
      mtDerVisi <- mtDerVis[jjVis] + zTimeDerVis[jjVis,] %*% trgrid  
      
      # Log-hazard at survival time (for all GH points for the random effects)
      logVisContr <- colSums(deltaVis[jjVis]*(loghazardVis[jjVis] + mtPredVisi*assocVis[1] +
                                                mtDerVisi*assocVis[2]))
      
      # Relevant time points for integration
      jjVisFI <- idsVisForInt[[i]]
      
      # True marker values at additional time points for integration
      # Rows correspond to time points (for integration) and columns to GH points for the random effects
      mtPredVisForInti <- mtPredVisForInt[jjVisFI] + zTimeVisForInt[jjVisFI,] %*% trgrid
      
      # Derivative of the true marker values at additional time points for integration
      # Rows correspond to time points (for integration) and columns to GH points for the random effects
      mtDerVisForInti <- mtDerVisForInt[jjVisFI] + zTimeDerVisForInt[jjVisFI,] %*% trgrid
      
      # Function to be integrated
      aaVis <- loghazardVisForInt[jjVisFI] + assocVis[1]*mtPredVisForInti + assocVis[2]*mtDerVisForInti
      
      xxVis <- p2Vis[jjVisFI]*(exp(aaVis)*wwsVisLong[jjVisFI])
      logVisContr <- logVisContr - colSums(xxVis)
      
      # ff <- function(xx)
      # {
      #   temp <- bSpline(xx,knots = IntKnotsVis,Boundary.knots = bKnotsVis,intercept = T)
      #   mt = betaMark[1] + rgrid[97,1] + (betaMark[2] + rgrid[97,2])*(xx+timesVis[idVis==i][4])
      #   mtDer = betaMark[2] + rgrid[97,2]
      #   exp( c(temp %*% psiVis) + c(xVis[idVis==i,][4,] %*% betasVis) +
      #        assocVis[1]*mt + assocVis[2]*mtDer)
      # }
      # integrate(ff,lower = tstartVis[idVis==i][4],upper = tstopVis[idVis==i][4])$val
      
      ###############################
      ### Dropout model (cause 1) ###
      ###############################
      
      # True marker values at survival time (for all GH points for the random effects)
      mtPredDrop1Finali <- mtPredDrop1Final[i] + c(rgrid %*% zTimeDrop1[i,])
      
      # Derivative of the true marker values at survival time (for all GH points for the random effects)
      mtDerDrop1Finali <- mtDerDrop1Final[i] + c(rgrid %*% zTimeDerDrop1[i,])
      
      # Log-hazard at survival time (for all GH points for the random effects)
      logDrop1Contr <- deltaDrop1Sum[i]*(loghazardDrop1[i] + mtPredDrop1Finali*assocDrop1[1] +
                                       mtDerDrop1Finali*assocDrop1[2])

      # Relevant time points for integration
      jjDrop1 <- idsDrop1ForInt[[i]]
      
      # True marker values at additional time points for integration
      # Rows correspond to time points (for integration) and columns to GH points for the random effects
      mtPredDrop1ForInti <- mtPredDrop1ForInt[jjDrop1] + zTimeDrop1ForInt[jjDrop1,] %*% trgrid
      
      # Derivative of the true marker values at additional time points for integration
      # Rows correspond to time points (for integration) and columns to GH points for the random effects
      mtDerDrop1ForInti <- mtDerDrop1ForInt[jjDrop1] + zTimeDerDrop1ForInt[jjDrop1,]%*% trgrid
      
      # Function to be integrated
      aaDrop1 <- loghazardDrop1ForInt[jjDrop1] + assocDrop1[1]*mtPredDrop1ForInti + assocDrop1[2]*mtDerDrop1ForInti
      
      xxDrop1 <- p2Drop1[jjDrop1]*(exp(aaDrop1)*wwsDrop1Long[jjDrop1])
      logDrop1Contr <- logDrop1Contr - colSums(xxDrop1)
      
      # ff <- function(xx)
      # {
      #   temp <- bSpline(xx,knots = IntKnotsDrop1,Boundary.knots = bKnotsDrop1,intercept = T)
      #   mt = betaMark[1] + rgrid[97,1] + (betaMark[2] + rgrid[97,2])*xx
      #   mtDer = betaMark[2] + rgrid[97,2]
      #   exp( c(temp %*% psiDrop1) + c(xDrop1[idDrop1==i,][7,] %*% betasDrop1) +
      #        assocDrop1[1]*mt + assocDrop1[2]*mtDer)
      # }
      # integrate(ff,lower = tstartDrop1[idDrop1==i][7],upper = tstopDrop1[idDrop1==i][7])$val
      
      ###############################
      ### Dropout model (cause 2) ###
      ###############################
      
      # True marker values at survival time (for all GH points for the random effects)
      mtPredDrop2Finali <- mtPredDrop2Final[i] + c(rgrid %*% zTimeDrop2[i,])
      
      # Derivative of the true marker values at survival time (for all GH points for the random effects)
      mtDerDrop2Finali <- mtDerDrop2Final[i] + c(rgrid %*% zTimeDerDrop2[i,])
      
      # Log-hazard at survival time (for all GH points for the random effects)
      logDrop2Contr <- deltaDrop2Sum[i]*(loghazardDrop2[i] + mtPredDrop2Finali*assocDrop2[1] +
                                           mtDerDrop2Finali*assocDrop2[2])
      
      # Relevant time points for integration
      jjDrop2 <- idsDrop2ForInt[[i]]
      
      # True marker values at additional time points for integration
      # Rows correspond to time points (for integration) and columns to GH points for the random effects
      mtPredDrop2ForInti <- mtPredDrop2ForInt[jjDrop2] + zTimeDrop2ForInt[jjDrop2,] %*% trgrid
      
      # Derivative of the true marker values at additional time points for integration
      # Rows correspond to time points (for integration) and columns to GH points for the random effects
      mtDerDrop2ForInti <- mtDerDrop2ForInt[jjDrop2] + zTimeDerDrop2ForInt[jjDrop2,]%*% trgrid
      
      # Function to be integrated
      aaDrop2 <- loghazardDrop2ForInt[jjDrop2] + assocDrop2[1]*mtPredDrop2ForInti + assocDrop2[2]*mtDerDrop2ForInti
      
      xxDrop2 <- p2Drop2[jjDrop2]*(exp(aaDrop2)*wwsDrop2Long[jjDrop2])
      logDrop2Contr <- logDrop2Contr - colSums(xxDrop2)
      
      # ff <- function(xx)
      # {
      #   temp <- bSpline(xx,knots = IntKnotsDrop2,Boundary.knots = bKnotsDrop2,intercept = T)
      #   mt = betaMark[1] + rgrid[97,1] + (betaMark[2] + rgrid[97,2])*xx
      #   mtDer = betaMark[2] + rgrid[97,2]
      #   exp( c(temp %*% psiDrop2) + c(xDrop2[idDrop2==i,][7,] %*% betasDrop2) +
      #        assocDrop2[1]*mt + assocDrop2[2]*mtDer)
      # }
      # integrate(ff,lower = tstartDrop2[idDrop2==i][7],upper = tstopDrop2[idDrop2==i][7])$val
      
      logVisDropCond <- logVisContr+logDrop1Contr+logDrop2Contr
      mm <- median(logVisDropCond)
      logVisDrop[i] = log(pi^(-qLong/2)*sum(weights*exp(logVisDropCond-mm)  )) + mm
      
    }
    return(-sum(logl) - sum(logVisDrop))
    
  }
  
  # Function approximating the log-likelihood
  logMarker <- function(param)
  {
    # Parameters of the marker model
    betaMark <- param[1:pLong] # Fixed-effects of the LMM
    sigma2 <- exp(param[(pLong+1)]) # within-individual variance
    vechLogDre <- param[(pLong + 2):nParLong] # vech(log(Dre)) 
    
    # Matrix logarithm of Dre
    logDre <- matrix(Gqq %*% vechLogDre,nrow = qLong,ncol = qLong,byrow = F)
    decom <- eigen(logDre)
    
    # Variance-covariance of the random effects
    Dre <- decom$vectors %*% diag(exp(decom$values)) %*% t(decom$vectors)
    DreInv <- solve(Dre)
    
    # Get cholesky factor of the variance-covariance matrix of the random effects
    LDre <- chol(Dre)
    LDreT <- t(LDre)
    
    for (i in 1:N)
    {
      ####################
      ### Marker model ###
      ####################
      zi = zLongList[[i]]
      zzi = zzLongList[[i]]
      
      # Residuals
      ri <- yList[[i]] - xLongList[[i]] %*% betaMark
      nii <- ni[i]
      
      # Vi: variance-covariance of Yi
      Vi <- diag(ni[i])*sigma2 + tcrossprod(zi %*% LDreT )
      Li = .Internal(La_chol(Vi, pivot, tol)) # Cholesky factor of Vi
      
      # Quadratic form, ri^{\top}Vi^{-1}ri
      tmp = .Internal(backsolve(Li, ri,nii,TRUE,TRUE))
      rss = sum(tmp^2)
      
      logl[i] <- - sum(log(diag(Li))) - 0.5*nii*log(2*pi) - 0.5*rss
      
    }
    return( -sum(logl) )
    
  }
  
  # Score vector of model parameters for each individual
  SCbetasVis <- matrix(NA,nr = N,nc = ncol(xVis))
  SCassocVis <- matrix(NA,nr = N,nc = lassocVis)
  SCpsiVis <- matrix(NA,nr = N,nc = lpsiVis)
  SCbetasDrop1 <- matrix(NA,nr = N,nc = ncol(xDrop1Final))
  SCassocDrop1 <- matrix(NA,nr = N,nc = lassocDrop1)
  SCpsiDrop1 <- matrix(NA,nr = N,nc = lpsiDrop1)
  
  SCbetasDrop2 <- matrix(NA,nr = N,nc = ncol(xDrop2Final))
  SCassocDrop2 <- matrix(NA,nr = N,nc = lassocDrop2)
  SCpsiDrop2 <- matrix(NA,nr = N,nc = lpsiDrop2)
  
  SClogSigma <- rep(NA,N)
  SCbetaMark <- matrix(NA,nr = N,nc = pLong)
  SCvecD <- matrix(NA,nr = N,nc = qLong^2)
  crossRgrid <- matrix(NA,nr = nrow(bgrid),nc = qLong^2)
  
  vecEYvechY = function(vechLogDre)
  {
    logD = matrix(Gqq %*% vechLogDre,nr = qLong,nc = qLong)
    eigLogD = eigen(logD)
    
    # Get the eigenvectors in a list
    ujs = split(t(eigLogD$vectors),1:qLong)
    
    # A vector with the eigenvalues
    ljs = eigLogD$values
    
    # Ejs[[j]] n*n matrix with all its entries zero except for identity in (j,j)
    Ejs = list()
    
    # ejs is a list with the rows of the identity matrix
    ejs = list()
    ttujs = list()
    
    zerom = matrix(0,nr = qLong,nc = qLong)
    zerov = rep(0,qLong)
    
    for (j in 1:qLong)
    {
      Ejs[[j]] = zerom
      Ejs[[j]][j,j] = 1
      
      ejs[[j]] = zerov
      ejs[[j]][j] = 1
      
      ttujs[[j]] = tcrossprod(ujs[[j]])
    }
    #################################################################################
    ### Suppose Y = U * ? * U', where ? is a diagonal matrix with the eigenvalues ###
    ### and U the matrix including the corresponding eigenvectors                 ###
    #################################################################################
    
    # ss evaluates the derivative dvec(exp(?))/dvec(Y)'
    ss = matrix(0,nr = qLong^2,nc = qLong^2)
    
    # ss2 evaluates the derivative dvec(U)/dvec(Y)'
    ss2 = ss
    
    for (j in 1:qLong)
    {
      temp = tcrossprod( c(Ejs[[j]]), c(ttujs[[j]]) )*exp(ljs[j])
      
      ss = ss + temp
      
      ss2 = ss2 +  ejs[[j]] %x% diag(qLong) %*% ( t(ujs[[j]]) %x% ginv(ljs[j]*diag(qLong) - logD) )
    }
    eUD = eigLogD$vectors %*% diag( exp(eigLogD$values) )
    
    # The jacobian of tranforming vec(Y) to vec(exp(Y)) through the product rule
    res = ( eUD %x% diag(qLong) ) %*% ss2
    
    res = res +  (eigLogD$vectors %x% eigLogD$vectors) %*% ss
    
    res = res + ( diag(qLong) %x% eUD ) %*% Kqq %*% ss2
    
    # Use the duplication matrix to get the jacobian of tranforming vech(Y) to vech(exp(Y))
    res %*% Gqq
  }
  
  # Function approximating the score vector of the model
  ScoreVec <- function(param)
  {
    # Parameters of the marker model
    betaMark <- param[1:pLong] # Fixed-effects of the LMM
    sigma2 <- exp(param[(pLong+1)]) # within-individual variance
    vechLogDre <- param[(pLong + 2):nParLong] # vech(log(Dre)) 
    
    # Parameters of the visiting model
    paramVis <- param[(nParLong + 1):(nParLong + lVis)]
    betasVis <- paramVis[1:ncol(xVis)]
    assocVis <- paramVis[(ncol(xVis)+1):lbVis]
    psiVis <- paramVis[-(1:lbVis)]
    
    # Parameters of the dropout model (cause 1)
    paramDrop1 <- param[(nParLong + lVis+1):(nParLong + lVis + lDrop1)]
    betasDrop1 <- paramDrop1[1:ncol(xDrop1)]
    assocDrop1 <- paramDrop1[(ncol(xDrop1)+1):lbDrop1]
    psiDrop1 <- paramDrop1[-(1:lbDrop1)]
    
    # Parameters of the dropout model (cause 2)
    paramDrop2 <- param[-(1:(nParLong + lVis + lDrop1))]
    betasDrop2 <- paramDrop2[1:ncol(xDrop2)]
    assocDrop2 <- paramDrop2[(ncol(xDrop2)+1):lbDrop2]
    psiDrop2 <- paramDrop2[-(1:lbDrop2)]
    
    # Matrix logarithm of Dre
    logDre <- matrix(Gqq %*% vechLogDre,nrow = qLong,ncol = qLong,byrow = F)
    decom <- eigen(logDre)
    
    # Variance-covariance of the random effects
    Dre <- decom$vectors %*% diag(exp(decom$values)) %*% t(decom$vectors)
    DreInv <- solve(Dre)
    DreInvKron <- DreInv %x% DreInv
    
    # Get cholesky factor of the variance-covariance matrix of the random effects
    LDre <- chol(Dre)
    LDreT <- t(LDre)
    
    ######################
    ### Visiting model ###
    ######################
    # Log-hazard at visit times (tstop)
    loghazardVis <- c(xVis %*% betasVis + BstVis %*% psiVis)
    
    # Log-hazard at additional time points for integration (ignoring the association parameters)
    loghazardVisForInt <- c((xVisForInt %*% betasVis) + (BstVisForInt %*% psiVis))
    
    # True marker value at the observed survival time (ignoring random effects)
    mtPredVis <- c(xTimeVis %*% betaMark)
    
    # True marker value at additional time points for integration (ignoring random effects)
    mtPredVisForInt <- c(xTimeVisForInt %*% betaMark)
    
    # Derivative of the true marker value at the observed survival time (ignoring random effects)
    mtDerVis <- c(xTimeDerVis %*% betaMark)
    
    # Derivative of the true marker value at additional time points for integration (ignoring random effects)
    mtDerVisForInt <- c(xTimeDerVisForInt %*% betaMark)
    
    ###############################
    ### Dropout model (cause 1) ###
    ###############################
    
    # Log-hazard at the observed survival time (ignoring the association parameters)
    loghazardDrop1 <- c(xDrop1Final %*% betasDrop1 + BstDrop1Final %*% psiDrop1)
    
    # Log-hazard at additional time points for integration (ignoring the association parameters)
    loghazardDrop1ForInt <- c((xDrop1ForInt %*% betasDrop1) + (BstDrop1ForInt %*% psiDrop1))
    
    # True marker value at the observed survival time (ignoring random effects)
    mtPredDrop1Final <- c(xTimeDrop1 %*% betaMark)
    
    # True marker value at additional time points for integration (ignoring random effects)
    mtPredDrop1ForInt <- c(xTimeDrop1ForInt %*% betaMark)
    
    # Derivative of the true marker value at the observed survival time (ignoring random effects)
    mtDerDrop1Final <- c(xTimeDerDrop1 %*% betaMark)
    
    # Derivative of the true marker value at additional time points for integration (ignoring random effects)
    mtDerDrop1ForInt <- c(xTimeDerDrop1ForInt %*% betaMark)
    
    ###############################
    ### Dropout model (cause 2) ###
    ###############################
    
    # Log-hazard at the observed survival time (ignoring the association parameters)
    loghazardDrop2 <- c(xDrop2Final %*% betasDrop2 + BstDrop2Final %*% psiDrop2)
    
    # Log-hazard at additional time points for integration (ignoring the association parameters)
    loghazardDrop2ForInt <- c((xDrop2ForInt %*% betasDrop2) + (BstDrop2ForInt %*% psiDrop2))
    
    # True marker value at the observed survival time (ignoring random effects)
    mtPredDrop2Final <- c(xTimeDrop2 %*% betaMark)
    
    # True marker value at additional time points for integration (ignoring random effects)
    mtPredDrop2ForInt <- c(xTimeDrop2ForInt %*% betaMark)
    
    # Derivative of the true marker value at the observed survival time (ignoring random effects)
    mtDerDrop2Final <- c(xTimeDerDrop2 %*% betaMark)
    
    # Derivative of the true marker value at additional time points for integration (ignoring random effects)
    mtDerDrop2ForInt <- c(xTimeDerDrop2ForInt %*% betaMark)
    
    for (i in 1:N)
    {
      ####################
      ### Marker model ###
      ####################
      zi = zLongList[[i]]
      zzi = zzLongList[[i]]
      
      # Residuals
      xi <- xLongList[[i]]
      ri <- yList[[i]] - xi %*% betaMark
      nii <- ni[i]
      
      # Vi: variance-covariance of Yi
      Vi <- diag(ni[i])*sigma2 + tcrossprod(zi %*% LDreT )
      Li = .Internal(La_chol(Vi, pivot, tol)) # Cholesky factor of Vi
      
      # Quadratic form, ri^{\top}Vi^{-1}ri
      tmp = .Internal(backsolve(Li, ri,nii,TRUE,TRUE))
      rss = sum(tmp^2)
      
      logl[i] <- - sum(log(diag(Li))) - 0.5*nii*log(2*pi) - 0.5*rss
      
      ################################################################
      ### Contribution of the random effects in the visiting model ###
      ################################################################
      # Posterior distribution of b_{i}|y_{i}
      
      # Inverse covariance matrix of b_{i}|y_{i}
      Linf = DreInv + zzi/sigma2
      
      # Cholesky of Linf 
      Bi = .Internal(La_chol(Linf, pivot, tol))
      temp = crossprod(zi,ri/sigma2)
      
      # Mean of b_{i}|y_{i}
      mub = c(.Internal( backsolve(Bi, .Internal(backsolve(Bi,temp,qLong,TRUE,TRUE)) ,qLong,TRUE,FALSE )  ))
      
      # Approximation of the integral using quasi Monte Carlo
      # br = mnormL(1500,mu = mub,L = Bi)
      # assocVisBr <- c(br %*% assocVis)
      # logVis[i] = log(mean(exp(logVis1[i] + (nii-1)*assocVisBr - logVis2[i]*exp(assocVisBr))))
      
      ####################################################################
      ### Approximation of the integral using Gauss-Hermite quadrature ###
      ####################################################################
      # Quadrature points
      rgrid = t(sqrt(2)*.Internal(backsolve(Bi,tbgrid,qLong,T,F)) + c(mub))
      trgrid <- t(rgrid)
      
      # Indices for measurements of individual i in the visiting dataset
      jjVis <- idsVis[[i]]
      
      # True marker values at survival time (for all GH points for the random effects)
      mtPredVisi <- mtPredVis[jjVis] + zTimeVis[jjVis,] %*% trgrid
      
      # Derivative of the true marker values at survival time (for all GH points for the random effects)
      mtDerVisi <- mtDerVis[jjVis] + zTimeDerVis[jjVis,] %*% trgrid  
      
      # Log-hazard at survival time (for all GH points for the random effects)
      logVisContr <- colSums(deltaVis[jjVis]*(loghazardVis[jjVis] + mtPredVisi*assocVis[1] +
                                                mtDerVisi*assocVis[2]))
      
      # Relevant time points for integration
      jjVisFI <- idsVisForInt[[i]]
      
      # True marker values at additional time points for integration
      # Rows correspond to time points (for integration) and columns to GH points for the random effects
      mtPredVisForInti <- mtPredVisForInt[jjVisFI] + zTimeVisForInt[jjVisFI,] %*% trgrid
      
      # Derivative of the true marker values at additional time points for integration
      # Rows correspond to time points (for integration) and columns to GH points for the random effects
      mtDerVisForInti <- mtDerVisForInt[jjVisFI] + zTimeDerVisForInt[jjVisFI,] %*% trgrid
      
      # Function to be integrated
      aaVis <- loghazardVisForInt[jjVisFI] + assocVis[1]*mtPredVisForInti + assocVis[2]*mtDerVisForInti
      
      xxVis <- p2Vis[jjVisFI]*(exp(aaVis)*wwsVisLong[jjVisFI])
      logVisContr <- logVisContr - colSums(xxVis)
      
      ###############################
      ### Dropout model (cause 1) ###
      ###############################
      
      # True marker values at survival time (for all GH points for the random effects)
      mtPredDrop1Finali <- mtPredDrop1Final[i] + c(rgrid %*% zTimeDrop1[i,])
      
      # Derivative of the true marker values at survival time (for all GH points for the random effects)
      mtDerDrop1Finali <- mtDerDrop1Final[i] + c(rgrid %*% zTimeDerDrop1[i,])
      
      # Log-hazard at survival time (for all GH points for the random effects)
      logDrop1Contr <- deltaDrop1Sum[i]*(loghazardDrop1[i] + mtPredDrop1Finali*assocDrop1[1] +
                                         mtDerDrop1Finali*assocDrop1[2])
      
      # Relevant time points for integration
      jjDrop1 <- idsDrop1ForInt[[i]]
      
      # True marker values at additional time points for integration
      # Rows correspond to time points (for integration) and columns to GH points for the random effects
      mtPredDrop1ForInti <- mtPredDrop1ForInt[jjDrop1] + zTimeDrop1ForInt[jjDrop1,] %*% trgrid
      
      # Derivative of the true marker values at additional time points for integration
      # Rows correspond to time points (for integration) and columns to GH points for the random effects
      mtDerDrop1ForInti <- mtDerDrop1ForInt[jjDrop1] + zTimeDerDrop1ForInt[jjDrop1,]%*% trgrid
      
      # Function to be integrated
      aaDrop1 <- loghazardDrop1ForInt[jjDrop1] + assocDrop1[1]*mtPredDrop1ForInti + assocDrop1[2]*mtDerDrop1ForInti
      
      xxDrop1 <- p2Drop1[jjDrop1]*(exp(aaDrop1)*wwsDrop1Long[jjDrop1])
      logDrop1Contr <- logDrop1Contr - colSums(xxDrop1)
    
      ###############################
      ### Dropout model (cause 2) ###
      ###############################
      
      # True marker values at survival time (for all GH points for the random effects)
      mtPredDrop2Finali <- mtPredDrop2Final[i] + c(rgrid %*% zTimeDrop2[i,])
      
      # Derivative of the true marker values at survival time (for all GH points for the random effects)
      mtDerDrop2Finali <- mtDerDrop2Final[i] + c(rgrid %*% zTimeDerDrop2[i,])
      
      # Log-hazard at survival time (for all GH points for the random effects)
      logDrop2Contr <- deltaDrop2Sum[i]*(loghazardDrop2[i] + mtPredDrop2Finali*assocDrop2[1] +
                                           mtDerDrop2Finali*assocDrop2[2])
      
      # Relevant time points for integration
      jjDrop2 <- idsDrop2ForInt[[i]]
      
      # True marker values at additional time points for integration
      # Rows correspond to time points (for integration) and columns to GH points for the random effects
      mtPredDrop2ForInti <- mtPredDrop2ForInt[jjDrop2] + zTimeDrop2ForInt[jjDrop2,] %*% trgrid
      
      # Derivative of the true marker values at additional time points for integration
      # Rows correspond to time points (for integration) and columns to GH points for the random effects
      mtDerDrop2ForInti <- mtDerDrop2ForInt[jjDrop2] + zTimeDerDrop2ForInt[jjDrop2,]%*% trgrid
      
      # Function to be integrated
      aaDrop2 <- loghazardDrop2ForInt[jjDrop2] + assocDrop2[1]*mtPredDrop2ForInti + assocDrop2[2]*mtDerDrop2ForInti
      
      xxDrop2 <- p2Drop2[jjDrop2]*(exp(aaDrop2)*wwsDrop2Long[jjDrop2])
      logDrop2Contr <- logDrop2Contr - colSums(xxDrop2)
      
      # Marginal likelihood Ui,Ti,di|Yi
      logVisDropCond <- logVisContr+logDrop1Contr+logDrop2Contr
      mm <- median(logVisDropCond)
      logVisDrop[i] = log(pi^(-qLong/2)*sum(weights*exp(logVisDropCond-mm)  )) + mm
      
      # Ratio of conditional likelihood Ti,di,ui|bi,yi to marginal likelihood Ti,di,ui|yi
      RCondMarg <- exp(logVisDropCond-logVisDrop[i])
      
      ########################
      ### Visiting process ###
      ########################
      
      # Score vector for betasVis
      xVisi <- xVis[jjVis,];if (ncol(xVis)==1){xVisi <- as.matrix(xVisi)};if (length(jjVis)==1){xVisi <- t(xVisi)}
      deltaVisi <- deltaVis[jjVis]
      xVisForInti <- xVisForInt[jjVisFI,];if (ncol(xVis)==1){xVisForInti <- as.matrix(xVisForInti)}
      for (j in 1:length(betasVis))
      {
        cSCbetasVis <- sum(deltaVisi*xVisi[,j]) - colSums(xVisForInti[,j]*xxVis)
        SCbetasVis[i,j] <- pi^(-qLong/2)*sum(weights*(cSCbetasVis*RCondMarg))
      }
      
      # Score vector for assocVis
      cSCassocVis1 <- colSums(deltaVisi*mtPredVisi) - colSums(mtPredVisForInti*xxVis)
      SCassocVis[i,1] <- pi^(-qLong/2)*sum(weights*(cSCassocVis1*RCondMarg))
      
      cSCassocVis2 <- colSums(deltaVisi*mtDerVisi) - colSums(mtDerVisForInti*xxVis)
      SCassocVis[i,2] <- pi^(-qLong/2)*sum(weights*(cSCassocVis2*RCondMarg))
      
      # Score vector for psiVis
      BstVisForInti <- BstVisForInt[jjVisFI,]
      BstVisi <- BstVis[jjVis,];if (length(jjVis)==1){BstVisi <- t(BstVisi)}
      for (j in 1:lpsiVis)
      {
        cSCpsiVis <- sum(deltaVisi*BstVisi[,j]) - colSums(BstVisForInti[,j]*xxVis)
        SCpsiVis[i,j] <- pi^(-qLong/2)*sum(weights*(cSCpsiVis*RCondMarg))
      }
      
      #######################
      ### Dropout process ###
      #######################
      
      # Score vector for betasDrop1
      xDrop1ForInti <- xDrop1ForInt[jjDrop1,];if (ncol(xDrop1)==1){xDrop1ForInti <- as.matrix(xDrop1ForInti)}
      for (j in 1:length(betasDrop1))
      {
        cSCbetasDrop1 <- deltaDrop1Sum[i]*xDrop1Final[i,j] - colSums(xDrop1ForInti[,j]*xxDrop1)
        SCbetasDrop1[i,j] <- pi^(-qLong/2)*sum(weights*(cSCbetasDrop1*RCondMarg))
      }
      
      # Score vector for betasDrop2
      xDrop2ForInti <- xDrop2ForInt[jjDrop2,];if (ncol(xDrop2)==1){xDrop2ForInti <- as.matrix(xDrop2ForInti)}
      for (j in 1:length(betasDrop2))
      {
        cSCbetasDrop2 <- deltaDrop2Sum[i]*xDrop2Final[i,j] - colSums(xDrop2ForInti[,j]*xxDrop2)
        SCbetasDrop2[i,j] <- pi^(-qLong/2)*sum(weights*(cSCbetasDrop2*RCondMarg))
      }
      
      # Score vector for assocDrop1
      cSCassocDrop11 <- deltaDrop1Sum[i]*mtPredDrop1Finali - colSums(mtPredDrop1ForInti*xxDrop1)
      SCassocDrop1[i,1] <- pi^(-qLong/2)*sum(weights*(cSCassocDrop11*RCondMarg))
      
      cSCassocDrop12 <- deltaDrop1Sum[i]*mtDerDrop1Finali - colSums(mtDerDrop1ForInti*xxDrop1)
      SCassocDrop1[i,2] <- pi^(-qLong/2)*sum(weights*(cSCassocDrop12*RCondMarg))
      
      # Score vector for assocDrop2
      cSCassocDrop21 <- deltaDrop2Sum[i]*mtPredDrop2Finali - colSums(mtPredDrop2ForInti*xxDrop2)
      SCassocDrop2[i,1] <- pi^(-qLong/2)*sum(weights*(cSCassocDrop21*RCondMarg))
      
      cSCassocDrop22 <- deltaDrop2Sum[i]*mtDerDrop2Finali - colSums(mtDerDrop2ForInti*xxDrop2)
      SCassocDrop2[i,2] <- pi^(-qLong/2)*sum(weights*(cSCassocDrop22*RCondMarg))
      
      # Score vector for psiDrop1
      BstDrop1ForInti <- BstDrop1ForInt[jjDrop1,]
      BstDrop1Finali <- BstDrop1Final[i,]
      for (j in 1:lpsiDrop1)
      {
        cSCpsiDrop1 <- deltaDrop1Sum[i]*BstDrop1Finali[j] - colSums(BstDrop1ForInti[,j]*xxDrop1)
        SCpsiDrop1[i,j] <- pi^(-qLong/2)*sum(weights*(cSCpsiDrop1*RCondMarg))
      }
      
      # Score vector for psiDrop2
      BstDrop2ForInti <- BstDrop2ForInt[jjDrop2,]
      BstDrop2Finali <- BstDrop2Final[i,]
      for (j in 1:lpsiDrop2)
      {
        cSCpsiDrop2 <- deltaDrop2Sum[i]*BstDrop2Finali[j] - colSums(BstDrop2ForInti[,j]*xxDrop2)
        SCpsiDrop2[i,j] <- pi^(-qLong/2)*sum(weights*(cSCpsiDrop2*RCondMarg))
      }
      
      # Score vector for log(sigma2)
      ei <- c(ri) - zi %*% trgrid
      cSClogSigma <- -0.5*nii + 0.5*colSums(ei^2)/sigma2
      SClogSigma[i] <- pi^(-qLong/2)*sum(weights*(cSClogSigma*RCondMarg))
      
      # Score vector for betaMark
      xTimeVisi <- xTimeVis[jjVis,]
      xTimeDerVisi <- xTimeDerVis[jjVis,];if (length(jjVis)==1){xTimeVisi <- t(xTimeVisi);xTimeDerVisi <- t(xTimeDerVisi)}
      tempBetaMark <- crossprod(xi,ei)/sigma2 + colSums(deltaVisi*(assocVis[1]*xTimeVisi+
                                                                   assocVis[2]*xTimeDerVisi))
      tempBetaMark <- tempBetaMark + deltaDrop1Sum[i]*(assocDrop1[1]*xTimeDrop1[i,]+assocDrop1[2]*xTimeDerDrop1[i,])
      tempBetaMark <- tempBetaMark + deltaDrop2Sum[i]*(assocDrop2[1]*xTimeDrop2[i,]+assocDrop2[2]*xTimeDerDrop2[i,])
      
      xTimeVisForInti <- xTimeVisForInt[jjVisFI,]
      xTimeDerVisForInti <- xTimeDerVisForInt[jjVisFI,]
      xTimeDrop1ForInti <- xTimeDrop1ForInt[jjDrop1,]
      xTimeDerDrop1ForInti <- xTimeDerDrop1ForInt[jjDrop1,]
      xTimeDrop2ForInti <- xTimeDrop2ForInt[jjDrop2,]
      xTimeDerDrop2ForInti <- xTimeDerDrop2ForInt[jjDrop2,]
      
      for (j in 1:pLong)
      {
        cSCbetaMark <- tempBetaMark[j,] - colSums((assocVis[1]*xTimeVisForInti[,j] + 
                                                   assocVis[2]*xTimeDerVisForInti[,j])*xxVis)
        cSCbetaMark <- cSCbetaMark - colSums((assocDrop1[1]*xTimeDrop1ForInti[,j] + 
                                              assocDrop1[2]*xTimeDerDrop1ForInti[,j])*xxDrop1)
        cSCbetaMark <- cSCbetaMark - colSums((assocDrop2[1]*xTimeDrop2ForInti[,j] + 
                                              assocDrop2[2]*xTimeDerDrop2ForInti[,j])*xxDrop2)
        SCbetaMark[i,j] <- pi^(-qLong/2)*sum(weights*(cSCbetaMark*RCondMarg))
      }
      
      # Score vector for vecD
      # matrix of cross-products for each quadrature point
      iter <- 0
      for (j in 1:qLong)
      {
        for (k in 1:qLong)
        {
          iter <- iter + 1
          crossRgrid[,iter] <- rgrid[,j]*rgrid[,k]
        }
      }
      cSCvecD <- 0.5*(crossRgrid %*% DreInvKron)
      SCvecD[i,] <- pi^(-qLong/2)*colSums(weights*(cSCvecD*RCondMarg)) 
      SCvecD[i,] <- SCvecD[i,] - 0.5*pi^(-qLong/2)*sum(weights*(RCondMarg))*c(DreInv)
    }
    matSCvecD <- matrix(colSums(SCvecD),nr = qLong,nc = qLong,byrow = F)
    SCvechLogD <- crossprod(vecEYvechY(vechLogDre),c(matSCvecD))
    
    # Score vectors of all model parameters
    score <- -c(colSums(SCbetaMark),sum(SClogSigma),SCvechLogD,
                colSums(cbind(SCbetasVis,SCassocVis,SCpsiVis,
                              SCbetasDrop1,SCassocDrop1,SCpsiDrop1,
                              SCbetasDrop2,SCassocDrop2,SCpsiDrop2)))
    score
    
  }
  
  # Gradient/Hessian of log-likelihood
  fdGrad <- function(x,n=length(x)){fd(x,logMarg,n)}
  fdHess <- function(x){fd.vec(x,ScoreVec)}
  
  cdGrad <- function(x){cd(x,logMarg)}
  cdHess <- function(x){cd.vec(x,cdGrad)}
  
  # Difference between successive parameter values and iteration index
  diffCoef = Inf
  it = 0
  
  # Info for each iteration of the BFGS algorithm
  BFGShist = matrix(nrow = maxit + 1,ncol = 3)
  colnames(BFGShist) = c("Iteration","Diff_coef","logl")
  BFGShist[1,] = c(it,diffCoef,-logMarg(paramStart))
  
  # Parameter values for the longitudinal model during the BFGS
  coefHistLong = matrix(nrow = maxit + 1,ncol = 1 + nParLong )
  namesCoefLong = c(names(fixed.effects(fitlme)),"log(sigma^2)",
                    paste0("Vech(logD)",1:(qLong*(qLong + 1)/2))) 
  colnames(coefHistLong) = c("Iteration",namesCoefLong)
  coefHistLong[1,] = c(it,paramStart[1:nParLong])
  
  # History of coefficients for the dropout model
  coefHistVis = matrix(nrow = maxit + 1,ncol = 1 + lVis)
  namesCoefVis = c( c(names(coef(fitCoxVis)),paste0("AssocVis",1:lassocVis),
                      paste0("PsiVis",1:lpsiVis) )) 
  colnames(coefHistVis) = c("Iteration",namesCoefVis)
  coefHistVis[1,] = c(it,paramVis)
  
  # History of coefficients for the dropout model (cause 1)
  coefHistDrop1 = matrix(nrow = maxit + 1,ncol = 1 + lDrop1)
  namesCoefDrop1 = c( c(names(coef(fitCoxDrop1)),paste0("AssocDrop1",1:lassocDrop1),
                       paste0("PsiDrop1",1:lpsiDrop1) )) 
  colnames(coefHistDrop1) = c("Iteration",namesCoefDrop1)
  coefHistDrop1[1,] = c(it,paramDrop1)
  
  # History of coefficients for the dropout model (cause 2)
  coefHistDrop2 = matrix(nrow = maxit + 1,ncol = 1 + lDrop2)
  namesCoefDrop2 = c( c(names(coef(fitCoxDrop2)),paste0("AssocDrop2",1:lassocDrop2),
                        paste0("PsiDrop2",1:lpsiDrop2) )) 
  colnames(coefHistDrop2) = c("Iteration",namesCoefDrop2)
  coefHistDrop2[1,] = c(it,paramDrop2)
  
  # Initial approximation for the Hessian matrix
  approxhess <- matrix(0,length(paramStart),length(paramStart))
  approxhess[c(1:nParLong),c(1:nParLong)] <- fd.vec(paramStart[1:nParLong],function(x) fd(x,logMarker) )
  
  # Hessian for the visiting process
  hessVis <- fd.vec(paramVis,fdGradVis)
  aa <- try(chol(hessVis),silent = T) # Make it diagonal if not positive-definite
  if ("try-error" %in% class(aa)){hessVis <- diag(abs(diag(hessVis)))}
  
  # Hessian for the dropout process
  hessDrop1 <- fd.vec(paramDrop1,fdGradDrop1)
  aa <- try(chol(hessDrop1),silent = T) # Make it diagonal if not positive-definite
  if ("try-error" %in% class(aa)){hessDrop1 <- diag(abs(diag(hessDrop1)))}
  
  # Hessian for the dropout process (cause 2)
  hessDrop2 <- fd.vec(paramDrop2,fdGradDrop2)
  aa <- try(chol(hessDrop2),silent = T) # Make it diagonal if not positive-definite
  if ("try-error" %in% class(aa)){hessDrop2 <- diag(abs(diag(hessDrop2)))}
  
  # Marker block unchanged, set the off-diagonal block equal to zero
  approxhess[-c(1:nParLong),-c(1:nParLong)] <- 0
  approxhess[-c(1:nParLong),] <- 0
  approxhess[,-c(1:nParLong)] <- 0
  approxhess[(nParLong + 1):(nParLong + lVis),(nParLong + 1):(nParLong + lVis)] <- hessVis
  approxhess[(nParLong + lVis+1):(nParLong + lVis + lDrop1),
             (nParLong + lVis+1):(nParLong + lVis + lDrop1)] <- hessDrop1
  approxhess[-(1:(nParLong + lVis + lDrop1)),
             -(1:(nParLong + lVis + lDrop1))] <- hessDrop2
  
  paramNext <- paramStart
  gradientNext = ScoreVec(paramNext)
  
  # Begin the NR
  while (diffCoef > epsilon & it < maxit)
  {
    # Update iteration index
    it = it + 1
    param = paramNext
    print(param)
    
    # Update the parameters
    paramNext = try(param - solve(approxhess,gradientNext),silent = T)
    
    if ("try-error" %in% class(paramNext))
    {
      break;
    }
    
    # Evaluate the gradient
    gradient = gradientNext
    gradientNext = ScoreVec(paramNext)
    
    if ( any(!is.finite(gradientNext)) )
    {
      break;
    }
    
    s = (paramNext - param)
    psi = (gradientNext - gradient)
    
    # Update the Hessian matrix
    approxhess = approxhess + tcrossprod(psi)/sum(psi*s) - (approxhess%*%tcrossprod(s)%*%approxhess)/c(crossprod(s,approxhess) %*% s)
    
    # Difference in coefficients
    diffCoef = sqrt( sum( (param-paramNext)^2 ) )
    
    BFGShist[it+1,] = c(it,diffCoef,logMarg(paramNext))
    coefHistLong[it+1,] = c(it,paramNext[1:nParLong])
    coefHistVis[it+1,] <- c(it,paramNext[(nParLong+1):(nParLong+lVis)])
    coefHistDrop1[it+1,] <- c(it,paramNext[(nParLong + lVis+1):(nParLong + lVis + lDrop1)])
    coefHistDrop2[it+1,] <- c(it,paramNext[-(1:(nParLong + lVis + lDrop1))])
  }
  
  if (any(!is.finite(gradientNext)) | diffCoef > epsilon | "try-error" %in% class(paramNext))
  {
    fitConstr <- optim(paramStart,fn = logMarg,gr = ScoreVec,control = list(trace = 100),method = "BFGS")
    paramNext <- fitConstr$par 
    
    param <- paramNext
    # Hessian matrix at final values
    HessianFinal <- fdHess(param)
    VarCovFinal <- solve(HessianFinal)
  } else{
    param <- paramNext
    # Hessian matrix at final values
    HessianFinal <- fdHess(param)
    VarCovFinal <- solve(HessianFinal)
  }
  
  # History of iterations
  BFGShist = BFGShist[1:(it+1),]
  coefHistVis = coefHistVis[1:(it+1),]
  coefHistDrop1 = coefHistDrop1[1:(it+1),]
  coefHistDrop2 = coefHistDrop2[1:(it+1),]
  coefHistLong = coefHistLong[1:(it+1),]
  
  # Matrix logarithm of Dre
  vechLogDre <- param[(pLong + 2):nParLong] # vech(log(Dre)) 
  logDre <- matrix(Gqq %*% vechLogDre,nrow = qLong,ncol = qLong,byrow = F)
  decom <- eigen(logDre)
  
  # Variance-covariance of the random effects
  Dre <- decom$vectors %*% diag(exp(decom$values)) %*% t(decom$vectors)
  
  # Results for the longitudinal model
  sumLong = matrix(nr = nParLong,nc = 4)
  rownames(sumLong) = namesCoefLong
  colnames(sumLong) = c("Estimate","SE","LB","UB")
  
  sumLong[,1] = param[c(1:nParLong)]
  sumLong[,2] = sqrt(diag(VarCovFinal[1:nParLong,1:nParLong]))
  sumLong[,3] = sumLong[,1] - qnorm(1-alpha/2)*sumLong[,2]
  sumLong[,4] = sumLong[,1] + qnorm(1-alpha/2)*sumLong[,2]
  
  # Results for the visiting model
  sumVis = matrix(nr = lVis,nc = 4)
  rownames(sumVis) = namesCoefVis
  colnames(sumVis) = c("Estimate","SE","LB","UB")
  
  sumVis[,1] = param[c(nParLong+1):(nParLong+lVis)]
  sumVis[,2] = sqrt(diag(VarCovFinal[(nParLong+1):(nParLong+lVis),(nParLong+1):(nParLong+lVis)]))
  sumVis[,3] = sumVis[,1] - qnorm(1-alpha/2)*sumVis[,2]
  sumVis[,4] = sumVis[,1] + qnorm(1-alpha/2)*sumVis[,2]
  
  # Results for the dropout model (cause 1)
  sumDrop1 = matrix(nr = lDrop1,nc = 4)
  rownames(sumDrop1) = namesCoefDrop1
  colnames(sumDrop1) = c("Estimate","SE","LB","UB")
  
  sumDrop1[,1] = param[(nParLong + lVis+1):(nParLong + lVis + lDrop1)]
  sumDrop1[,2] = sqrt(diag(VarCovFinal[(nParLong + lVis+1):(nParLong + lVis + lDrop1),
                                       (nParLong + lVis+1):(nParLong + lVis + lDrop1)]))
  sumDrop1[,3] = sumDrop1[,1] - qnorm(1-alpha/2)*sumDrop1[,2]
  sumDrop1[,4] = sumDrop1[,1] + qnorm(1-alpha/2)*sumDrop1[,2]
  
  # Results for the dropout model (cause 2)
  sumDrop2 = matrix(nr = lDrop2,nc = 4)
  rownames(sumDrop2) = namesCoefDrop2
  colnames(sumDrop2) = c("Estimate","SE","LB","UB")
  
  sumDrop2[,1] = param[-(1:(nParLong + lVis + lDrop1))]
  sumDrop2[,2] = sqrt(diag(VarCovFinal[-(1:(nParLong + lVis + lDrop1)),-(1:(nParLong + lVis + lDrop1))]))
  sumDrop2[,3] = sumDrop2[,1] - qnorm(1-alpha/2)*sumDrop2[,2]
  sumDrop2[,4] = sumDrop2[,1] + qnorm(1-alpha/2)*sumDrop2[,2]
  
  # Prepare output
  out = list()
  
  if (diffCoef > epsilon)
  {
    out$note = paste("Convergence failed at",epsilon,
                     "Euclidean distance for coefficients",
                     ", may need more iterations than",maxit)
  } else {out$note = "Successful convergence"}
  
  out$it = it
  out$BFGShist = BFGShist
  out$coefHistVis = coefHistVis
  out$coefHistDrop1 = coefHistDrop1
  out$coefHistDrop2 = coefHistDrop2
  out$coefHistLong = coefHistLong
  out$sumVis = sumVis
  out$sumDrop1 = sumDrop1
  out$sumDrop2 = sumDrop2
  out$sumLong = sumLong
  out$var = VarCovFinal
  out$hessian = HessianFinal
  out$approxhess = approxhess
  
  out$IntKnotsVis = IntKnotsVis
  out$IntKnotsDrop1 = IntKnotsDrop1
  out$IntKnotsDrop2 = IntKnotsDrop2
  out$bKnotsVis = bKnotsVis
  out$bKnotsDrop1 = bKnotsDrop1
  out$bKnotsDrop2 = bKnotsDrop2
  out$param = param
  out$gradient = gradientNext
  out$Dre <- Dre
  
  return(out)
  
}