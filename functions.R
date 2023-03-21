## Functions for replicating smart-meter analysis
# Alex Gibberd

evalForecast <- function(pred,true,type="RMSE"){
  p = dim(pred)[2]
  n = dim(pred)[1]
  d = true-pred
  res = vector()
  for (j in 1:p){
    notna = !is.na(d[,j])
    if(type=="RMSE"){
      res[j] = sqrt((1/n)*sum(d[notna,j]^2))
    }else if(type=="MAE"){
      res[j] = (1/n)*sum(abs(d[notna,j]))
    }else if(type=="MAX"){
      res[j] = max(abs(d[notna,j]))
    }else if(type=="MAPE"){
      res[j] = 0
      for (i in 1:n){
        if( !is.na(d[i,j]) & (true[i,j]!=0)){
          res[j] = res[j] + abs(100*d[i,j]/true[i,j])
        }
      }
      res[j] = res[j]/n
    }
  }
  return(res)
}


fillNA2 <-function(X){
  # Taken from sparseDFM
  n <- dim(X)[1]
  p <- dim(X)[2]
  k <- 3
  idx.na <- is.na(X)
  
  for (i in 1:p){  
    x = X[,i]
    na_x = is.na(x)
    t1 = min(which(!na_x))
    t2 = max(which(!na_x))
    
    x1 <- stats::spline(x[t1:t2],xout = 1:(t2-t1+1))
    xx <- x1$y
    x[t1:t2] <- x1$y
    na_x <- is.na(x)
    x[na_x] <- median(x,na.rm = T)
    
    x_MA3 <- stats::filter(x = c(rep(x[1],k),x,rep(x[length(x)],k)),filter = rep(1,2*k+1)/(2*k+1),sides = 1)
    
    x_MA3 = x_MA3[(2*k+1):length(x_MA3)]
    x[idx.na[,i]] = x_MA3[idx.na[,i]]
    X[,i] = x
  }
  
  return(X)
}

arimaForecast <- function(data_train,data_test,par,qar,h){
  
  ntest = dim(data_test)[1]
  p = dim(data_test)[2]
  n = dim(data_train)[1]
  xhat = matrix(NA,nrow=ntest,ncol=p)
  XComb = rbind(data_train,data_test)
  XComb = fillNA2(XComb)
  for (i in 1:(ntest-h)){
    for (j in 1:p){
      mod = arima(XComb[1:(n+i),j],order=c(par,0,qar))
      xhat[i,j] = as.numeric(predict(mod,n.ahead=h)$pred[h])
    }
  }
  return(xhat)
}

stepForecast <- function(mod,data_train,data_test,h){
  
  # Takes in DFM model parameters and produces h step ahead forecasts
  # mod - sparseDFM object
  # data_train - nxp training dataset
  # data_test - nxp test dataset
  # h - step ahead
  
  n = dim(data_train)[1]
  p = dim(data_train)[2]
  
  a0 = mod$params$a0_0
  P0_0 = mod$params$P0_0
  A = mod$params$A
  Lambda = mod$params$Lambda
  Sigma_u = mod$params$Sigma_u
  Sigma_eps = diag(mod$params$Sigma_epsilon)
  
  mu = mod$data$X.mean
  std = mod$data$X.sd
  # Scale to standard of training data
  for (i in 1:p){
    data_train[,i] = (data_train[,i]- mu[i])/std[i]
    data_test[,i] = (data_test[,i]- mu[i])/std[i]
  }
  ntest = dim(data_test)[1]
  Xhat = array(NA, dim=c(ntest,p))
  XCov = array(NA, dim=c(ntest,p,p))
  Xcomb = rbind(data_train,data_test)
  for (i in 1:(ntest-h)){
    q = n+i+h
    Xtest = Xcomb[1:(n+i),]
    for (j in 1:h){
      Xtest = rbind(Xtest,rep(NA,p))
    }
    Xtest = as.matrix(Xtest)
    res = kalmanUnivariate(Xtest,a0,P0_0,A,Lambda,Sigma_eps,Sigma_u)
    Ft = res$at_n[,q]
    FtCov = res$Pt_n[,,q]
    Xhat[i,] = Lambda %*% Ft 
    XCov[i,,] = Lambda %*% FtCov %*% t(Lambda) + Sigma_eps
  }
  # Scale to original scale
  for (i in 1:p){
    Xhat[,i] = (Xhat[,i] * std[i]) + mu[i]
    for (j in 1:p){
      XCov[,i,j] = XCov[,i,j] * std[i] * std[j]
    }
  }
  return(list("Xhat"=Xhat,"XCov"=XCov,"test"=data_test,"train"=data_train))
}

plotDayFactor <- function(data,period,newlength=24,avg=FALSE,col="black",lty=1,ylim=NULL){
  xtick = newlength*(1:period)/period
  n = dim(data)[1]
  nperiod = (n/period)
  p = dim(data)[2]
  tmp = array(data=NA,dim=c(p,nperiod,period))
  for (j in 1:p){
    for(i in 1:nperiod){
      s = (i-1)*(period)+1
      e = i*period
      tmp[j,i,]=data[s:e,j]
    }
  }
  favg = apply(tmp,c(1,3),mean)
  fstd = apply(tmp,c(1,3),sd)
  favg.ts = ts(t(favg),deltat =24/period)
  fcu = ts(t(favg)+1.96*t(fstd),deltat=24/period)
  fcl = ts(t(favg)-1.96*t(fstd),deltat=24/period)
  
  oldpar <- par(mar = c(0, 5, 0, 2), oma = c(6, 0, 5, 0), mfrow = c(p, 1L))
  on.exit(par(oldpar))
  if(is.null(ylim)){
    ylim = c(min(fcl),max(fcu))
  }
  for(i in 1:p) {
    
    plot(xtick,favg.ts[,i], col = col, type="l", lty=lty, lwd=2,axes = FALSE, xlab = "", ylab = "",ylim=ylim)
    #lines(fcu[,i], type = 'l', col = col, lwd = 1,lty=3)
    #lines(fcl[,i], type = 'l', col = col, lwd = 1, lty=3)
    # Polygon
    x <- c(xtick, rev(xtick))
    y <- c(fcu[,i], rev(fcl[,i]))
    polygon(x, y, col =  adjustcolor(col, alpha.f = 0.10), border = NA)
    
    box()
    axis(2, cex.axis = 1.2)
    if(i == p) axis(1, cex.axis = 1.2)
    mtext(paste("Factor", i), 2, line = 3, cex = 1.2)
    if(i == p) mtext('Time' , side = 1, line = 3, cex = 1.2)
    
    
  }
  
  return(list("favg"=favg,"fstd"=fstd))
  
}

downSample <- function(data, period=6){
  # Averages across intervals of period
  # Assumes first column is timestamp
  p = dim(data)[2]-1
  n = dim(data)[1]
  n2 = floor(n/period)
  data2 = array(NA,dim=c(n2,p+1))
  data2 = as.data.frame(data2)
  colnames(data2) = colnames(data)
  for (t in 1:n2){
    s = (t-1)*(period)+1
    e = (t*period)
    data2[t,1] = data[s,1]
    for (i in 1:p){
      # Multiply by period to get kw/period
      data2[t,i+1] = period*mean(as.numeric(data[s:e,i+1]))
    }
  }
  return(data2)
}

weekendIdx <- function(n,dayobs){
  wknd = 1:n
  for (i in 1:n){
    day = floor(i/dayobs)+1
    if((day %% 7) %in% c(6,0) ){
      wknd[i] = TRUE
    }else{
      wknd[i] = FALSE
    }
  }
  return(wknd)
}

trimNA <- function(X,thresh,exclude=NULL){
  # Trims some columns based on missingness and non-unique entries
  p = dim(X)[2]
  missing = rep(FALSE,p)
  R = c()
  for (i in 1:p){
    uvals = X[which(!is.na(X[,i])),i]
    uvals = sum(unique(uvals))
    exc.cond = (exclude %in% colnames(data)[i])
    if((sum(is.na(X[,i])) <= thresh) & (uvals > 1) & !exc.cond){
      R = cbind(R,X[,i])
    }else{
      missing[i] = TRUE
    }
  }
  return(list("X"=R,"missing"=missing))
}

svarForecast <- function(lambda,data_train,data_test,h){
  ntest = dim(data_test)[1]
  p = dim(data_test)[2]
  n = dim(data_train)[1]
  Xhat = matrix(NA,nrow=ntest,ncol=p)
  XComb = rbind(data_train,data_test)
  XComb = as.matrix(fillNA2(XComb))
  A <- BigVAR.fit(XComb[1:n,], p=1, "Basic", lambda=lambda, intercept=TRUE)
  
  for (i in 1:(ntest-h)){
    lastX = XComb[(n+i-1),]
    Ahat <- matrix(A[,2:(p+1),1],nrow=p,ncol=p)
    # Project forwards
    pred <- Ahat %*% lastX
    if(h>1){
      for (k in 1:(h-1)){
        Ahat %*% pred
      }
    }
    # Add mean
    pred = pred + A[,1,1]
    Xhat[i,] = pred
  }
  return(Xhat)
  
}