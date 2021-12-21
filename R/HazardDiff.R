#' @title HazardDiff
#'
#' @description This function estimates the conditional treatment effect for competing risks data in observational studies. While it is described as a constant difference between the hazard functions given the covariates, we do not assume specific functional forms for the covariates.
#'
#' @param data=data frame of survival data. First column needs to be the censored failure time and needs to be named 'time', second column needs to be the failure indicator (1=observed failure, 2=observed competing risk, 0=censored) and needs to be named 'status', third column needs to be the exposure (1 or 0) and the other columns are covariates,
#' @param type=type if type is 1, score 1 is used, if type is 2, score 2 is used.
#' @param CL=CL confidence level of the confidence interval
#' @param ps=ps a vector of propensity score for each observations
#' @param G=G a nxn data.frame with the censoring survival function. n is the number of observations in data. Every row of the data.frame corresponds to a different observation. every column corresponds to a different time. Times need to be the same censored failure time of data. If not specified, the censoring survival function is not considered. If specified, type needs to be equal to 1.
#'
#' @return beta=estimated treatment effect, se=standard error, CI=confidence interval
#'
#' @examples
#'
#' @export HazardDiff

HazardDiff<-function(data, type, ps, CL=.95, G='NonSpecified')
{#sorting
  train<-data[order(data$time),]

  #estimating gamma
  classical1<-ahaz::ahaz(survival::Surv(train$time,(train$status==1)),train[,-c(1:2)])
  gamma1hat<-stats::coefficients(classical1)[-1]
  classical2<-ahaz::ahaz(survival::Surv(train$time,(train$status==2)),train[,-c(1:2)])
  gamma2hat<-stats::coefficients(classical2)[-1]

  #if score 1 is chosen
  if (!is.numeric(G))
  {n=nrow(train)
  G=matrix(rep(1,n^2),n,n)}

  if (type==1)
  {w=train[,3]*(1-ps)
  beta<-Score1(train,w,ps,gamma1hat,gamma2hat,G)
  beta1<-beta[1]
  beta2<-beta[2]
  sd1<-beta[3]
  sd2<-beta[4]


  c1r<-beta1-stats::qnorm((1-CL)/2)*sd1
  c1l<-beta1+stats::qnorm((1-CL)/2)*sd1
  c2r<-beta2-stats::qnorm((1-CL)/2)*sd2
  c2l<-beta2+stats::qnorm((1-CL)/2)*sd2

  }

  #if score 2 is chosen
  if (type==2)
  {fun<-function(x){Score2(x,train,ps,gamma1hat,gamma2hat)}
  beta<-rootSolve::multiroot(fun, c(0,0))
  beta1<-beta$root[1]
  beta2<-beta$root[2]

  sd<-Score2sd(train,c(beta1,beta2),ps)
  sd1<-sd[1]
  sd2<-sd[2]

  c1r<-beta1-stats::qnorm((1-CL)/2)*sd1
  c1l<-beta1+stats::qnorm((1-CL)/2)*sd1
  c2r<-beta2-stats::qnorm((1-CL)/2)*sd2
  c2l<-beta2+stats::qnorm((1-CL)/2)*sd2
  }

  CI=data.frame(lower=c(c1l,c2l),upper=c(c1r,c2r))
  row.names(CI)=c('beta1','beta2')

  return(list(beta=c(beta1,beta2),se=c(sd1,sd2),CI=CI))
}

#score 1
Score1<-function(train,w,pred,gamma1,gamma2,sc)
{n=nrow(train)
Z=data.frame(train[,-c(1:3)])
Zbar=matrix(n*ncol(Z),n,ncol(Z))
for (i in 1:ncol(Z))
{Zbar[,i]=rev(cumsum(rev(w*Z[,i])))/rev(cumsum(rev(w)))}
Zbar[is.nan(Zbar)]<-0
denNbar=1/rev(cumsum(rev(w)))
denNbar[(denNbar==Inf)]<-0
scbeta1<-1:nrow(train)
scbeta2<-1:nrow(train)
den=1:nrow(train)
difft<-diff(c(0,train$time))
for (i in 1:nrow(train))
{scbeta1[i]<--((1-train$A[i])*pred[i]*((train$status[i]==1)/sc[i,i]-cumsum(1/sc[i,]*difft)[i]*as.numeric(gamma1%*%t(Z[i,]))+cumsum(as.numeric(gamma1%*%t(Zbar))/sc[i,]*difft)[i]-cumsum((train$status==1)*w*denNbar/sc[i,])[i]))
scbeta2[i]<--((1-train$A[i])*pred[i]*((train$status[i]==2)/sc[i,i]-cumsum(difft/sc[i,])[i]*as.numeric(gamma2%*%t(Z[i,]))+cumsum(as.numeric(gamma2%*%t(Zbar))*difft/sc[i,])[i]-cumsum((train$status==2)*w*denNbar/sc[i,])[i]))
den[i]<-cumsum((1-train$A[i])*pred[i]/sc[i,]*difft)[i]}

beta1<-sum(scbeta1)/sum(den)
beta2<-sum(scbeta2)/sum(den)
beta<-c(beta1,beta2)
sigmat<-1:nrow(train)
timetemp=c(0,train$time[-length(train$time)])
for (i in 1:nrow(train))
{sigmat[i]=(train$A[i]*(train$A[i]-pred[i])*cumsum(1/(beta1+beta2)*(exp((beta1+beta2)*train$time)-(exp((beta1+beta2)*timetemp)))/sc[i,])[i])}
sigma=sum(sigmat)
v1=sigma^{-1}*sqrt(sum((train$status==1)*exp(2*(beta1+beta2)*train$time*train$A)*(train$A-pred)^2*(diag(sc))^(-2)))
v2=sigma^{-1}*sqrt(sum((train$status==2)*exp(2*(beta1+beta2)*train$time*train$A)*(train$A-pred)^2*(diag(sc))^(-2)))
return(c(beta,v1,v2))}

#score 2
Score2<-function(beta,train,pred,gamma1,gamma2)
{cov<-train[,-c(1:3)]
n=nrow(train)
beta1<-beta[1]
beta2<-beta[2]
Abar<-rev(cumsum(rev(train$A)))/n:1
Eps<-(exp(-(beta1+beta2)*train$time)*pred)/(exp(-(beta1+beta2)*train$time)*pred+(1-pred))

sc1t<-1:nrow(train)
sc2t<-1:nrow(train)
difft<-diff(c(0,train$time))

temptime<-c(0,train$time)
temptemp<-matrix(1:nrow(train)^2,nrow(train),nrow(train))
for (l in 1:nrow(train))
{temptemp[l,]<-1/(n-l+1)*(log(abs(exp(-(beta1+beta2)*train$time[l])*pred+1-pred))-log(abs(exp(-(beta1+beta2)*temptime[l])*pred+1-pred)))}

if (beta1+beta2==0)
{temptime<-c(0,train$time)
temptemp<-matrix(1:nrow(train)^2,nrow(train),nrow(train))
for (l in 1:nrow(train))
{temptemp[l,]<-pred/(n-l+1)*(train$time[l]-temptime[l])}

for (i in 1:nrow(train))
{Epstemp<-(exp(-(beta1+beta2)*train$time[i])*pred)/(exp(-(beta1+beta2)*train$time[i])*pred+(1-pred))
Epsbar<-rev(cumsum(rev(Epstemp)))/n:1

funtemp<-function(j) {return(sum(temptemp[1:min(i,j),j]))}
temp<-c(unlist(lapply(1:nrow(train),funtemp)))

sc1t[i]<-(train$status[i]==1)*train$A[i]-train$A[i]*train$time[i]*(beta1*train$A[i]+c(gamma1%*%t(cov[i,])))-(train$status[i]==1)*Abar[i]+(beta1*train$A[i]+c(gamma1%*%t(cov[i,])))*cumsum(Abar*difft)[i]-(train$status[i]==1)*Eps[i]+(beta1*train$A[i]+c(gamma1%*%t(cov[i,])))*(pred[i]*train$time[i])+(train$status[i]==1)*Epsbar[i]-(beta1*train$A[i]+c(gamma1%*%t(cov[i,])))*sum(temp)
sc2t[i]<-(train$status[i]==2)*train$A[i]-train$A[i]*train$time[i]*(beta2*train$A[i]+c(gamma2%*%t(cov[i,])))-(train$status[i]==2)*Abar[i]+(beta2*train$A[i]+c(gamma2%*%t(cov[i,])))*cumsum(Abar*difft)[i]-(train$status[i]==2)*Eps[i]+(beta2*train$A[i]+c(gamma2%*%t(cov[i,])))*(pred[i]*train$time[i])+(train$status[i]==2)*Epsbar[i]-(beta2*train$A[i]+c(gamma2%*%t(cov[i,])))*sum(temp)}}


if (beta1+beta2!=0)
{for (i in 1:nrow(train))
{Epstemp<-(exp(-(beta1+beta2)*train$time[i])*pred)/(exp(-(beta1+beta2)*train$time[i])*pred+(1-pred))
Epsbar<-rev(cumsum(rev(Epstemp)))/n:1

funtemp<-function(j) {return(sum(temptemp[1:min(i,j),j]))}
temp<-c(unlist(lapply(1:nrow(train),funtemp)))

sc1t[i]<-(train$status[i]==1)*train$A[i]-train$A[i]*train$time[i]*(beta1*train$A[i]+c(gamma1%*%t(cov[i,])))-(train$status[i]==1)*Abar[i]+(beta1*train$A[i]+c(gamma1%*%t(cov[i,])))*cumsum(Abar*difft)[i]-(train$status[i]==1)*Eps[i]-(beta1*train$A[i]+c(gamma1%*%t(cov[i,])))*(1/(beta1+beta2)*log(abs(exp(-(beta1+beta2)*train$time[i])*pred[i]+1-pred[i])))+(train$status[i]==1)*Epsbar[i]+(beta1*train$A[i]+c(gamma1%*%t(cov[i,])))*sum(temp)/(beta1+beta2)
sc2t[i]<-(train$status[i]==2)*train$A[i]-train$A[i]*train$time[i]*(beta2*train$A[i]+c(gamma2%*%t(cov[i,])))-(train$status[i]==2)*Abar[i]+(beta2*train$A[i]+c(gamma2%*%t(cov[i,])))*cumsum(Abar*difft)[i]-(train$status[i]==2)*Eps[i]-(beta2*train$A[i]+c(gamma2%*%t(cov[i,])))*(1/(beta1+beta2)*log(abs(exp(-(beta1+beta2)*train$time[i])*pred[i]+1-pred[i])))+(train$status[i]==2)*Epsbar[i]+(beta2*train$A[i]+c(gamma2%*%t(cov[i,])))*sum(temp)/(beta1+beta2)}}




sc1<-mean(sc1t)
sc2<-mean(sc2t)
score<-c(sc1,sc2)
return(score)}


#SE for score 2
Score2sd<-function(train,beta,pred)
{n=nrow(train)
beta1<-beta[1]
beta2<-beta[2]
Abar<-rev(cumsum(rev(train$A)))/n:1
Eps<-(exp(-(beta1+beta2)*train$time)*pred)/(exp(-(beta1+beta2)*train$time)*pred+(1-pred))

v1=sum((train$status==1)*(train$A-Eps)^2)
v2=sum((train$status==2)*(train$A-Eps)^2)
sigma=sum(train$A*(train$A*train$time+1/(beta1+beta2)*log(abs(exp(-(beta1+beta2)*train$time)*pred+1-pred))))

sd1=sigma^{-1}*sqrt(v1)
sd2=sigma^{-1}*sqrt(v2)
return(c(sd1,sd2))
}








