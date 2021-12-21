\name{HazardDiff}
\alias{HazardDiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Conditional Treatment Effect for Competing Risks}
\description{
This function estimates the conditional treatment effect for competing risks data in observational studies. The effect is described as a constant difference between the hazard functions given the covariates. The specific functional form for the covariates does not need to be assumed.}
\usage{
HazardDiff(data,type,ps,CL=.95,G="NonSpecified")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{DataFrame of survival data. First column needs to be the censored failure time and needs to be named 'time', second column needs to be the type of failure/event indicator (1=observed failure from primary risk, 2=observed failure from competing risk, 0=censored) and needs to be named 'status', third column needs to be the exposure (1 or 0) and the other columns are covariates.}
  \item{type}{If type=1, score 1 is used, if type=2, score 2 is used.}
  \item{CL}{Confidence level of the confidence interval.}
  \item{ps}{A vector of propensity scores, one for each observation.}
  \item{G}{a nxn dataframe with the censoring survival function. n is the number of observations in data. Every row of G corresponds to a different observation. Every column of G corresponds to a different censored failure time of the column 'time' of data. If not specified, the censoring survival function is not considered. If specified, type needs to be equal to 1.}
}
\details{
The conditional treatment effect on the primary risk
\deqn{\beta_1=\lambda_1(t | A=1, Z)-\lambda_1(t | A=0, Z)}{\beta~1=\lambda_1(t | A=1, Z)-\lambda_1(t | A=0, Z)} and the conditional treatment effect on the competing risk \deqn{\beta_2=\lambda_2(t | A=1, Z)-\lambda_2(t | A=0, Z)}{\beta_2=\lambda_2(t | A=1, Z)-\lambda_2(t | A=0, Z)} are estimated.
}
\value{
\item{beta}{A vector with two components. First component is the estimate of the conditional treatment effect on the primary risk. Second component is the estimate of the conditional treatment effect on the competing risk.}
\item{se}{A vector with two components. First component is the standard error of the estimator of the conditional treatment effect on the primary risk. Second component is the standard error of the estimator of the conditional treatment effect on the competing risk.}
\item{CI}{A matrix with two rows. First row contains the bounds of the confidence interval for the conditional treatment effect on the risk. Second row contains the bounds of the confidence interval for the conditional treatment effect on the competing risk.}
}
\references{
Rava, D., Xu, R. "Doubly Robust Estimation of the Hazard Difference for Competing Risks Data", December 2021, arXiv:2112.09535}

\author{Denise Rava and Ronghui Xu}

\examples{
#simulating data
n=1000
Z1=1:n
Z2=Z1
A=Z1
time=Z1
status=Z1
x=status
c=status
beta1=0.1
beta2=0.1
for (i in 1:n)
{Z1[i]=runif(1,0,0.5)
Z2[i]=runif(1,0,0.5)
predt<-1/(1+exp(-Z1[i]-Z2[i]))
A[i]=rbinom(1,1,predt)
u=runif(1,0,1)
time[i]=-log(u)/(1.5+1.5*((Z1[i]))+1.5*((Z2[i]))+(beta1+beta2)*A[i])
status[i]=rbinom(1,1,1-((.5+.5*((Z1[i]))+.5*(Z2[i]))+(beta1)*A[i])/(1.5
+1.5*((Z1[i]))+1.5*((Z2[i]))+(beta1+beta2)*A[i]))+1
u=runif(1,0,1)
c[i]=runif(1,0,3)
x[i]=apply(cbind(time[i],c[i]),1,min)
status[i]<-ifelse(time[i]<c[i],status[i],0)}

data<-data.frame(time=x,status=status,A=A,Z1=Z1,Z2=Z2)

#estimating propensity score
mylogit <- glm(A ~ 1+Z1+Z2,data,family = binomial(link="logit"))
ps=predict(mylogit,newdata = data,type='response')

#score 1
HazardDiff(data=data, type=1, ps=ps)
#score 2
\donttest{HazardDiff(data=data, type=2, ps=ps)}


}


\keyword{ ~HazardDifference }
\keyword{ ~TreatmentEffect }
