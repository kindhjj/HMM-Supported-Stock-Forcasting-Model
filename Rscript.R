#---------------R code for stock price analysis and prediction -- written by He Junjie-----------
library(readxl)
library(tseries)
library(forecast)
library(rugarch)
library(depmixS4)
library(quantmod)
library(rmgarch)
library(ramify)
library(MTS)
library(vars)

#load data
data=read_excel("case_study_data.xlsx",sheet=1)
attach(data)

#analyze data
plot(close,type='l')
diff_close=diff(close)
log_close=diff(log(close))
plot(log_close,type='l')
adf.test(log_close)
acf(log_close)
Box.test(log_close,type = 'Ljung')
adf.test(log_close^2)
Box.test(log_close^2,type = 'Ljung')
summary(log_close)
summary(log_close^2)
plot(density(log_close),main="Daily Return")
xplot=seq(min(log_close),max(log_close),length=length(log_close))
yplot=dnorm(xplot,mean=mean(log_close),sd=sd(log_close))
lines(xplot,yplot,col="blue")
qqnorm(log_close,datax=T)
qqline(log_close,datax=T)
frac_change=(close-open)/open
frac_high=(high-open)/open
frac_low=(open-low)/open
frac_vec=as.matrix(cbind(frac_change,frac_high,frac_low))
ccm(frac_vec) #cross-correlation matrix

length=length(log_close)
test1=log_close[1:(length-5)]
test_all=log_close

#HMM segmentation
##fit a HMM model
close_x=xts(data$close, order.by=strptime(data$date, format="%Y-%m-%d", tz=""))
close_x=na.omit(cbind(close_x,diff(log(close_x))))
colnames(close_x)=c('close','return')

ATRindicator=ATR(data[,2:5],n=14)
ATR=ATRindicator[,2]
ModelData=data.frame(log_close,ATR[-1])
ModelData=ModelData[-c(1:13),]
colnames(ModelData)=c("LogReturns","ATR")

##show fitted HMM states
model_data=data.frame(close_x[1:length])
log_close_test=log_close[1:length]
HMM = depmix(list(log_close_test~1), data = model_data, nstates = 3,family=list(gaussian()))
HMMfit=fit(HMM)
HMMpost=posterior(HMMfit)
plot.xts(close_x[1:length]$close,col='white')
points(close_x[1:length]$close[HMMpost$state==1],col='black')
points(close_x[1:length]$close[HMMpost$state==2],col='red')
points(close_x[1:length]$close[HMMpost$state==3],col='green')
#points(close_x[1:length]$close[HMMpost$state==4],col='blue')

##transfer HMM states to dummy variables
log_close_x_pre=log_close[(length-5):length]
model_data=data.frame(close_x)
model_data_test=data.frame(close_x[(length-5):length])

close_x_test=close_x[1:(length-5)]
log_close_x_test=log_close[1:(length-5)]
model_data=data.frame(close_x_test)
log_close_x_pre=log_close[(length-5):length]
model_data_test=data.frame(close_x[(length-5):length])

HMM = depmix(list(log_close_x_test~1), data = model_data, nstates = 3, family=list(gaussian()))
HMMfit=fit(HMM)
HMMpost=posterior(HMMfit)
HMMfor=forwardbackward(setpars(depmix(list(log_close_x_pre~1),data=model_data_test,nstates = 3,family = list(gaussian())),getpars(HMMfit)))
pred_stat=argmax(HMMfor$alpha)
pred_stat=pred_stat[-1]
HMMpost_state=c(HMMpost$state,pred_stat)

#HMM + gjrGARCH model to fit close price
dummy1=(HMMpost_state==1)
dummy2=(HMMpost_state==2)
dummy3=(HMMpost_state==3)
#dummy4=(HMMpost_state==4)
dummy=as.matrix(cbind(as.numeric(dummy1),as.numeric(dummy2),as.numeric(dummy3)))
test_dummy=dummy
gspec=ugarchspec(variance.model=list(model='gjrGARCH',garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0),include.mean = TRUE,external.regressors=test_dummy[,1:2]), distribution.model = "sstd")
garch_fit=ugarchfit(spec=gspec,data=test_all,out.sample = 5)
garch_forcast=log_close[length-5]
for (i in (1:5))
{
  garch_forcast=c(garch_forcast,garch_fit@model$pars[1,1]+garch_fit@model$pars[6,1]*test_dummy[length+i-5,1]+garch_fit@model$pars[7,1]*test_dummy[length+i-5,2])
}
errors=as.numeric()
pred_price=as.numeric(close[length+1-5])
for (i in (1:5))
{
  pred=pred_price[length(pred_price)]*exp(garch_forcast[1+i])
  pred_price=c(pred_price,pred)
  errors=c(errors,(close[length+1+i-5]-pred)/close[length+1+i-5])
}
errors

compute_errors<-function(errors)
{
  avg_error=as.numeric()
  exp_error=as.numeric()
  for (i in (1:3))
  {
    total_avg=0
    total_exp=0
    for (j in (1:3))
    {
      total_avg=total_avg+abs(errors[i+j-1])
      total_exp=total_exp+abs(errors[i+j-1]*j)
    }
    avg_error=c(avg_error,total_avg/3)
    exp_error=c(exp_error,total_exp/6)
  } 
  return(cbind(avg_error,exp_error))
}

ae_error=compute_errors(errors)
ae_error

#HMM + VAR model to fit open/high/low
frac_change=(close-open)/open
frac_high=(high-open)/open
frac_low=(open-low)/open
frac_vec=as.matrix(cbind(frac_change,frac_high,frac_low))
frac_vec_test=frac_vec[2:(length-4),]

var_dummy=test_dummy[1:(length-5),]
var_dummy_test=test_dummy[(length-4):length,]
colnames(var_dummy)=c("exo1","exo2","exo3")
var_fit=VAR(frac_vec_test,exogen = var_dummy[,1:2])
var_pred=predict(var_fit,n.ahead=5,dumvar=var_dummy_test[,1:2])

open_pred=pred_price[-1]/(1+var_pred$fcst$frac_change[,1])
high_pred=open_pred*(1+var_pred$fcst$frac_high[,1])
low_pred=open_pred*(1-var_pred$fcst$frac_low[,1])
errors_open=(open[(length+1-4):(length+1)]-open_pred)/open[(length+1-4):(length+1)]
errors_high=(high[(length+1-4):(length+1)]-high_pred)/high[(length+1-4):(length+1)]
errors_low=(low[(length+1-4):(length+1)]-low_pred)/low[(length+1-4):(length+1)]
errors_open[-1:-2]
errors_high[-1:-2]
errors_low[-1:-2]
ae_errors_open=compute_errors(errors_open)
ae_errors_high=compute_errors(errors_high)
ae_errors_low=compute_errors(errors_low)
ae_errors_open
ae_errors_high
ae_errors_low

