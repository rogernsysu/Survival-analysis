######Packages######
{
library(splines)
library(car)
library(KMsurv)
library(survival)
library(dplyr)
library(ggplot2)
library(plotly)
library(GGally)
library(survminer)
library(leaps)
library(MASS)
library(gbm)
}
detach("package:MASS")
data(colon)
#######data cleaning#######
data(colon)
colon %>% head
na.omit(colon)->colon
colonnew = select(colon,-c(age,nodes,id,study,time,status))
colonnew1 = select(colon,c(age,nodes,id,study,time,status))
apply(colonnew,2,as.factor) -> colonnew
colonnew %>% as.data.frame -> colonnew
colon = cbind(colonnew1,colonnew)
summary(colon)
#######Build model#######
coxph(Surv(time,status)~.-id-study,data=colon) -> modell
modell
summary(modell)
library(MASS)
stepAIC(modell,direction="backward")
coxph(Surv(time,status)~ rx + nodes + obstruct +
        adhere + extent + surg + node4, data=colon, method="breslow")->coxmodel
coxmodel %>% summary
coxmodel$coefficients
plot(basehaz(coxmodeltd)$time,basehaz(coxmodeltd)$hazard)
########time-dependent data and model#######
colon %>% arrange(id,time)-> colon
mutate(colon,start=0,end=time)->colondata ###新增start跟end
filter(colondata,status==0)->nonevent
filter(colondata,status==1)->eventcolon
eventcolon$id %>% table %>% data.frame ->ids
colnames(ids)=c("id","Freq")
eventcolon[eventcolon$id %in% ids$ids[which(ids$Freq==1)],]->lonely
eventcolon[eventcolon$id %in% ids$ids[which(ids$Freq==2)],]->together
for(i in 1:393){
  together[(2*i),17]=together[(2*i)-1,18]##end to start
  together[(2*i-1),6]=0
}
rbind(nonevent,lonely,together) %>% arrange(id,etype)->colonn
colontime=colonn
#######################
colontime[-which(colontime$start>=colontime$end),] %>% arrange(id)->colontime
#######################
coxph(Surv(start,end,status)~.-time-id-study,data=colontime) -> modeltd
stepAIC(modeltd)
coxmodeltd=coxph(formula = Surv(start, end, status) ~ rx + nodes + obstruct +
        adhere + extent + surg + node4 , data = colontime)
coxmodeltd %>% summary
######confidence interval#####
confint(coxmodel)
confint(coxmodeltd)
exp(confint(coxmodel))
exp(confint(coxmodeltd))
hazmod=basehaz(coxmodel)
hazmodtd=basehaz(coxmodeltd)
AIC(coxmodel)
AIC(coxmodeltd)



#######Visualization##############
kmmodel = survfit(Surv(time,status)~rx,data=colon)#本資料主要檢驗項目
ggsurvplot(kmmodel,colontime,conf.int=T,pval=T,ylab="K-M survival function",title="K-M curves of different rx approaches",legend.title="rx")
###############################
kmmodel2 = survfit(Surv(time,status)~extent,data=colon)#本資料主要檢驗項目
ggsurvplot(kmmodel2,colon,conf.int=T,pval=T,ylab="K-M survival function",title="K-M curves of different extent situations",legend.title="extent")
###############################
kmmodel3 = survfit(Surv(time,status)~adhere,data=colon)
ggsurvplot(kmmodel3,colon,conf.int=T,pval=T,ylab="K-M survival function",title="K-M curves of different adhere situations",legend.title="adhere")
###############################
kmmodel4 = survfit(Surv(time,status)~obstruct,data=colon)
ggsurvplot(kmmodel4,colon,conf.int=T,pval=T,ylab="K-M survival function",title="K-M curves of different obstruct situations",legend.title="obstruct")
###############################
kmmodel5 = survfit(Surv(time,status)~node4,data=colon)
ggsurvplot(kmmodel5,colon,conf.int=T,pval=T,ylab="K-M survival function",title="K-M curves of different node4",legend.title="whether node4")
###############################
kmmodel6 = survfit(Surv(time,status)~surg,data=colon)
ggsurvplot(kmmodel6,colon,conf.int=T,pval=T,ylab="K-M survival function",title="K-M curves of surgery",legend.title="whether surgery")
#########comprehensive comparsion###########
merg=survfit(Surv(time,status)~ node4 + rx + extent,data=colon)
ggsurvt <- ggsurvplot(merg, fun = "event", conf.int = F,
                      ggtheme = theme_bw(),title="Different situations comprehensive comparsion")
ggsurvt$plot +theme_bw() + 
  theme (legend.position = "right")+
  facet_grid(rx ~ extent)


#########local test###########

##wald test##
beta1hat=coxmodel$coefficients[1:2]
var11=coxmodel$var[1:2,1:2]
chi=(beta1hat %>% t)%*%solve(var11)%*% beta1hat # test-statistic
1-pchisq(chi,2) #chi-square distribution with df 2
##likelihood ratio test##
fit_model.reduced=coxph(Surv(time,status)~ nodes + obstruct + adhere + extent + surg + node4 ,data=colon)
LR=2*(coxmodel$loglik[2]-fit_model.reduced$loglik[2])
1-pchisq(LR,2) #chi-square distribution with df 2
##score test##
fit0=coxph(Surv(time,status) ~ rx + nodes + obstruct + adhere + extent + surg + node4,data=colon,init=c(0,0,fit_model.reduced$coefficients),iter=0)->coxmodel
score.vector=colSums(coxph.detail(fit0)$score)
chiSC=t(score.vector[1:2])%*%fit0$var[1:2,1:2]%*%score.vector[1:2] 
# test-statistic
1-pchisq(chiSC,2) #chi-square distribution with df 2

############time dependent local test##########################

##wald test##
tdhat=coxmodeltd$coefficients[1:2]
vartd=coxmodeltd$var[1:2,1:2]
chitd=(tdhat %>% t) %*% solve(vartd) %*% tdhat # test-statistic
1-pchisq(chitd,2) #chi-square distribution with df 2
##likelihood ratio test##
td.reduced=coxph(Surv(start,time,status)~ nodes + obstruct + adhere + extent + surg + node4 ,data=colontime)
LRtd=2*(modeltd$loglik[2]-td.reduced$loglik[2])
1-pchisq(LRtd,2) #chi-square distribution with df 2
##score test##
fittd=coxph(Surv(start,end,status) ~ rx + nodes + obstruct + adhere + extent + surg + node4,data=colontime,init=c(0,0,td.reduced$coefficients),iter=0)->coxmodel
score.vector.td=colSums(coxph.detail(fittd)$score)
chiSC.td=t(score.vector.td[1:2])%*%fittd$var[1:2,1:2]%*%score.vector.td[1:2] 
# test-statistic
1-pchisq(chiSC.td,2) #chi-square distribution with df 2


#######original and td comparsion#######

cox.zph(coxmodel) -> coxz
coxz
coxz %>% ggcoxzph
ggcoxzph(cox.zph(coxmodel,transform = "log"))
cox.zph(coxmodeltd) %>% ggcoxzph
coxph(Surv(time,status)~rx + bs(nodes) + obstruct + adhere + extent + node4,data = colon) -> modelfit2
cox.zph(modelfit2)
coxmodel
coxmodel %>% summary
#######residuals#######
resd = residuals(coxmodel,"martingale")
fixedres=residuals(modelfit2,"martinagle")
resdtd = residuals(coxmodeltd,"martingale")
colon$status-resd -> rom
colontime$status-resdtd -> rtd
rom = sort(rom)
rtd = sort(rtd) 
fit = coxph(Surv(rom,colon$status)~1,data=colon,ties="efron")
tdfit=coxph(Surv(rtd,colontime$status)~1,data=colontime,ties="efron")
basehaz(fit,centered=F)->baseha
baseha
basehaz(tdfit)->basehatd
snel=qplot(x = unique(rom),y = baseha$hazard,xlab="cox-snell residuals",ylab="Estimated cumulative hazard function",main="Cox-snell residuals plot")+geom_abline(slope=1,color="red")
sneltd=qplot(x = unique(rtd),y = basehatd$hazard,xlab="cox-snell residuals",ylab="Estimated cumulative hazard function",main="Cox-snell residuals plot with time dependent data")+geom_abline(slope=1,color="red")
snel 
sneltd
qplot(colon$nodes,resd,xlab="nodes",ylab="martingale residuals",main="martingale residuals Plot")+geom_smooth(method='lm')
qplot(colontime$nodes,resdtd,xlab="nodes",ylab="martingale residuals",main="martingale residuals Plot with time-dependent data")+geom_smooth(method='lm')




