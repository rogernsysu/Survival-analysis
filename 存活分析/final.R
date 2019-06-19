######packages######
{
library(KMsurv)
library(survival)
library(dplyr)
library(ggplot2)
library(plotly)
library(GGally)
library(survminer)
}
###################
clinic=read.csv("clinica.csv")
clinic$delta=1
clinic$delta[which(clinic$Alive=="Yes")]=0
#######Build model#######
coxph(Surv(Survival_length,delta)~factor(PCR)+age,data=clinic)->modell
modell %>% summary
kmmodel=survfit(Surv(Survival_length,delta)~factor(PCR),data=clinic)
#######Visualization#######
plot(kmmodel,col=c("blue","green","gold","red"),xlab="Survival_length",ylab="K-M survival function",main="K-M Curves")
legend("bottomright",legend=c("0級","I類(RCB<=1.36)","II類(1.36<RCB<=3.28)","III類(RCB>3.28)"),lwd=2,col=c("blue","green","gold","red"),text.col = c("blue","green","gold","red"))
ggsurvplot(kmmodel,clinic,title="K-M curves of different RCB indicators",legend.title="Group",legend.labs=c("0級","I類(RCB<=1.36)","II類(1.36<RCB<=3.28)","III類(RCB>3.28)"))
ggsurv(kmmodel,main="K-M curves") %>% ggplotly
#######log-rank-test#####
survdiff(Surv(Survival_length,delta)~factor(PCR),data=clinic)
cox.zph(modell)
summary(clinic)
#######residuals#######
resd=residuals(modell,"martingale")
clinic$delta-resd->rom
rom = sort(rom)
names(clinic)
fit = survfit(Surv(rom,clinic$delta)~1,type='fh2')
fit$surv
plot(resd)
plot(rom,1-fit$surv)
fit$surv
