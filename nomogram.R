getwd() #获取当前工作路径
setwd("/Users/liangweiwei/Desktop/2020ACC")

library(survival)
library(rms)


nomo_dat<-read.table("./nomo_dat.txt",header=TRUE,sep="\t")
nomo_dat<-as.data.frame(nomo_dat)
dd<-datadist(nomo_dat)

options(datadist="dd")


f<-cph(Surv(OS,Event) ~  T.status + N.status +M.status +Path.status + riskScore, data = nomo_dat,
       x=TRUE, y=TRUE,surv=TRUE)
survival<-Survival(f)
survival1<-function(x)survival(1,x)
survival2<-function(x)survival(3,x)
nom<-nomogram(f,fun = list(survival1,survival2),
              fun.at = c(0.99,0.5,0.01),
              funlabel = c('1 year survival','3 years survival'))
plot(nom)
