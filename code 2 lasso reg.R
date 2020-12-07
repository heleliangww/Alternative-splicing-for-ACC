##A demo for pridictive model 
## Ref: Prognostic Alternative mRNA Splicing in Adrenocortical Carcinoma 
##Copyright Weiwei Liang

## Start clean
rm(list = ls())

## Make sure current working directory 
getwd()

###
library("glmnet") 
library("survival")

###load data
###data was availabe in the sup. material of the reference
rt_AA<-read.table("train.txt",header=T,row.names = 1)
rt_AA$OS=rt_AA$OS/365 

##Lasso analysis
set.seed(20)    
x=data.matrix(rt_AA[,c(3:ncol(rt_AA))]) 
y=data.matrix(Surv(rt_AA$OS,rt_AA$Event)) 
fit=glmnet(x, y, family = "cox", maxit = 5000) 
plot(fit, xvar = "lambda", label = TRUE)  ##figure 4a

cvfit = cv.glmnet(x, y, family="cox", maxit = 50000) 
plot(cvfit) ##figure 4b
abline(v = log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")

#call result
fit
c(cvfit$lambda.min,cvfit$lambda.1se)

###coef
coef = coef(fit, s = cvfit$lambda.min) 
index = which(coef != 0) 
actCoef = coef[index] 
lassoGene = row.names(coef)[index] 
geneCoef = cbind(Gene=lassoGene,Coef=actCoef) 
geneCoef  ##table 1 
write.table(geneCoef,"./gene-coef.txt",row.names=TRUE,sep="\t",quote=FALSE)

##calculate riskScore
FinalGeneExp = rt_AA[,lassoGene] 
myFun = function(x){crossprod(as.numeric(x),actCoef)} 
riskScore = apply(FinalGeneExp,1,myFun) 
outCol = c("OS", "Event", lassoGene) 
output = c("OS", "Event")
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low")) 
dat = cbind(rt_AA[,outCol], riskScore=as.vector(riskScore), risk)
dat_cli =cbind(rt_AA[,output], riskScore=as.vector(riskScore), risk)
write.table(dat,"./risk-score_train.txt",row.names=TRUE,sep="\t",quote=FALSE)

#plot
#install.packages("ggpubr") 
library(ggpubr)   

p <- ggboxplot(dat, x = "Event", y = "riskScore", 
               color = "Event", palette = "jco", 
               add = "jitter") 
p <- p + stat_compare_means()   #  Add p-value 
p   ##figure 5a

#calculate riskSCore for testing data
#load data
testingdata<-read.table("test.txt",header=T,row.names = 1)
testingdata$OS=testingdata$OS/365 
FinalGeneExp_testing = testingdata[,lassoGene] 
myFun = function(x){crossprod(as.numeric(x),actCoef)} 
riskScore = apply(FinalGeneExp_testing,1,myFun) 
outCol = c("OS", "Event", lassoGene) 
output = c("OS", "Event")
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low")) 
dat_test = cbind(testingdata[,outCol], riskScore=as.vector(riskScore), risk)
dat_cli_test =cbind(testingdata[,output], riskScore=as.vector(riskScore), risk)
write.table(dat_test,"./risk-score_test.txt",row.names=TRUE,sep="\t",quote=FALSE)

##plot
p1 <- ggboxplot(dat_test, x = "Event", y = "riskScore", 
               color = "Event", palette = "jco", 
               add = "jitter") 
p1 <- p1 + stat_compare_means()   #  Add p-value 
p1  ##figure 5b

##ROC curve
#install.packages("ROCR") 
library(ROCR)   
library(glmnet) 
library(caret)
library(pROC)

pred <- prediction(dat$riskScore, dat$Event) 
perf <- performance(pred,"tpr","fpr") 
performance(pred,"auc")   # shows calculated AUC for model 
plot(perf,colorize=FALSE, col="red")   #绘制ROC曲线 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )   #figure 6

roc1<-roc(dat$Event,dat$riskScore);roc1

pred1 <- prediction(dat_test$riskScore, dat_test$Event) 
perf1 <- performance(pred1,"tpr","fpr") 
performance(pred1,"auc")   # shows calculated AUC for model 
plot(perf1,colorize=FALSE, col="red")   
lines(c(0,1),c(0,1),col = "gray", lty = 4 )  #figure 6b

roc2<-roc(dat_test$Event,dat_test$riskScore);roc2

##KM curve 
library("survival")
library("survminer")
dat$OS=dat$OS*365
dat_test$OS=dat_test$OS*365

fit <- survfit(Surv(OS,Event) ~ risk, data = dat_test)
#firuge 7b
ggsurvplot(fit,
           conf.int = TRUE,
           linetype = "strata", 
           legend.title = "", 
           # surv.median.line = "hv", 
           ggtheme = theme_bw(), 
           palette = c("#E7B800", "#2E9FDF"),
           pval = TRUE,
           risk.table = TRUE,
           xlim = c(0,2000)  ) 

fit2 <- survfit(Surv(OS,Event) ~ risk, data = dat)
#firuge 7a
ggsurvplot(fit2,
           conf.int = TRUE,
           linetype = "strata", 
           legend.title = "", 
           # surv.median.line = "hv",
           ggtheme = theme_bw(), 
           palette = c("#E7B800", "#2E9FDF"),
           pval = TRUE,
           risk.table = TRUE,
           xlim = c(0,5000)  ) 


###time dependent ROC curve 1-3-5-years

library(survivalROC) # 载入程序包

dev<-dat
nobs<- NROW(dev) 
cutoff1<-365 
cutoff2<- 1095 
cutoff3<- 1825

SROC3= survivalROC(Stime = dev$OS, status = dev$Event, marker = dev$riskScore, predict.time =cutoff2, method= "KM" )
cut.op3= SROC3$cut.values[which.max(SROC3$TP-SROC3$FP)]  
cut.op3 

plot(SROC3$FP,SROC3$TP, type="l", xlim=c(0,1), ylim=c(0,1),  
     xlab = paste( "FP","\n", "AUC = ",round(SROC3$AUC,3)),
     ylab = "TP", col="blue")
abline(0,1)
legend("bottomright",c("3-year OS"),col="blue",lty=c(1,1))

SROC1= survivalROC(Stime = dev$OS, status = dev$Event, marker = dev$riskScore, predict.time =cutoff1, method= "KM" ) 
cut.op1= SROC1$cut.values[which.max(SROC1$TP-SROC1$FP)] 
cut.op1 

plot(SROC1$FP,SROC1$TP, type="l", xlim=c(0,1), ylim=c(0,1),  
     xlab = paste( "FP","\n", "AUC = ",round(SROC1$AUC,3)),
     ylab = "TP", col="red")
abline(0,1)
legend("bottomright",c("1-year OS"),col="blue",lty=c(1,1))

SROC5= survivalROC(Stime = dev$OS, status = dev$Event, marker = dev$riskScore, predict.time =cutoff3, method= "KM" ) #构建生存函数
cut.op5= SROC5$cut.values[which.max(SROC5$TP-SROC5$FP)]  
cut.op5 

plot(SROC5$FP,SROC5$TP, type="l", xlim=c(0,1), ylim=c(0,1),  
     xlab = paste( "FP","\n", "AUC = ",round(SROC5$AUC,3)),
     ylab = "TP", col="red")
abline(0,1)
legend("bottomright",c("5-year OS"),col="blue",lty=c(1,1))


