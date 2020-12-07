rm(list = ls())
library(survival)
clinical<-read.table("./clinical_finial.txt",header=TRUE,sep="\t")
rownames(clinical)<-clinical[,1]
clinical<-clinical[,-1]

AS<-read.table("./AS.txt",header=TRUE,sep="\t")
rownames(AS)<-AS[,2]
AS<-AS[,-2]
AS<-as.matrix(AS)

intersection<-intersect(rownames(clinical),colnames(AS))
AS_a<-AS[,intersection]
AS_b<-AS[,1:9]
AS_clinical<-cbind(AS_a,AS_b)
cli<-clinical[intersection,]

y<-Surv(as.numeric(cli[,1]),as.numeric(cli[,2]))

ACCmarker <- c() #ID
ACCgene <- c() #gene
ACCpvalue<-c() 
ACCHR<-c() 
ACClow<-c() #Low 95%CI
ACChigh<-c() #High 95%CI

for(x in rownames(AS_a)){
  cox.fit_dan <- coxph(y~as.numeric(AS_a[x,]))
  coxresult<-summary(cox.fit_dan)
  pvalue=coxresult$coefficients[5]
  gene<- AS_b[x,1]
  
  ACCmarker=c(ACCmarker,x)
  ACCgene<-c(ACCgene,gene)
  ACCpvalue<-c(ACCpvalue,pvalue)
  ACCHR<-c(ACCHR,coxresult$conf.int[1])
  ACClow<-c(ACClow,coxresult$conf.int[3])
  ACChigh<-c(ACChigh,coxresult$conf.int[4])
  
}

ACC_cox_single_all<-cbind(ACCmarker,ACCgene,ACCpvalue,ACCHR,ACClow,ACChigh)
name<-c("marker id","gene symbol","P value","HR","Low 95%CI","High 95%CI")
