rm (list=ls())
library(survival)
library(randomForestSRC)
library(dplyr)
library(tidyr)
library(janitor)

setwd("C:/Users/Jiwei.He/Documents/tech application/survival project")
data<-read.csv("survival_data.csv", header=TRUE)
data<-as_tibble(data)

data2<-data%>%filter(exclusion.criteria.==0 & !is.na(PatientID))

#create status column
data3<-data2%>%mutate(status=case_when(Local_Failure. %in% c(0, 5)~0,
                                Local_Failure. %in% c(1, 4)~1))

#convert categorical variables to dummy variables for regression model
data4<-data3%>%mutate(DivType=ifelse(DivType %in% c(6, 7, 9, 10, NA), 99, DivType))%>%
  mutate(Hydro=na_if(Hydro, 9), LVI=na_if(LVI, 9), Margin=na_if(Margin, 99), Node=na_if(Node, 9),
         DivType=na_if(DivType, 99))%>%
  mutate(hydro.indicator=C(as.factor(Hydro), treatment),
         type.indicator=C(as.factor(DivType), treatment),
         neoadj.indicator=C(as.factor(Neoadj_Chemo), treatment))

####################################################################
#model 1: include all significant factors from univariate analysis #
####################################################################
m1<-coxph(Surv(TTE.Local.Failure, status)~
                    Node+Margin+pT3plus.+LVI+hydro.indicator+type.indicator+neoadj.indicator+HistologyCode+Age.at.Surg, data=data4)
summary(m1)

#################################
#Model 2, Remove Type.indicator #
#################################

m2<-coxph(Surv(TTE.Local.Failure, status)~
                    Node+Margin+pT3plus.+LVI+hydro.indicator+neoadj.indicator+HistologyCode+Age.at.Surg, data=data4)
summary(m2)

#######################################
#Model 3, Further remove LVI and neoadj
#######################################

m3<-coxph(Surv(TTE.Local.Failure, status)~
            Node+Margin+pT3plus.+hydro.indicator+HistologyCode+Age.at.Surg, data=data4)
summary(m3)

###############
#RPA approach #
###############

#combine Hydro 0 and 1 (dichotomize at bilateral Hydro)
#convert to categorical variables for rpart()
data5<-data4%>%
  mutate(Hydro.binary=case_when(Hydro %in% c(0,1)~0, Hydro==2~1))%>%
  mutate(Margin=as.factor(Margin), Node=as.factor(Node), pT3plus.=as.factor(pT3plus.),
         Hydro.binary=as.factor(Hydro.binary), HistologyCode=as.factor(HistologyCode),
         TotalNodes.dich=ifelse(TotalNodes>=10, 1, 0))

library(rpart)

r1<-rpart(Surv(TTE.Local.Failure, status)~Node+Margin+pT3plus.+Hydro.binary+HistologyCode+Age.at.Surg, data=data5)
print(r1$cptable)
plot(r1, uniform = TRUE, margin = 0.1, branch = 0.5,compress = TRUE)
text(r1, use.n=T, cex=1)

#include total nodes
r1<-rpart(Surv(TTE.Local.Failure, status)~TotalNodes+Node+Margin+pT3plus.+Hydro.binary+HistologyCode+Age.at.Surg, data=data5)
plot(r1, uniform = TRUE, margin = 0.1, branch = 0.5,compress = TRUE)
text(r1, use.n=T, cex=1)

#prune at the leftmost value for which the mean lies below the horizontal line
printcp(r1)
plotcp(r1)

## plot unprune tree ##
grouping<-r1$where
newfit<-survfit(Surv(TTE.Local.Failure, status)~factor(grouping), data=data5)
plot(newfit, mark.time=F, col=2: 7)
title(main="Kaplan Meier curves for unpruned tree", xlab='Time to event',ylab='Probablity of freedom from local failure')
legend(170,.37, legend=paste('node',1: 5), lty=1, col=2: 7)

## pruned tree ##
cp <- r1$cptable[3,"CP"]
r1<- prune(r1, cp = cp)
printcp(r1)
plot(r1, uniform = TRUE, margin = 0.1, branch = 0.5,compress = TRUE)
text(r1, use.n=T)

grouping<-r1$where
newfit<-survfit(Surv(TTE.Local.Failure, status)~factor(grouping), data=data5)
plot(newfit, mark.time=F, col=2:4)
title(main="Kaplan Meier curves for pruned tree", xlab='Time to event',ylab='Probablity of freedom from local failure')
legend(160,.37, legend=paste('node',1: 3), lty=1, col=2:4)

## explortatory analysis ##

###################
## random forest ##
###################
set.seed(123)
#mtry: Number of variables randomly selected as candidates for splitting a node
model_rf<- rfsrc(Surv(TTE.Local.Failure, status)~TotalNodes+Node+Margin+pT3plus.+Hydro.binary+HistologyCode+Age.at.Surg, 
                 data=data5, ntree=5000, mtry=2, importance=TRUE)
model_rf
plot(model_rf)

vimp(model_rf, importance="permute")
vimp(model_rf, importance="permute")$importance

var.select(model_rf, method="md")
var.select(model_rf, method="vh")

find.interaction(model_rf, nvar=7, method="vimp")

#######################
## gradient boosting ##
#######################

library(gbm)

g1<-gbm(Surv(TTE.Local.Failure, status)~Margin+pT3plus.+Hydro.binary+HistologyCode+Age.at.Surg+TotalNodes+Node, 
        data=data5)

summary(g1)

best.iter <- gbm.perf(g1, method = "OOB")
print(best.iter)


