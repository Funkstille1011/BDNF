################################################################################
## Exploratory analysis of BDNF response to training 
## in the Phys-Stroke data
## Date: 03.07.2023
## Author: Konrad Neumann, Torsten Rackoll
################################################################################

## packages
library(mice)
library(openxlsx)
library(tidyverse)


## Helper function to create variance-covariance matrix
Rubin <- function(CO, VAR, NAMES=NULL){
  if(any(is.na(CO[1,]))){
    ind <- which(is.na(CO[1,]))
    CO <- CO[,-ind]
    VAR <- VAR[,,-ind]
  }
  m <- dim(VAR)[3]
  co <- apply(CO, MARGIN=1, FUN=mean, na.rm=TRUE)
  B <- (CO-co) %*% t(CO-co)/(m-1)
  
  Var <- VAR[,,1]
  for(j in 2:m) Var <- Var+VAR[,,j]
  Var <- Var/m
  
  var <- Var + (1+1/m)*B
  
  if(!is.null(NAMES)){
    names(co) <- NAMES
    colnames(var) <- NAMES
    rownames(var) <- NAMES
  }
  
  return(list(coef=co,var=var))
}

####################### data #########
## load data
load("Data.RData")



####################### imputation ###########
Dat.imp <- data.frame(
  Geschlecht=Master$Geschlecht,
  Alter=Master$Bas_Alter,
  Treatment=Master$Treatment,
  DBSD1=Master$DBSD1,
  Center=Master$Center,
  Smoking=Master$Bas_RF_Niko_bis,
  NIHSS = Master$Bas_NIHSS_su,
  
  bdnf.baseline=Master$bdnf.baseline_visite_arm_1,
  bdnf.v1=Master$bdnf.v1_visite_arm_1,
  bdnf.v2=Master$bdnf.v2_visite_arm_1,
  bdnf.v3=Master$bdnf.v3_visite_arm_1,
  
  CESD = Master$Bas_CESD_Sum,
  
  thrombozyten.baseline=Master$thrombozyten.baseline_visite_arm_1, 
  thrombozyten.v1=Master$thrombozyten.v1_visite_arm_1, 
  thrombozyten.v2=Master$thrombozyten.v2_visite_arm_1, 
  thrombozyten.v3=Master$thrombozyten.v3_visite_arm_1
  
)

m <- 1000 # Number of imputation samples.
ini <- mice(Dat.imp, maxit = 0)
meth <- ini$meth
meth
pred <- ini$pred

pred[,] <- 1
pred[,"DBSD1"] <- 0
diag(pred) <- 0

MI <- mice(Dat.imp, method=meth, predictorMatrix=pred, m=m, seed=7545713, maxit=5)

# store number of iterations
m <- MI$m

# prepare matrices and arrays for helper Function Rubin
VAR0 <- array(dim=c(4,4,m))
CO0 <- matrix(as.numeric(rep(NA,m*4)),ncol=m, nrow=4)# size of array needs to be adjusted for number of additional confounding factors

VAR1 <- array(dim=c(5,5,m))
CO1 <- matrix(as.numeric(rep(NA,m*5)),ncol=m, nrow=5)

VAR2 <- array(dim=c(5,5,m))
CO2 <- matrix(as.numeric(rep(NA,m*5)),ncol=m, nrow=5)

VAR3 <- array(dim=c(5,5,m))
CO3 <- matrix(as.numeric(rep(NA,m*5)),ncol=m, nrow=5)

VAR4 <- array(dim=c(5,5,m))
CO4 <- matrix(as.numeric(rep(NA,m*5)),ncol=m, nrow=5)

VAR5 <- array(dim=c(5,5,m))
CO5 <- matrix(as.numeric(rep(NA,m*5)),ncol=m, nrow=5)

VAR6 <- array(dim=c(5,5,m))
CO6 <- matrix(as.numeric(rep(NA,m*5)),ncol=m, nrow=5)

VAR7 <- array(dim=c(5,5,m))
CO7 <- matrix(as.numeric(rep(NA,m*5)),ncol=m, nrow=5)

VAR8 <- array(dim=c(5,5,m))
CO8 <- matrix(as.numeric(rep(NA,m*5)),ncol=m, nrow=5)

VAR9 <- array(dim=c(5,5,m))
CO9 <- matrix(as.numeric(rep(NA,m*5)),ncol=m, nrow=5)


# code loops m-times over the data, creates a df (s) and fills NA values with values from MI.
for(i in 1:m){
  S <- complete(MI, action=i)
  S$Center1 <- ifelse(S$Center=="Rehab",1,ifelse(!(S$Center =="Rehab"),0,NA))
  S$Center2 <- ifelse(S$Center== "Geriatric",1,ifelse(!(S$Center== "Geriatric"),0,NA))
  S$Center3 <- ifelse(S$Center == "Early-rehab",1,ifelse(!(S$Center == "Early-rehab"),0,NA))
  S$Center1 <- S$Center1 - mean(S$Center1, na.rm=TRUE)
  S$Center2 <- S$Center2 - mean(S$Center2, na.rm=TRUE)
  S$Center3 <- S$Center3 - mean(S$Center3, na.rm=TRUE)
  
  # sets bdnf values from visits
  S$X0 <- S$bdnf.baseline
  S$X1 <- S$bdnf.v1
  S$X2 <- S$bdnf.v2
  S$X3 <- S$bdnf.v3
  
  # calculates differences between measurement time points
  S$X4 <- S$X1 - S$X0
  S$X5 <- S$X2 - S$X0
  S$X6 <- S$X3 - S$X0
  S$X7 <- S$X2 - S$X1
  S$X8 <- S$X3 - S$X1
  S$X9 <- S$X3 - S$X2
  
  # sets thrombocyztes values from visits
  S$T0 <- S$thrombozyten.baseline
  S$T1 <- S$thrombozyten.v1
  S$T2 <- S$thrombozyten.v2
  S$T3 <- S$thrombozyten.v3
  
  # calculates differences between measurement time points of thrombocytes
  S$T4 <- S$T1 - S$T0
  S$T5 <- S$T2 - S$T0
  S$T6 <- S$T3 - S$T0
  S$T7 <- S$T2 - S$T1
  S$T8 <- S$T3 - S$T1
  S$T9 <- S$T3 - S$T2
  
  # centers baseline variables of BDNF and thrombozytes
  S$BL <- S$X0 - mean(S$X0, na.rm=TRUE)
  S$T.BL <- S$T0 - mean(S$T0, na.rm=TRUE)
  S$T.V1 <- S$T1 - mean(S$T1, na.rm=TRUE)
  S$T.V2 <- S$T2 - mean(S$T2, na.rm=TRUE)
  S$T.V3 <- S$T3 - mean(S$T3, na.rm=TRUE)
  S$T4.c <- S$T4 - mean(S$T4, na.rm=TRUE)
  S$T5.c <- S$T5 - mean(S$T5, na.rm=TRUE)
  S$T6.c <- S$T6 - mean(S$T6, na.rm=TRUE)
  S$T7.c <- S$T7 - mean(S$T7, na.rm=TRUE)
  S$T8.c <- S$T8 - mean(S$T8, na.rm=TRUE)
  S$T9.c <- S$T9 - mean(S$T9, na.rm=TRUE)
  S$NIHSS.center <- S$NIHSS - mean(S$NIHSS, na.rm = TRUE)
  S$Smoking.center <- S$Smoking - mean(S$Smoking, na.rm=T)
  S$CESD.center <- S$CESD - mean(S$CESD, na.rm = T)
  S$Geschlecht.center <- as.numeric(as.character(S$Geschlecht)) - mean(as.numeric(as.character(S$Geschlecht)), na.rm = T)
  S$Alter.center <- S$Alter - mean(S$Alter, na.rm = T)
  ## need to center confounding factors 
  # factors that are confounders need to be centered as well
  
  ################################################
  # calculates linear models for all time differences
  M0 <- try(lm(X0~DBSD1+Alter.center+Geschlecht.center+NIHSS.center+T.BL+CESD.center+Smoking.center+Center1+Center2,data=S), silent=TRUE) #Zentrum in 3 dummy variablen plus Zentrierung: Reha, Früh-Reha und Geriatrie
  M1 <- try(lm(X1~DBSD1+BL+Alter.center+Geschlecht.center+NIHSS.center+T.V1+CESD.center+Smoking.center+Center1+Center2,data=S), silent=TRUE)
  M2 <- try(lm(X2~DBSD1+BL+Alter.center+Geschlecht.center+NIHSS.center+T.V2+CESD.center+Smoking.center+Center1+Center2,data=S), silent=TRUE)
  M3 <- try(lm(X3~DBSD1+BL+Alter.center+Geschlecht.center+NIHSS.center+T.V3+CESD.center+Smoking.center+Center1+Center2,data=S), silent=TRUE)
  M4 <- try(lm(X4~DBSD1+BL+Alter.center+Geschlecht.center+NIHSS.center+T4.c+CESD.center+Smoking.center+Center1+Center2,data=S), silent=TRUE)
  M5 <- try(lm(X5~DBSD1+BL+Alter.center+Geschlecht.center+NIHSS.center+T5.c+CESD.center+Smoking.center+Center1+Center2,data=S), silent=TRUE)
  M6 <- try(lm(X6~DBSD1+BL+Alter.center+Geschlecht.center+NIHSS.center+T6.c+CESD.center+Smoking.center+Center1+Center2,data=S), silent=TRUE)
  M7 <- try(lm(X7~DBSD1+BL+Alter.center+Geschlecht.center+NIHSS.center+T7.c+CESD.center+Smoking.center+Center1+Center2,data=S), silent=TRUE)
  M8 <- try(lm(X8~DBSD1+BL+Alter.center+Geschlecht.center+NIHSS.center+T8.c+CESD.center+Smoking.center+Center1+Center2,data=S), silent=TRUE)
  M9 <- try(lm(X9~DBSD1+BL+Alter.center+Geschlecht.center+NIHSS.center+T9.c+CESD.center+Smoking.center+Center1+Center2,data=S), silent=TRUE)
  # add an additional model with additional covariants
  
  
  # if linear models do not create an error, results will be written 
  # in variance-covaricance matrix
  if(class(M0)[1]!="try-error"){		
    CO0[,i] <- coef(M0)
    VAR0[,,i] <- vcov(M0)
  }
  
  if(class(M1)[1]!="try-error"){		
    CO1[,i] <- coef(M1)
    VAR1[,,i] <- vcov(M1)
  }
  
  if(class(M2)[1]!="try-error"){		
    CO2[,i] <- coef(M2)
    VAR2[,,i] <- vcov(M2)
  }
  
  if(class(M3)[1]!="try-error"){		
    CO3[,i] <- coef(M3)
    VAR3[,,i] <- vcov(M3)
  }
  
  if(class(M4)[1]!="try-error"){		
    CO4[,i] <- coef(M4)
    VAR4[,,i] <- vcov(M4)
  }
  
  if(class(M5)[1]!="try-error"){		
    CO5[,i] <- coef(M5)
    VAR5[,,i] <- vcov(M5)
  }
  
  if(class(M6)[1]!="try-error"){		
    CO6[,i] <- coef(M6)
    VAR6[,,i] <- vcov(M6)
  }
  
  if(class(M7)[1]!="try-error"){		
    CO7[,i] <- coef(M7)
    VAR7[,,i] <- vcov(M7)
  }
  
  if(class(M8)[1]!="try-error"){		
    CO8[,i] <- coef(M8)
    VAR8[,,i] <- vcov(M8)
  }
  
  if(class(M9)[1]!="try-error"){		
    CO9[,i] <- coef(M9)
    VAR9[,,i] <- vcov(M9)
    
  }
}


# 9 lists of imputed values are generated
RU0 <- Rubin(CO0,VAR0, NAMES=names(coef(M0)))
Var0 <- RU0$var
B0 <- RU0$coef

RU1 <- Rubin(CO1,VAR1, NAMES=names(coef(M1)))
Var1 <- RU1$var
B1 <- RU1$coef

RU2 <- Rubin(CO2,VAR2, NAMES=names(coef(M2)))
Var2 <- RU2$var
B2 <- RU2$coef

RU3 <- Rubin(CO3,VAR3, NAMES=names(coef(M3)))
Var3 <- RU3$var
B3 <- RU3$coef

RU4 <- Rubin(CO4,VAR4, NAMES=names(coef(M4)))
Var4 <- RU4$var
B4 <- RU4$coef

RU5 <- Rubin(CO5,VAR5, NAMES=names(coef(M5)))
Var5 <- RU5$var
B5 <- RU5$coef

RU6 <- Rubin(CO6,VAR6, NAMES=names(coef(M6)))
Var6 <- RU6$var
B6 <- RU6$coef

RU7 <- Rubin(CO7,VAR7, NAMES=names(coef(M7)))
Var7 <- RU7$var
B7 <- RU7$coef

RU8 <- Rubin(CO8,VAR8, NAMES=names(coef(M8)))
Var8 <- RU8$var
B8 <- RU8$coef

RU9 <- Rubin(CO9,VAR9, NAMES=names(coef(M9)))
Var9 <- RU9$var
B9 <- RU9$coef

# add one step for M10

# helper variables for later comparisons are created

# helper variables for later comparisons are created
n <- nrow(Master)
v1 <- c(1,0,0,0,0,0,0,0,0,0,0)
v2 <- c(1,1,0,0,0,0,0,0,0,0,0)
v3 <- c(1,0,1,0,0,0,0,0,0,0,0)

# calculates means and sd for different timepoints and the changes
erg0 <- c(
  B0%*%v1,
  B0%*%v1-sqrt(t(v1)%*%Var0%*%v1)*qt(0.975,n-11),
  B0%*%v1+sqrt(t(v1)%*%Var0%*%v1)*qt(0.975,n-11),
  
  B0%*%v2,
  B0%*%v2-sqrt(t(v2)%*%Var0%*%v2)*qt(0.975,n-11),
  B0%*%v2+sqrt(t(v2)%*%Var0%*%v2)*qt(0.975,n-11),
  
  B0%*%v3,
  B0%*%v3-sqrt(t(v3)%*%Var0%*%v3)*qt(0.975,n-11),
  B0%*%v3+sqrt(t(v3)%*%Var0%*%v3)*qt(0.975,n-11),
  
  1- pf((t(B0[c(2,3)])%*%solve(Var0[c(2,3),c(2,3)])%*%B0[c(2,3)])/2, df1=2,df2=n-11),
  
  2*pt(-abs(B0%*%(v2-v1)/sqrt(t(v2-v1)%*%Var0%*%(v2-v1))),df=n-11),
  2*pt(-abs(B0%*%(v3-v1)/sqrt(t(v3-v1)%*%Var0%*%(v3-v1))),df=n-11),
  2*pt(-abs(B0%*%(v3-v2)/sqrt(t(v3-v2)%*%Var0%*%(v3-v2))),df=n-11)
)


w1 <- c(1,0,0,0,0,0,0,0,0,0,0,0)
w2 <- c(1,1,0,0,0,0,0,0,0,0,0,0)
w3 <- c(1,0,1,0,0,0,0,0,0,0,0,0)

erg1 <- c(
  B1%*%w1,
  B1%*%w1-sqrt(t(w1)%*%Var1%*%w1)*qt(0.975,n-12),
  B1%*%w1+sqrt(t(w1)%*%Var1%*%w1)*qt(0.975,n-12),
  
  B1%*%w2,
  B1%*%w2-sqrt(t(w2)%*%Var1%*%w2)*qt(0.975,n-12),
  B1%*%w2+sqrt(t(w2)%*%Var1%*%w2)*qt(0.975,n-12),
  
  B1%*%w3,
  B1%*%w3-sqrt(t(w3)%*%Var1%*%w3)*qt(0.975,n-12),
  B1%*%w3+sqrt(t(w3)%*%Var1%*%w3)*qt(0.975,n-12),
  
  1- pf((t(B1[c(2,3)])%*%solve(Var1[c(2,3),c(2,3)])%*%B1[c(2,3)])/2, df1=2,df2=n-12),
  
  2*pt(-abs(B1%*%(w2-w1)/sqrt(t(w2-w1)%*%Var1%*%(w2-w1))),df=n-12),
  2*pt(-abs(B1%*%(w3-w1)/sqrt(t(w3-w1)%*%Var1%*%(w3-w1))),df=n-12),
  2*pt(-abs(B1%*%(w3-w2)/sqrt(t(w3-w2)%*%Var1%*%(w3-w2))),df=n-12)
)

erg2 <- c(
  B2%*%w1,
  B2%*%w1-sqrt(t(w1)%*%Var2%*%w1)*qt(0.975,n-12),
  B2%*%w1+sqrt(t(w1)%*%Var2%*%w1)*qt(0.975,n-12),
  
  B2%*%w2,
  B2%*%w2-sqrt(t(w2)%*%Var2%*%w2)*qt(0.975,n-12),
  B2%*%w2+sqrt(t(w2)%*%Var2%*%w2)*qt(0.975,n-12),
  
  B2%*%w3,
  B2%*%w3-sqrt(t(w3)%*%Var2%*%w3)*qt(0.975,n-12),
  B2%*%w3+sqrt(t(w3)%*%Var2%*%w3)*qt(0.975,n-12),
  
  1- pf((t(B2[c(2,3)])%*%solve(Var2[c(2,3),c(2,3)])%*%B2[c(2,3)])/2, df1=2,df2=n-12),
  
  2*pt(-abs(B2%*%(w2-w1)/sqrt(t(w2-w1)%*%Var2%*%(w2-w1))),df=n-12),
  2*pt(-abs(B2%*%(w3-w1)/sqrt(t(w3-w1)%*%Var2%*%(w3-w1))),df=n-12),
  2*pt(-abs(B2%*%(w3-w2)/sqrt(t(w3-w2)%*%Var2%*%(w3-w2))),df=n-12)
)


erg3 <- c(
  B3%*%w1,
  B3%*%w1-sqrt(t(w1)%*%Var3%*%w1)*qt(0.975,n-12),
  B3%*%w1+sqrt(t(w1)%*%Var3%*%w1)*qt(0.975,n-12),
  
  B3%*%w2,
  B3%*%w2-sqrt(t(w2)%*%Var3%*%w2)*qt(0.975,n-12),
  B3%*%w2+sqrt(t(w2)%*%Var3%*%w2)*qt(0.975,n-12),
  
  B3%*%w3,
  B3%*%w3-sqrt(t(w3)%*%Var3%*%w3)*qt(0.975,n-12),
  B3%*%w3+sqrt(t(w3)%*%Var3%*%w3)*qt(0.975,n-12),
  
  1- pf((t(B3[c(2,3)])%*%solve(Var3[c(2,3),c(2,3)])%*%B3[c(2,3)])/2, df1=2,df2=n-12),
  
  2*pt(-abs(B3%*%(w2-w1)/sqrt(t(w2-w1)%*%Var3%*%(w2-w1))),df=n-12),
  2*pt(-abs(B3%*%(w3-w1)/sqrt(t(w3-w1)%*%Var3%*%(w3-w1))),df=n-12),
  2*pt(-abs(B3%*%(w3-w2)/sqrt(t(w3-w2)%*%Var3%*%(w3-w2))),df=n-12)
)

erg4 <- c(
  B4%*%w1,
  B4%*%w1-sqrt(t(w1)%*%Var4%*%w1)*qt(0.975,n-12),
  B4%*%w1+sqrt(t(w1)%*%Var4%*%w1)*qt(0.975,n-12),
  
  B4%*%w2,
  B4%*%w2-sqrt(t(w2)%*%Var4%*%w2)*qt(0.975,n-12),
  B4%*%w2+sqrt(t(w2)%*%Var4%*%w2)*qt(0.975,n-12),
  
  B4%*%w3,
  B4%*%w3-sqrt(t(w3)%*%Var4%*%w3)*qt(0.975,n-12),
  B4%*%w3+sqrt(t(w3)%*%Var4%*%w3)*qt(0.975,n-12),
  
  1- pf((t(B4[c(2,3)])%*%solve(Var4[c(2,3),c(2,3)])%*%B4[c(2,3)])/2, df1=2,df2=n-12),
  
  2*pt(-abs(B4%*%(w2-w1)/sqrt(t(w2-w1)%*%Var4%*%(w2-w1))),df=n-12),
  2*pt(-abs(B4%*%(w3-w1)/sqrt(t(w3-w1)%*%Var4%*%(w3-w1))),df=n-12),
  2*pt(-abs(B4%*%(w3-w2)/sqrt(t(w3-w2)%*%Var4%*%(w3-w2))),df=n-12)
)

erg5 <- c(
  B5%*%w1,
  B5%*%w1-sqrt(t(w1)%*%Var5%*%w1)*qt(0.975,n-12),
  B5%*%w1+sqrt(t(w1)%*%Var5%*%w1)*qt(0.975,n-12),
  
  B5%*%w2,
  B5%*%w2-sqrt(t(w2)%*%Var5%*%w2)*qt(0.975,n-12),
  B5%*%w2+sqrt(t(w2)%*%Var5%*%w2)*qt(0.975,n-12),
  
  B5%*%w3,
  B5%*%w3-sqrt(t(w3)%*%Var5%*%w3)*qt(0.975,n-12),
  B5%*%w3+sqrt(t(w3)%*%Var5%*%w3)*qt(0.975,n-12),
  
  1- pf((t(B5[c(2,3)])%*%solve(Var5[c(2,3),c(2,3)])%*%B5[c(2,3)])/2, df1=2,df2=n-12),
  
  2*pt(-abs(B5%*%(w2-w1)/sqrt(t(w2-w1)%*%Var5%*%(w2-w1))),df=n-12),
  2*pt(-abs(B5%*%(w3-w1)/sqrt(t(w3-w1)%*%Var5%*%(w3-w1))),df=n-12),
  2*pt(-abs(B5%*%(w3-w2)/sqrt(t(w3-w2)%*%Var5%*%(w3-w2))),df=n-12)
)

erg6 <- c(
  B6%*%w1,
  B6%*%w1-sqrt(t(w1)%*%Var6%*%w1)*qt(0.975,n-12),
  B6%*%w1+sqrt(t(w1)%*%Var6%*%w1)*qt(0.975,n-12),
  
  B6%*%w2,
  B6%*%w2-sqrt(t(w2)%*%Var6%*%w2)*qt(0.975,n-12),
  B6%*%w2+sqrt(t(w2)%*%Var6%*%w2)*qt(0.975,n-12),
  
  B6%*%w3,
  B6%*%w3-sqrt(t(w3)%*%Var6%*%w3)*qt(0.975,n-12),
  B6%*%w3+sqrt(t(w3)%*%Var6%*%w3)*qt(0.975,n-12),
  
  1- pf((t(B6[c(2,3)])%*%solve(Var6[c(2,3),c(2,3)])%*%B6[c(2,3)])/2, df1=2,df2=n-12),
  
  2*pt(-abs(B6%*%(w2-w1)/sqrt(t(w2-w1)%*%Var6%*%(w2-w1))),df=n-12),
  2*pt(-abs(B6%*%(w3-w1)/sqrt(t(w3-w1)%*%Var6%*%(w3-w1))),df=n-12),
  2*pt(-abs(B6%*%(w3-w2)/sqrt(t(w3-w2)%*%Var6%*%(w3-w2))),df=n-12)
)

erg7 <- c(
  B7%*%w1,
  B7%*%w1-sqrt(t(w1)%*%Var7%*%w1)*qt(0.975,n-12),
  B7%*%w1+sqrt(t(w1)%*%Var7%*%w1)*qt(0.975,n-12),
  
  B7%*%w2,
  B7%*%w2-sqrt(t(w2)%*%Var7%*%w2)*qt(0.975,n-12),
  B7%*%w2+sqrt(t(w2)%*%Var7%*%w2)*qt(0.975,n-12),
  
  B7%*%w3,
  B7%*%w3-sqrt(t(w3)%*%Var7%*%w3)*qt(0.975,n-12),
  B7%*%w3+sqrt(t(w3)%*%Var7%*%w3)*qt(0.975,n-12),
  
  1- pf((t(B7[c(2,3)])%*%solve(Var7[c(2,3),c(2,3)])%*%B7[c(2,3)])/2, df1=2,df2=n-12),
  
  2*pt(-abs(B7%*%(w2-w1)/sqrt(t(w2-w1)%*%Var7%*%(w2-w1))),df=n-12),
  2*pt(-abs(B7%*%(w3-w1)/sqrt(t(w3-w1)%*%Var7%*%(w3-w1))),df=n-12),
  2*pt(-abs(B7%*%(w3-w2)/sqrt(t(w3-w2)%*%Var7%*%(w3-w2))),df=n-12)
)

erg8 <- c(
  B8%*%w1,
  B8%*%w1-sqrt(t(w1)%*%Var8%*%w1)*qt(0.975,n-12),
  B8%*%w1+sqrt(t(w1)%*%Var8%*%w1)*qt(0.975,n-12),
  
  B8%*%w2,
  B8%*%w2-sqrt(t(w2)%*%Var8%*%w2)*qt(0.975,n-12),
  B8%*%w2+sqrt(t(w2)%*%Var8%*%w2)*qt(0.975,n-12),
  
  B8%*%w3,
  B8%*%w3-sqrt(t(w3)%*%Var8%*%w3)*qt(0.975,n-12),
  B8%*%w3+sqrt(t(w3)%*%Var8%*%w3)*qt(0.975,n-12),
  
  1- pf((t(B8[c(2,3)])%*%solve(Var8[c(2,3),c(2,3)])%*%B8[c(2,3)])/2, df1=2,df2=n-12),
  
  2*pt(-abs(B8%*%(w2-w1)/sqrt(t(w2-w1)%*%Var8%*%(w2-w1))),df=n-12),
  2*pt(-abs(B8%*%(w3-w1)/sqrt(t(w3-w1)%*%Var8%*%(w3-w1))),df=n-12),
  2*pt(-abs(B8%*%(w3-w2)/sqrt(t(w3-w2)%*%Var8%*%(w3-w2))),df=n-12)
)

erg9 <- c(
  B9%*%w1,
  B9%*%w1-sqrt(t(w1)%*%Var9%*%w1)*qt(0.975,n-12),
  B9%*%w1+sqrt(t(w1)%*%Var9%*%w1)*qt(0.975,n-12),
  
  B9%*%w2,
  B9%*%w2-sqrt(t(w2)%*%Var9%*%w2)*qt(0.975,n-12),
  B9%*%w2+sqrt(t(w2)%*%Var9%*%w2)*qt(0.975,n-12),
  
  B9%*%w3,
  B9%*%w3-sqrt(t(w3)%*%Var9%*%w3)*qt(0.975,n-12),
  B9%*%w3+sqrt(t(w3)%*%Var9%*%w3)*qt(0.975,n-12),
  
  1- pf((t(B9[c(2,3)])%*%solve(Var9[c(2,3),c(2,3)])%*%B9[c(2,3)])/2, df1=2,df2=n-12),
  
  2*pt(-abs(B9%*%(w2-w1)/sqrt(t(w2-w1)%*%Var9%*%(w2-w1))),df=n-12),
  2*pt(-abs(B9%*%(w3-w1)/sqrt(t(w3-w1)%*%Var9%*%(w3-w1))),df=n-12),
  2*pt(-abs(B9%*%(w3-w2)/sqrt(t(w3-w2)%*%Var9%*%(w3-w2))),df=n-12)
  # need to add one step for M10
)

#combines distribution functions
ERG <- rbind(
  erg0,
  erg1,
  erg2,
  erg3,
  erg4,
  erg5,
  erg6,
  erg7,
  erg8,
  erg9
)


rownames(ERG) <- c(
  "Baseline",
  "V1",
  "V2",
  "V3",
  "V1-BL",
  "V2-BL",
  "V3-BL",
  "V2-V1",
  "V3-V1",
  "V3-V2"
)

colnames(ERG) <- c(
  "<=15 Mean",
  "CI95-lower",
  "CI95-upper",
  "16-30 Mean",
  "CI95-lower",
  "CI95-upper",
  ">=31 Mean",
  "CI95-lower",
  "CI95-upper",
  "p.value global",
  "p.value <=15 vs. 16-30",
  "p.value <=15 vs. >=31",
  "p.value 16-30 vs. >=31"
)