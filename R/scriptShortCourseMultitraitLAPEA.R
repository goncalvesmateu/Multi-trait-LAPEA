################################################################################
#
# Script of chapter 5
#
# Multivariate Statistical machine Learning methods for genomic prediction
#
# by Osval Antonio Montesinos López, Abelardo Montesinos López & Jose Crossa 
#
################################################################################
# Installing packages
# Github version
install.packages('devtools');
library(devtools);
install_github('covaruber/sommer')
# CRAN version
install.packages('sommer',dependencies = TRUE)
library(sommer)
################################################################################
################################################################################
# Cross-validation Single-trait models (rrBLUP) using sommer
# Datasets
data(DT_cpdata)
# Phenotypic data
DT <- DT_cpdata
head(DT)
# Marker data
GT <- GT_cpdata
GT_cpdata[1:5,1:5]
# Chromossome positions
#MP <- MP_cpdata
#head(MP_cpdata)
# Genomic relationship matrix
G <- A.mat(GT) # additive relationship matrix
colnames(G) <- rownames(G) <- rownames(DT)
G[1:5,1:5]
# Creating new dataset
data.1.2 <- DT
# Scaling phenotypes
data.1.2$Yield <- as.vector(scale(data.1.2$Yield, center = T, scale = T)) # Yield
data.1.2$FruitAver <- as.vector(scale(data.1.2$FruitAver, center = T, scale = T)) # fruit avarage
data.1.2$Firmness <- as.vector(scale(data.1.2$Firmness, center = T, scale = T)) # firmness
data.1.2$color <- as.vector(scale(data.1.2$color, center = T, scale = T)) # color
######################################################################################
# Fitting models using sommer package to compute narrow-sense heritability estimates for Yield
# No  marker information  
mod.1 = mmer(Yield ~ 1,
             random= ~ vsr(id,Gtc=unsm(1)) + Rowf,
             rcov = ~ units,
             data = data.1.2,verbose = F)
# Model summary
(summary(mod.1)$varcomp)
# Narrow-sense heritability
h2b <- (mod.1$sigma$`u:id`+mod.1$sigma$Rowf) /(mod.1$sigma$`u:id`+mod.1$sigma$Rowf+mod.1$sigma$units); h2b #  broad-sense manually
h2a <- (mod.1$sigma$`u:id`) /(mod.1$sigma$`u:id`+mod.1$sigma$units); h2a # narrow sense manually
vpredict(mod.1, h2 ~ (V1+V2) / ( V1+V2+V3) ) # broad-sense function
vpredict(mod.1, h2 ~ (V1) / ( V1+V3) ) # narrow sense function
# Genomic based  
mod.1 = mmer(Yield ~ 1,
             random= ~ vsr(id, Gu=G,Gtc=unsm(1)) + Rowf,
             rcov = ~ units,
             data = data.1.2,verbose = F)
# Model summary
(summary(mod.1)$varcomp)
# Narrow-sense heritability
h2b <- (mod.1$sigma$`u:id`+mod.1$sigma$Rowf) /(mod.1$sigma$`u:id`+mod.1$sigma$Rowf+mod.1$sigma$units); h2b #  broad-sense manually
h2a <- (mod.1$sigma$`u:id`) /(mod.1$sigma$`u:id`+mod.1$sigma$units); h2a # narrow sense manually
vpredict(mod.1, h2 ~ (V1+V2) / ( V1+V2+V3) ) # broad-sense function
vpredict(mod.1, h2 ~ (V1) / ( V1+V3) ) # narrow sense function
######################################################################################
# Fitting models using sommer package to compute narrow-sense heritability estimates for Fruit average
# No  marker information  
mod.1 = mmer(FruitAver ~ 1,
             random= ~ vsr(id,Gtc=unsm(1)) + Rowf,
             rcov = ~ units,
             data = data.1.2,verbose = F)
# Model summary
(summary(mod.1)$varcomp)
# Narrow-sense heritability
h2b <- (mod.1$sigma$`u:id`+mod.1$sigma$Rowf) /(mod.1$sigma$`u:id`+mod.1$sigma$Rowf+mod.1$sigma$units); h2b #  broad-sense manually
h2a <- (mod.1$sigma$`u:id`) /(mod.1$sigma$`u:id`+mod.1$sigma$units); h2a # narrow sense manually
vpredict(mod.1, h2 ~ (V1+V2) / ( V1+V2+V3) ) # broad-sense function
vpredict(mod.1, h2 ~ (V1) / ( V1+V3) ) # narrow sense function
# Genomic based  
mod.1 = mmer(FruitAver ~ 1,
             random= ~ vsr(id, Gu=G,Gtc=unsm(1)) + Rowf,
             rcov = ~ units,
             data = data.1.2,verbose = F)
# Model summary
(summary(mod.1)$varcomp)
# Narrow-sense heritability
h2b <- (mod.1$sigma$`u:id`+mod.1$sigma$Rowf) /(mod.1$sigma$`u:id`+mod.1$sigma$Rowf+mod.1$sigma$units); h2b #  broad-sense manually
h2a <- (mod.1$sigma$`u:id`) /(mod.1$sigma$`u:id`+mod.1$sigma$units); h2a # narrow sense manually
vpredict(mod.1, h2 ~ (V1+V2) / ( V1+V2+V3) ) # broad-sense function
vpredict(mod.1, h2 ~ (V1) / ( V1+V3) ) # narrow sense function
######################################################################################
# Cross-validation Single-trait models with no marker information using sommer for Yield
# Indentity matrix
I <- diag(nrow(data.1.2))
colnames(I) <- rownames(I) <- rownames(DT)
#30 Replicates
rep <- 3
#5 random partitions
fold <- 5
n <- dim(data.1.2)[1]
Replicates <- data.frame(Replicate = 1:rep,MSEP = NA,Cor = NA) 
for (j in 1:rep){ 
  set.seed(j)
  Partitions <- replicate(fold,sample(n,0.20*n,replace = F))
  Table <- data.frame(Partitions = 1:fold,MSEP = NA,Cor = NA)
  for(i in 1:fold) {
    tst <- Partitions[,i]
    data.1.2$y_NA <- data.1.2$Yield
    data.1.2$y_NA[tst] <- NA
    #Fit model
    fit = mmer(y_NA~ 1,
               random=~ vsr(id, Gu=I, Gtc=unsm(1)),
               rcov=~units,
               data=data.1.2,verbose = FALSE)
    #BLUPs
    bv <- fit$U$`u:id`$y_NA
    #Prediction of testing
    yp_tst <- bv[tst]
    #MSEP and Cor
    Table$MSEP[i] <- mean((data.1.2$Yield[tst]-yp_tst)^2)
    Table$Cor[i] <- cor(data.1.2$Yield[tst],yp_tst)
  }
  Replicates$MSEP[j] <- mean((na.omit(Table$MSEP))) # mean of folds
  Replicates$Cor[j] <-mean((na.omit(Table$Cor))) # mean of folds
}

mean(Replicates$MSEP) # mean of replicates
mean(Replicates$Cor) # mean of replicates
######################################################################################
######################################################################################
# Cross-validation Single-trait models (GBLUP) using sommer for Yield
#30 Replicates
rep <- 3
#5 random partitions
fold <- 5
n <- dim(data.1.2)[1]
Replicates <- data.frame(Replicate = 1:rep,MSEP = NA,Cor = NA) 
for (j in 1:rep){ 
  set.seed(j)
  Partitions <- replicate(fold,sample(n,0.20*n,replace = F))
  Table <- data.frame(Partitions = 1:fold,MSEP = NA,Cor = NA)
  for(i in 1:fold) {
    tst <- Partitions[,i]
    data.1.2$y_NA <- data.1.2$Yield
    data.1.2$y_NA[tst] <- NA
    #Fit model
    fit = mmer(y_NA ~ 1,
               random=~vsr(id,Gu=G,Gtc=unsm(1)) + vsr(Rowf),
               rcov=~units,
               data=data.1.2,verbose = FALSE)
    #BLUPs
    bv <- fit$U$`u:id`$y_NA
    #Prediction of testing
    yp_tst <- bv[tst]
    #MSEP and Cor
    Table$MSEP[i] <- mean((data.1.2$Yield[tst]-yp_tst)^2)
    Table$Cor[i] <- cor(data.1.2$Yield[tst],yp_tst)
  }
  Replicates$MSEP[j] <- mean((na.omit(Table$MSEP))) # mean of folds
  Replicates$Cor[j] <-mean((na.omit(Table$Cor))) # mean of folds
}

mean(Replicates$MSEP) # mean of replicates
mean(Replicates$Cor) # mean of replicates
######################################################################################
######################################################################################
# Cross-validation Single-trait models (GBLUP) using sommer for FruitAver
#30 Replicates
rep <- 3
#5 random partitions
fold <- 5
n <- dim(data.1.2)[1]
Replicates <- data.frame(Replicate = 1:rep,MSEP = NA,Cor = NA) 
for (j in 1:rep){ 
  set.seed(j)
  Partitions <- replicate(fold,sample(n,0.20*n,replace = F))
  Table <- data.frame(Partitions = 1:fold,MSEP = NA,Cor = NA)
  for(i in 1:fold) {
    tst <- Partitions[,i]
    data.1.2$y_NA <- data.1.2$FruitAver
    data.1.2$y_NA[tst] <- NA
    #Fit model
    fit = mmer(y_NA ~ 1,
               random=~vsr(id,Gu=G,Gtc=unsm(1)) + Rowf,
               rcov=~units,
               data=data.1.2,verbose = FALSE)
    #BLUPs
    bv <- fit$U$`u:id`$y_NA
    #Prediction of testing
    yp_tst <- bv[tst]
    #MSEP and Cor
    Table$MSEP[i] <- mean((data.1.2$Yield[tst]-yp_tst)^2)
    Table$Cor[i] <- cor(data.1.2$Yield[tst],yp_tst)
  }
  Replicates$MSEP[j] <- mean((na.omit(Table$MSEP))) # mean of folds
  Replicates$Cor[j] <-mean((na.omit(Table$Cor))) # mean of folds
}

mean(Replicates$MSEP) # mean of replicates
mean(Replicates$Cor) # mean of replicates
################################################################################
################################################################################
# Cross-validation Multi-trait models (GBLUP) using sommer for Yield and FruitAver
#30 Replicates
rep <- 3
#5 random partitions
fold <- 5
n <- dim(data.1.2)[1]
Replicates <- data.frame(Replicate = 1:rep,MSEP_T1 = NA,Cor_T1 = NA,
                         MSEP_T2 = NA,Cor_T2 = NA) 

for (j in 1:rep){ 
  set.seed(j)
  Partitions <- replicate(fold,sample(n,0.20*n,replace = F))
  Table <- data.frame(Partitions = 1:fold,MSEP_T1 = NA,Cor_T1 = NA,
                      MSEP_T2 = NA,Cor_T2 = NA)
  for(i in 1:fold) {
    tst <- Partitions[,i]
    data.1.2$y_NA.1 <- data.1.2$Yield
    data.1.2$y_NA.1[tst] <- NA
    data.1.2$y_NA.2 <- data.1.2$FruitAver
    data.1.2$y_NA.2[tst] <- NA
    #Fit model
    fit <- mmer(cbind(y_NA.1,y_NA.2)~ 1,
                random= ~ vsr(id, Gu=G,Gtc=unsm(2)),
                rcov= ~ vsr(units,Gtc=diag(2)), data=data.1.2,
                verbose=F, method = "AI")
    #BLUPs
    b_ls <- fit$U$`u:id`
    b_T1 <- b_ls$y_NA.1 # Trait 1
    b_T2 <- b_ls$y_NA.2 # Trait 2
    b_mat <- cbind(b_T1,b_T2)
    #Prediction of testing for both traits
    yp_tst = b_mat[tst,1] # Trait 1
    y2p_tst = b_mat[tst,2] # Trait 2
    #MSEP and Cor
    #Trait 1
    Table$MSEP_T1[i] <- mean((data.1.2$Yield[tst]-yp_tst)^2)
    Table$Cor_T1[i] <- cor(data.1.2$Yield[tst],yp_tst)
    #Trait 2
    Table$MSEP_T2[i] <- mean((data.1.2$FruitAver[tst]-y2p_tst)^2)
    Table$Cor_T2[i] <- cor(data.1.2$FruitAver[tst],y2p_tst)
  }
  Replicates$MSEP_T1[j] <- mean(na.omit(Table$MSEP_T1))
  Replicates$Cor_T1[j] <-mean(na.omit(Table$Cor_T1))
  Replicates$MSEP_T2[j] <- mean(na.omit(Table$MSEP_T2))
  Replicates$Cor_T2[j] <-mean(na.omit(Table$Cor_T2))
}

mean(Replicates$MSEP_T1)
mean(Replicates$Cor_T1)
mean(Replicates$MSEP_T2)
mean(Replicates$Cor_T2)
################################################################################