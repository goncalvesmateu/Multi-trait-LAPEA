################################################################################
#
# Script of multitrait analysis using sommer
#
# by Mateus Gonçalves
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
# Dataset
data(DT_cpdata)
# Phenotypic data
DT <- DT_cpdata
head(DT)
# Marker data
GT <- GT_cpdata
GT_cpdata[1:5,1:5]
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
######################################################################################
# No  marker information  
mod.1 = mmer(Yield ~ 1,
             random= ~ vsr(id,Gtc=unsm(1)) + Rowf, # 1:estimated and constrained to be positive
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
# Marker information    
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
######################################################################################
# No marker information  
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
######################################################################################
# Marker information    
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
# NOT WORKING NOT WORKING NOT WORKING NOT WORKING
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
               random=~ vsr(id, Gu=I, Gtc=unsm(2)), # unstructured model (each pair of random effects or residuals to have its own distinct covariance value)
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
######################################################################################
# Fitting models using sommer package to compute genetic correlation estimates for Yield
######################################################################################
#Fit model
fit <- mmer(cbind(Yield,FruitAver)~ 1,
            random= ~ vsr(id, Gu=G,Gtc=unsm(2)), # 2: estimated and unconstrained
            rcov= ~ vsr(units,Gtc=diag(2)), data=data.1.2, # 2: estimated and unconstrained
            verbose=F)

# Model summary
summary(fit)$varcomp
fit$sigma$`u:id`

################################################################################