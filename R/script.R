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
# GBLUP model
A = mmer(y ~ 1,random= ~ vs(GID,Gu=G) , rcov= ~ vs(units),
         data=dat_F, verbose=FALSE)
################################################################################
# Sometimes is important to estimate genetic variance-covariance among traits–multi-reponse
# models. Let see an example with 3 traits (color, Yield, and Firmness) and a single random 
# effect (genotype; id) although multiple effects can be modeled as well. We need to use a
# variance covariance structure for the random effect to be able to obtain the genetic 
# covariance among traits.
data(DT_cpdata)
DT <- DT_cpdata
head(DT_cpdata)
GT <- GT_cpdata
GT_cpdata[1:5,1:5]
MP <- MP_cpdata
head(MP_cpdata)
A <- A.mat(GT)
A[1:5,1:5]
ans.m <- mmer(cbind(Yield,color)~1,
                random=~ vsr(id, Gu=A, Gtc=unsm(2))
                + vsr(Rowf,Gtc=diag(2))
                + vsr(Colf,Gtc=diag(2)),
                rcov=~ vsr(units, Gtc=unsm(2)), nIters=3,
                data=DT, verbose = FALSE)

# Now you can extract the BLUPs using randef(ans.m) or simply ans.m$U. Also, genetic 
# correlations and heritabilities can be calculated easily.
cov2cor(ans.m$sigma$`u:id`)

################################################################################
#fit an rrBLUP model in the multivariate setting, the DT_cpdata has a good example

librayr(sommer)
data(DT_cpdata)
mix.rrblup <- mmer(fixed=cbind(color,Yield)~1,
                   random=~vs(list(GT),Gtc=unsm(2)) + vs(Rowf,Gtc=diag(2))
                   rcov=~vs(units,Gtc=unsm(2)),
                   data=DT)
summary(mix.rrblup)

A <- A.mat(GT)
mix.gblup <- mmer(fixed=cbind(color,Yield)~1,
                  random=~vs(id,Gu=A, Gtc=unsm(2)) + vs(Rowf,Gtc=diag(2))
                  rcov=~vs(units,Gtc=unsm(2)),
                  data=DT)
summary(mix.gblup)

#the vs() function makes a variance structure for a given random effect, and the
#covariance structure for the univariate/multivariate setting is provided in the 
#Gtc argument as a matrix, where i.e. diagonal, unstructured or a customized 
#structure can be provided. When the user want to provide a customized matrix as
#a random effect such as a marker matrix GT to do rrBLUP it has to be provided 
#in a list() to internally help sommer to put it in the right format, whereas in
#the GBLUP version the random effect id which has the labels for individuals can
#have a covariance matrix provided in the Gu argument.

#################################################################################

#1) Univariate homogeneous variance models

#This type of model refers to single response models where a variable of interest (i.e. genotypes) needs to be
#analyzed as interacting with a 2nd random effect (i.e. environments), but you assume that across environments
#the genotypes have the same variance component. This is the so-called compound symmetry (CS) model.
data(DT_example)
DT <- DT_example
## solving for r records
ans1r <- mmer(Yield~Env,
              random= ~ Name + Env:Name,
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ans1r)$varcomp
## VarComp VarCompSE Zratio Constraint
## Name.Yield-Yield 3.681877 1.6909561 2.177394 Positive
## Env:Name.Yield-Yield 5.173062 1.4952313 3.459707 Positive
## units.Yield-Yield 4.366285 0.6470458 6.748031 Positive
## solving for c coefficients
ans1c <- mmec(Yield~Env,
              random= ~ Name + Env:Name,
              rcov= ~ units,
              data=DT, verbose = FALSE)
summary(ans1c)$varcomp
################################################################################

# Cross-validation Single-trait models (PBLUP) using rrBLUP
r2 <- array()
for (i in 1:5){
  fold <- 5
  partition <- sample(1:fold, size=nrow(Ag), replace = TRUE)
  r1 <- array()
  for (j in 1:fold){
    y.na <- as.matrix(cy) #cy
    tst <- which(partition == j) # index of the positions of individuals for each fold
    y.na[tst] <- NA
    fit <- mixed.solve(y=y.na, K=Ag)
    r1[j] <- cor(cy[tst], fit$u[tst])
  }
  r2[i] <- mean(r1)
}
mean(r2)

###############################################
# Fitting models using sommer package to compute narrow-sense heritability estimates
# Genomic based  
mod.5 <- mmer(BLUEs_TCHr ~ 1,
              random= ~ vs(Genotype_code, Gu=G),
              rcov = ~ vs(units),
              data = data.1.2,
              verbose = F)
# Narrow-sense heritability
h25 <- (as.double(mod.5$sigma$`u:Genotype_code`))/(as.double(mod.5$sigma$`u:Genotype_code`) + as.double(mod.5$sigma$`u:units`)); h25 
sigg5<- as.double(mod.5$sigma_scaled$`u:Genotype_code`)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
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
# Fitting models using sommer package to compute narrow-sense heritability estimates
# Genomic based  
mod.1 = mmer(Yield ~ 1,
              random= ~ vsr(id, Gu=G,Gtc=unsm(1)) + Rowf + Colf,
              rcov = ~ units,
              data = data.1.2,verbose = F)
# Narrow-sense heritability
h2 <- (as.double(mod.1$sigma$`u:id`))/(as.double(mod.1$sigma$`u:id`)) + (as.double(mod.1$sigma$units)); h2 
sigg<- as.double(mod.1$sigma_scaled$`u:id`)
######################################################################################
# Cross-validation Single-trait models (GBLUP) using sommer
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
         random=~vsr(id,Gu=G) + vsr(Rowf),
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
# Cross-validation Multi-trait models (GBLUP) using sommer: Scenario MT-1
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
    #Narrow sense heritability estimates
    #Trait 1 
    h21 <- (summary(fit)$varcomp[[1]][1])/(summary(fit)$varcomp[[1]][1] + summary(fit)$varcomp[[1]][4])
    #Trait 2 
    h22 <- (summary(fit)$varcomp[[1]][3])/(summary(fit)$varcomp[[1]][3] + summary(fit)$varcomp[[1]][5])
    #Genetic correlation
    rg <- (summary(fit)$varcomp[[1]][2])/sqrt((summary(fit)$varcomp[[1]][1] + summary(fit)$varcomp[[1]][3]))
    #cov2cor(fit$sigma$`u:Genotype_code`)
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
  Replicates$MSEP_T1[j] <- mean(Table$MSEP_T1)
  Replicates$Cor_T1[j] <-mean(Table$Cor_T1)
  Replicates$MSEP_T2[j] <- mean(Table$MSEP_T2)
  Replicates$Cor_T2[j] <-mean(Table$Cor_T2)
}

mean(Replicates$MSEP)
mean(Replicates$Cor)
################################################################################
#########################################################################
# DT_example
data(DT_example)
DT <- DT_example
DT$Yield <- as.vector(scale(DT$Yield))
DT$Weight <- as.vector(scale(DT$Weight))
# Single
ans3r <- mmer(Yield ~ 1,
              random= ~ vsr(Name),
              rcov= ~ vsr(units),
              data=DT, verbose = FALSE)
summary(ans3r)$varcomp
ans3r$sigma$`u:Name`
#Multi
ans4r <- mmer(cbind(Yield, Weight) ~ 1,
              random= ~ vsr(Name, Gtc=unsm(2)),
              rcov= ~ vsr(units, Gtc=diag(2)),
              data=DT, verbose = FALSE)
summary(ans4r)$varcomp
ans4r$sigma$`Name`
######################################################
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
#colnames(G) <- rownames(G) <- rownames(DT)
G[1:5,1:5]
# Creating new dataset
data.1.2 <- DT
# Scaling phenotypes
data.1.2$Yield <- as.vector(scale(data.1.2$Yield, center = T, scale = T)) # Yield
data.1.2$FruitAver <- as.vector(scale(data.1.2$FruitAver, center = T, scale = T)) # fruit avarage
data.1.2$Firmness <- as.vector(scale(data.1.2$Firmness, center = T, scale = T)) # firmness
data.1.2$color <- as.vector(scale(data.1.2$color, center = T, scale = T)) # color
##########################################################################################
# original
# DT_cpdata
data(DT_cpdata)
DT <- DT_cpdata
head(DT,3)
GT <- GT_cpdata
GT[1:5,1:5]
#### create the variance-covariance matrix
G <- A.mat(GT) # additive relationship matrix
DT$Yield <- as.vector(scale(DT$Yield))
DT$FruitAver <- as.vector(scale(DT$FruitAver))
data.1.2 <- DT # creating newdataset
# Single
ans3r <- mmer(Yield ~ 1,
              random= ~ vsr(id),
              rcov= ~ vsr(units),
              data=DT, verbose = FALSE)
summary(ans3r)$varcomp
ans3r$sigma$`u:id`
#Multi
ans4r <- mmer(cbind(Yield, FruitAver) ~ 1,
              random= ~ vsr(id, Gtc=unsm(2)),
              rcov= ~ vsr(units, Gtc=diag(2)),
              data=DT, verbose = FALSE)
summary(ans4r)$varcomp
ans4r$sigma$`u:id`
# Single + G
ans3r <- mmer(Yield ~ 1,
              random= ~ vsr(id, Gu=G),
              rcov= ~ vsr(units),
              data=DT, verbose = FALSE)
summary(ans3r)$varcomp
ans3r$sigma$`u:id`
#Multi + G
ans4r <- mmer(cbind(Yield, FruitAver) ~ 1,
              random= ~ vsr(id, Gu=G, Gtc=unsm(2)),
              rcov= ~ vsr(units, Gtc=diag(2)),
              data=DT, verbose = FALSE)
summary(ans4r)$varcomp
ans4r$sigma$`u:id`
#################################################################################
# Cross-validation Multi-trait models (GBLUP) using sommer: Scenario MT-1
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
