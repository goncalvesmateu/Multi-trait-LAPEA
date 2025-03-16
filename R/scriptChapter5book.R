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
fit an rrBLUP model in the multivariate setting, the DT_cpdata has a good example

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