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