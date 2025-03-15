################################################################################
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