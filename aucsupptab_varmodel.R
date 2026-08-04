# See how well Local Branching Index classifies a variant with elevated transmission rate
# Can vary taulbis over a different range to make sure we are choosing the most generous
# bandwidth possible

 
require(ape)
require(ggtree)
require(ggplot2)
require(lubridate)
require(ggnewscale)
require(cowplot)
require(xtable)

getauc <- function(seladv){
	# Which selection coefficient to use?
	#seladv <- 10

	### We load the rocdat dataframes to create our plots:
	#
	load(paste0('sims/varmodel/',seladv,'/roc/rocdat_10.Rdata'))
	load(paste0('sims/varmodel/',seladv,'/roc/rocdat_20.Rdata'))
	load(paste0('sims/varmodel/',seladv,'/roc/rocdat_30.Rdata'))
	load(paste0('sims/varmodel/',seladv,'/roc/rocdat_40.Rdata'))
	load(paste0('sims/varmodel/',seladv,'/roc/rocdat_50.Rdata'))


	# We need to calculate AUC from each rocdat by quadrature:

	getauc <- function(x,y){
		# given gridpoints at x, y, calculate area under the curve y=y(x)
		# by averaging the left and right-endpoint quadratures:
		left <- sum(diff(c(0,x,1))*c(0,y))
		right <- sum(diff(c(0,x,1))*c(y,1))
		auc <- mean(c(left,right))
		return(auc)
	}

	auc_10 <- with(as.data.frame(rocdat_10), getauc(rev(falsepos),rev(truepos)) )
	auc_20 <- with(as.data.frame(rocdat_20), getauc(rev(falsepos),rev(truepos)) )
	auc_30 <- with(as.data.frame(rocdat_30), getauc(rev(falsepos),rev(truepos)) )
	auc_40 <- with(as.data.frame(rocdat_40), getauc(rev(falsepos),rev(truepos)) )
	auc_50 <- with(as.data.frame(rocdat_50), getauc(rev(falsepos),rev(truepos)) )

	auc <-  c(auc_10, auc_20, auc_30, auc_40, auc_50)

	return(auc)
}

auc.10 <- getauc(10)
auc.15 <- getauc(15)
auc.20 <- getauc(20)
auc.25 <- getauc(25)

aucdat <- cbind(auc.10, auc.15, auc.20, auc.25)
aucdat <- as.data.frame(aucdat)
colnames(aucdat) <- c(0.10,0.15,0.20,0.25)
rownames(aucdat) <- paste0(c(10,20,30,40,50), '%')

xtable(aucdat)







