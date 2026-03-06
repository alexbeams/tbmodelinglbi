# See how well Local Branching Index classifies superspreaders
require(ape)

source('lbi.R')

# read in a simulated tree
#treepath <- 'sims/nullmodel/500/5yr/tree3.nwk'
#tree <- read.tree(treepath)

getmets <- function(tree){

	ntests <- Ntip(tree)

	# want to add in rows for nodes with times and LBIs
	crud <- data.frame(label = tree$tip.label)

	# the node labels have the times; extract these:
	m<- sapply(tree$node.label, function(z) substr(z, regexpr('=',z)[1]+1, regexpr(',re',z)[1]-1 ) ) 

	crud2 <- data.frame(label = names(m))

	# the node labels are super clunky, but we need to keep them to match with the tree
	crud <- rbind(crud,crud2)

	# calculate LBI for the tips and the nodes for several values of bandwidth, tau:
	taulbis <- 1:10

	for(taulbi in taulbis){crud[,paste0('lbi',taulbi)] <- lbi(tree, tau=taulbi)}

	# add in a column for the state of the node/tip:
	crud$state <- NA
	crud[grep('IH',crud$label[1:ntests]),'state'] <- 'IH'
	crud[grep('IL',crud$label[1:ntests]),'state'] <- 'IL'

	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], 
		function(z) substr(z, regexpr("S+",z)[1]+2, regexpr(".[+]=",z)[1]+0))

	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms

	# Calculate F statistics from analysis of variance for each tau:
	Fstats <- apply(crud[,paste0('lbi',taulbis)], 2, function(x) summary(aov(x~state,crud))[[1]][['F value']][1])

	lbidat <- crud[1:ntests,c('label','state')]
	lbidat$lbi <- crud[1:ntests,1+which(Fstats==max(Fstats))] 

	getPredMetrics <- function(thresh){
		# Now, try classifying stuff
		preds <- lbidat$state
		preds[lbidat$lbi > quantile(lbidat$lbi,thresh)] <- 'IH'
		preds[lbidat$lbi < quantile(lbidat$lbi,thresh)] <- 'IL'

		preds <- as.data.frame(cbind(lbidat$state, preds))
		colnames(preds) <- c('true','pred')

		# calculate true positives, false positives, precision, and accuracy:
		tab <- table(preds)
		tp <- tab[1,1]
		fp <- tab[2,1]
		tn <- tab[2,2]
		fn <- tab[1,2]


		prc <- tp/(tp+fp)
		acc <- (tp+tn)/(tp+tn+fp+fn)

		# turn false positives and true positives into probabilities:
		fpr <- fp/(fp+tn)
		tpr <- tp/(tp+fn)

		return(c('falsepos'=fpr,'truepos'=tpr,'precision'=prc,'accuracy'=acc))
	}

	thresholds <- seq(0.001,0.999,length=20)

	mets <- t(sapply(thresholds, getPredMetrics))
	mets <- cbind(thresholds, mets)
	colnames(mets)[1] <- 'thresh'

	return(mets)
}


## Commenting out the following portion, which calculates
## false positives, true positives, etc. from the trees 
## but which is now saved - so we load the rocdat's below
## to construct the plots

## What are the paths for all of the trees we want to load?
#
#
#
######## 5 years
########
#
## ntests=500, Nyr = 5
#treedir <- 'sims/nullmodel/500/5yr/'
#treenms <- paste0('tree',1:30,'.nwk')
#treenms <- paste0(treedir,treenms)
#
## make a big list where we will store all of the matrices to create our ROC plot
## One ROC curve for each configuration of tests (Nyr and ntests)
## ROC curves are averaged over 30 simulations
#metlist <- list()
#
#for(treenm in treenms){
#	tree <- read.tree(treenm)
#	met <- getmets(tree)
#	metlist[[treenm]] <- met	
#}
#
## average over the simulations:
#rocdat_500_5yr <- apply(simplify2array(metlist),1:2,mean)
#save(rocdat_500_5yr,file='sims/nullmodel/rocdat_500_5yr.Rdata')
#
## ntests=2000, Nyr = 5
#treedir <- 'sims/nullmodel/2000/5yr/'
#treenms <- paste0('tree',1:30,'.nwk')
#treenms <- paste0(treedir,treenms)
#
## make a big list where we will store all of the matrices to create our ROC plot
## One ROC curve for each configuration of tests (Nyr and ntests)
## ROC curves are averaged over 30 simulations
#metlist <- list()
#
#for(treenm in treenms){
#	tree <- read.tree(treenm)
#	met <- getmets(tree)
#	metlist[[treenm]] <- met	
#}
#
## average over the simulations:
#rocdat_2000_5yr <- apply(simplify2array(metlist),1:2,mean)
#save(rocdat_2000_5yr,file='sims/nullmodel/rocdat_2000_5yr.Rdata')
#
#
######## 10 years
########
#
## ntests=500, Nyr = 10
#treedir <- 'sims/nullmodel/500/10yr/'
#treenms <- paste0('tree',1:30,'.nwk')
#treenms <- paste0(treedir,treenms)
#
## make a big list where we will store all of the matrices to create our ROC plot
## One ROC curve for each configuration of tests (Nyr and ntests)
## ROC curves are averaged over 30 simulations
#metlist <- list()
#
#for(treenm in treenms){
#	tree <- read.tree(treenm)
#	met <- getmets(tree)
#	metlist[[treenm]] <- met	
#}
#
## average over the simulations:
#rocdat_500_10yr <- apply(simplify2array(metlist),1:2,mean)
#save(rocdat_500_10yr,file='sims/nullmodel/rocdat_500_10yr.Rdata')
#
#
#
## ntests=2000, Nyr = 10
#treedir <- 'sims/nullmodel/2000/10yr/'
#treenms <- paste0('tree',1:30,'.nwk')
#treenms <- paste0(treedir,treenms)
#
## make a big list where we will store all of the matrices to create our ROC plot
## One ROC curve for each configuration of tests (Nyr and ntests)
## ROC curves are averaged over 30 simulations
#metlist <- list()
#
#for(treenm in treenms){
#	tree <- read.tree(treenm)
#	met <- getmets(tree)
#	metlist[[treenm]] <- met	
#}
#
## average over the simulations:
#rocdat_2000_10yr <- apply(simplify2array(metlist),1:2,mean)
#save(rocdat_2000_10yr,file='sims/nullmodel/rocdat_2000_10yr.Rdata')
#
#
#
######## 15 years
########
#
## ntests=500, Nyr = 15
#treedir <- 'sims/nullmodel/500/15yr/'
#treenms <- paste0('tree',1:30,'.nwk')
#treenms <- paste0(treedir,treenms)
#
## make a big list where we will store all of the matrices to create our ROC plot
## One ROC curve for each configuration of tests (Nyr and ntests)
## ROC curves are averaged over 30 simulations
#metlist <- list()
#
#for(treenm in treenms){
#	tree <- read.tree(treenm)
#	met <- getmets(tree)
#	metlist[[treenm]] <- met	
#}
#
## average over the simulations:
#rocdat_500_15yr <- apply(simplify2array(metlist),1:2,mean)
#save(rocdat_500_15yr,file='sims/nullmodel/rocdat_500_15yr.Rdata')
#
#
#
## ntests=2000, Nyr = 15
#treedir <- 'sims/nullmodel/2000/15yr/'
#treenms <- paste0('tree',1:30,'.nwk')
#treenms <- paste0(treedir,treenms)
#
## make a big list where we will store all of the matrices to create our ROC plot
## One ROC curve for each configuration of tests (Nyr and ntests)
## ROC curves are averaged over 30 simulations
#metlist <- list()
#
#for(treenm in treenms){
#	tree <- read.tree(treenm)
#	met <- getmets(tree)
#	metlist[[treenm]] <- met	
#}
#
## average over the simulations:
#rocdat_2000_15yr <- apply(simplify2array(metlist),1:2,mean)
#save(rocdat_2000_15yr,file='sims/nullmodel/rocdat_2000_15yr.Rdata')
#
#
#
######## 20 years
########
#
## ntests=500, Nyr = 20
#treedir <- 'sims/nullmodel/500/20yr/'
#treenms <- paste0('tree',1:30,'.nwk')
#treenms <- paste0(treedir,treenms)
#
## make a big list where we will store all of the matrices to create our ROC plot
## One ROC curve for each configuration of tests (Nyr and ntests)
## ROC curves are averaged over 30 simulations
#metlist <- list()
#
#for(treenm in treenms){
#	tree <- read.tree(treenm)
#	met <- getmets(tree)
#	metlist[[treenm]] <- met	
#}
#
## average over the simulations:
#rocdat_500_20yr <- apply(simplify2array(metlist),1:2,mean)
#save(rocdat_500_20yr,file='sims/nullmodel/rocdat_500_20yr.Rdata')
#
#
#
## ntests=2000, Nyr = 20
#treedir <- 'sims/nullmodel/2000/20yr/'
#treenms <- paste0('tree',1:30,'.nwk')
#treenms <- paste0(treedir,treenms)
#
## make a big list where we will store all of the matrices to create our ROC plot
## One ROC curve for each configuration of tests (Nyr and ntests)
## ROC curves are averaged over 30 simulations
#metlist <- list()
#
#for(treenm in treenms){
#	tree <- read.tree(treenm)
#	met <- getmets(tree)
#	metlist[[treenm]] <- met	
#}
#
## average over the simulations:
#rocdat_2000_20yr <- apply(simplify2array(metlist),1:2,mean)
#save(rocdat_2000_20yr,file='sims/nullmodel/rocdat_2000_20yr.Rdata')


## We load the rocdat dataframes to create our plots:

load('sims/nullmodel/rocdat_500_5yr.Rdata')
load('sims/nullmodel/rocdat_500_10yr.Rdata')
load('sims/nullmodel/rocdat_500_15yr.Rdata')
load('sims/nullmodel/rocdat_500_20yr.Rdata')

load('sims/nullmodel/rocdat_2000_5yr.Rdata')
load('sims/nullmodel/rocdat_2000_10yr.Rdata')
load('sims/nullmodel/rocdat_2000_15yr.Rdata')
load('sims/nullmodel/rocdat_2000_20yr.Rdata')


# We need to calculate AUC from each rocdat by quadrature:

getauc <- function(x,y){
	# given gridpoints at x, y, calculate area under the curve y=y(x)
	# by averaging the left and right-endpoint quadratures:
	left <- sum(diff(c(0,x,1))*c(0,y))
	right <- sum(diff(c(0,x,1))*c(y,1))
	auc <- mean(c(left,right))
	return(auc)
}

auc_500_5yr <- with(as.data.frame(rocdat_500_5yr), getauc(rev(falsepos),rev(truepos)) )
auc_500_10yr <- with(as.data.frame(rocdat_500_10yr), getauc(rev(falsepos),rev(truepos)) )
auc_500_15yr <- with(as.data.frame(rocdat_500_15yr), getauc(rev(falsepos),rev(truepos)) )
auc_500_20yr <- with(as.data.frame(rocdat_500_20yr), getauc(rev(falsepos),rev(truepos)) )

auc_2000_5yr <- with(as.data.frame(rocdat_2000_5yr), getauc(rev(falsepos),rev(truepos)) )
auc_2000_10yr <- with(as.data.frame(rocdat_2000_10yr), getauc(rev(falsepos),rev(truepos)) )
auc_2000_15yr <- with(as.data.frame(rocdat_2000_15yr), getauc(rev(falsepos),rev(truepos)) )
auc_2000_20yr <- with(as.data.frame(rocdat_2000_20yr), getauc(rev(falsepos),rev(truepos)) )




pdf(file='figures/nullmodel_roc_2.pdf',height=6,width=12)
par(mfrow=c(1,2))
plot(truepos~falsepos,rocdat_500_5yr,type='l',lwd=2,
	xlab='False-positive rate',
	ylab='True-positive rate',
	main='500 infections sampled',
	cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(truepos~falsepos,rocdat_500_10yr,lwd=2,lty=2)
lines(truepos~falsepos,rocdat_500_15yr,lwd=2,lty=3)
lines(truepos~falsepos,rocdat_500_20yr,lwd=3,lty=4)
legend('bottomright',legend=c('5 yrs','10 yrs', '15 yrs', '20 yrs'),
	lty=c(1,2,3,4),lwd=2,bty='n',cex=2)
abline(a=0,b=1,lwd=1)

plot(truepos~falsepos,rocdat_2000_5yr,type='l',lwd=2,
	xlab='False-positive rate',
	ylab='True-positive rate',
	main='2000 infections sampled',
	cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(truepos~falsepos,rocdat_2000_10yr,lwd=2,lty=2)
lines(truepos~falsepos,rocdat_2000_15yr,lwd=2,lty=3)
lines(truepos~falsepos,rocdat_2000_20yr,lwd=3,lty=4)
abline(a=0,b=1,lwd=1)

dev.off()

# try plotting by duration:

pdf(file='figures/nullmodel_roc.pdf',height=10,width=10)
par(mfrow=c(2,2))

plot(truepos~falsepos,rocdat_500_5yr,type='l',lwd=2,
	xlab='False-positive rate',
	ylab='True-positive rate',
	main='Sampling over 5 years',
	cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(truepos~falsepos,rocdat_2000_5yr,lwd=2,lty=2)
abline(a=0,b=1,lwd=1,lty='dotted')
#legend('bottomright',
#	legend=c('500 tests', '2000 tests'),
#	lty=c(1,2),lwd=2,bty='n',cex=1.5)
legend('bottomright',
	legend=c(
		bquote('500 tests, AUC' ==.(round(auc_500_5yr,2))), 
		bquote('2000 tests, AUC' ==.(round(auc_2000_5yr,2)))
		),
	lty=c(1,2),lwd=2,bty='n',cex=1.5)


plot(truepos~falsepos,rocdat_500_10yr,type='l',lwd=2,
	xlab='False-positive rate',
	ylab='True-positive rate',
	main='Sampling over 10 years',
	cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(truepos~falsepos,rocdat_2000_10yr,lwd=2,lty=2)
abline(a=0,b=1,lwd=1,lty='dotted')
legend('bottomright',
	legend=c(
		bquote('500 tests, AUC' ==.(round(auc_500_10yr,2))), 
		bquote('2000 tests, AUC' ==.(round(auc_2000_10yr,2)))
		),
	lty=c(1,2),lwd=2,bty='n',cex=1.5)


plot(truepos~falsepos,rocdat_500_15yr,type='l',lwd=2,
	xlab='False-positive rate',
	ylab='True-positive rate',
	main='Sampling over 15 years',
	cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(truepos~falsepos,rocdat_2000_15yr,lwd=2,lty=2)
abline(a=0,b=1,lwd=1,lty='dotted')
legend('bottomright',
	legend=c(
		bquote('500 tests, AUC' ==.(round(auc_500_15yr,2))), 
		bquote('2000 tests, AUC' ==.(round(auc_2000_15yr,2)))
		),
	lty=c(1,2),lwd=2,bty='n',cex=1.5)



plot(truepos~falsepos,rocdat_500_20yr,type='l',lwd=2,
	xlab='False-positive rate',
	ylab='True-positive rate',
	main='Sampling over 20 years',
	cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(truepos~falsepos,rocdat_2000_20yr,lwd=2,lty=2)
abline(a=0,b=1,lwd=1,lty='dotted')
legend('bottomright',
	legend=c(
		bquote('500 tests, AUC' ==.(round(auc_500_20yr,2))), 
		bquote('2000 tests, AUC' ==.(round(auc_2000_20yr,2)))
		),
	lty=c(1,2),lwd=2,bty='n',cex=1.5)


#turn off plotting functionality in the pdf call:
dev.off()

