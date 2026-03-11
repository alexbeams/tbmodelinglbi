# See how well Local Branching Index classifies a variant with elevated transmission rate
# Can vary taulbis over a different range to make sure we are choosing the most generous
# bandwidth possible

 
require(ape)

source('lbi.R')




getmets <- function(tree){

	ntests <- Ntip(tree)

	# want to add in rows for nodes with times and LBIs
	crud <- data.frame(label = tree$tip.label)
	
	# Uncomment to use node information -- but it really only makes sense to maximize F statistics using tip data:
#	# the node labels have the times; extract these:
#	m<- sapply(tree$node.label, function(z) substr(z, regexpr('=',z)[1]+1, regexpr(',re',z)[1]-1 ) ) 
#
#	crud2 <- data.frame(label = names(m))
#
#	# the node labels are super clunky, but we need to keep them to match with the tree
#	crud <- rbind(crud,crud2)

	# calculate LBI for the tips and the nodes for several values of bandwidth, tau:
	taulbis <- seq(1,100,length=50)

	for(taulbi in taulbis){crud[,paste0('lbi',taulbi)] <- lbi(tree, tau=taulbi)[1:ntests]}

	# add in a column for the state of the node/tip:
	crud$state <- NA
	crud[grep('IH',crud$label[1:ntests]),'state'] <- 'I'
	crud[grep('IL',crud$label[1:ntests]),'state'] <- 'I'
	crud[grep('JH',crud$label[1:ntests]),'state'] <- 'J'
	crud[grep('JL',crud$label[1:ntests]),'state'] <- 'J'


#	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], 
#		function(z) substr(z, regexpr("S+",z)[1]+2, regexpr(".[+]=",z)[1]+0))
#
#	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms

	# Calculate F statistics from analysis of variance for each tau:
	Fstats <- apply(crud[,paste0('lbi',taulbis)], 2, function(x) summary(aov(x~state,crud))[[1]][['F value']][1])

	lbidat <- crud[1:ntests,c('label','state')]
	lbidat$lbi <- crud[1:ntests,1+which(Fstats==max(Fstats))] 

	getPredMetrics <- function(thresh){
		# Now, try classifying stuff
		preds <- lbidat$state
		preds[lbidat$lbi > quantile(lbidat$lbi,thresh)] <- 'J'
		preds[lbidat$lbi < quantile(lbidat$lbi,thresh)] <- 'I'

		preds <- as.data.frame(cbind(lbidat$state, preds))
		colnames(preds) <- c('true','pred')

		# calculate true positives, false positives, precision, and accuracy:
		tab <- table(preds)
		#tp <- tab[1,1]
		#fp <- tab[2,1]
		#tn <- tab[2,2]
		#fn <- tab[1,2]
	
		tn <- tab[1,1]
		fn <- tab[2,1]
		tp <- tab[2,2]
		fp <- tab[1,2]



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


# Commenting out the following portion, which calculates
# false positives, true positives, etc. from the trees 
# but which is now saved - so we load the rocdat's below
# to construct the plots

# What are the paths for all of the trees we want to load?



####### 10%
#######

treedir <- 'sims/varmodel/10/'
treenms <- paste0('tree',1:30,'.nwk')
treenms <- paste0(treedir,treenms)

# make a big list where we will store all of the matrices to create our ROC plot
# One ROC curve for each configuration of tests (Nyr and ntests)
# ROC curves are averaged over 30 simulations
metlist <- list()

for(treenm in treenms){
	tree <- read.tree(treenm)
	met <- getmets(tree)
	metlist[[treenm]] <- met	
}

# average over the simulations:
rocdat_10 <- apply(simplify2array(metlist),1:2,mean)
save(rocdat_10,file='sims/varmodel/roc/rocdat_10.Rdata')

####### 20%
#######

treedir <- 'sims/varmodel/20/'
treenms <- paste0('tree',1:30,'.nwk')
treenms <- paste0(treedir,treenms)

# make a big list where we will store all of the matrices to create our ROC plot
# One ROC curve for each configuration of tests (Nyr and ntests)
# ROC curves are averaged over 30 simulations
metlist <- list()

for(treenm in treenms){
	tree <- read.tree(treenm)
	met <- getmets(tree)
	metlist[[treenm]] <- met	
}

# average over the simulations:
rocdat_20 <- apply(simplify2array(metlist),1:2,mean)
save(rocdat_20,file='sims/varmodel/roc/rocdat_20.Rdata')

####### 30%
#######

treedir <- 'sims/varmodel/30/'
treenms <- paste0('tree',1:30,'.nwk')
treenms <- paste0(treedir,treenms)

# make a big list where we will store all of the matrices to create our ROC plot
# One ROC curve for each configuration of tests (Nyr and ntests)
# ROC curves are averaged over 30 simulations
metlist <- list()

for(treenm in treenms){
	tree <- read.tree(treenm)
	met <- getmets(tree)
	metlist[[treenm]] <- met	
}

# average over the simulations:
rocdat_30 <- apply(simplify2array(metlist),1:2,mean)
save(rocdat_30,file='sims/varmodel/roc/rocdat_30.Rdata')

####### 40%
#######

treedir <- 'sims/varmodel/40/'
treenms <- paste0('tree',1:30,'.nwk')
treenms <- paste0(treedir,treenms)

# make a big list where we will store all of the matrices to create our ROC plot
# One ROC curve for each configuration of tests (Nyr and ntests)
# ROC curves are averaged over 30 simulations
metlist <- list()

for(treenm in treenms){
	tree <- read.tree(treenm)
	met <- getmets(tree)
	metlist[[treenm]] <- met	
}

# average over the simulations:
rocdat_40 <- apply(simplify2array(metlist),1:2,mean)
save(rocdat_40,file='sims/varmodel/roc/rocdat_40.Rdata')

####### 50%
#######

treedir <- 'sims/varmodel/50/'
treenms <- paste0('tree',1:30,'.nwk')
treenms <- paste0(treedir,treenms)

# make a big list where we will store all of the matrices to create our ROC plot
# One ROC curve for each configuration of tests (Nyr and ntests)
# ROC curves are averaged over 30 simulations
metlist <- list()

for(treenm in treenms){
	tree <- read.tree(treenm)
	met <- getmets(tree)
	metlist[[treenm]] <- met	
}

# average over the simulations:
rocdat_50 <- apply(simplify2array(metlist),1:2,mean)
save(rocdat_50,file='sims/varmodel/roc/rocdat_50.Rdata')




### We load the rocdat dataframes to create our plots:
#
load('sims/varmodel/roc/rocdat_10.Rdata')
load('sims/varmodel/roc/rocdat_20.Rdata')
load('sims/varmodel/roc/rocdat_30.Rdata')
load('sims/varmodel/roc/rocdat_50.Rdata')
load('sims/varmodel/roc/rocdat_50.Rdata')



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




# Save plot:
pdf(file='figures/roc/varmodel_roc.pdf',height=7,width=7)
par(mfrow=c(1,1))
plot(truepos~falsepos,rocdat_10,type='l',lwd=2,
	xlab='False-positive rate',
	ylab='True-positive rate',
	main='Variant classification with LBI (n=500 samples)',
	cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(truepos~falsepos,rocdat_20,lwd=2,lty=2)
lines(truepos~falsepos,rocdat_30,lwd=2,lty=3)
lines(truepos~falsepos,rocdat_40,lwd=2,lty=4)
lines(truepos~falsepos,rocdat_50,lwd=2,lty=5)
legend('bottomright',
	legend=c('10%','20%','30%','40%','50%'),
	lty=c(1,2,3,4,5),lwd=2,bty='n',cex=1.5,
	title='Variant frequency')
abline(a=0,b=1,lwd=1)
# Turn of plotting:
dev.off()


