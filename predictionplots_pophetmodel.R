# See how well Local Branching Index classifies infections according to host
# subpopulation (one group has elevated susceptibility to infection)


# Let's try maximizing F statistics of lbi~state at tips (but not nodes), ditto but with nodes,
#	and also just try maximizing variance in LBI at the tips
 
require(ape)

source('lbi.R')

# Calculate prediction metrics (true positives, false positives, etc.) at the lbi bandwidth
# that maximizes the F-statistic in aov(lbi~state,lbidat) at various percentiles of lbi 
# to classify as "high LBI" or "low LBI"


getmets <- function(tree){

	ntests <- Ntip(tree)

	# want to add in rows for nodes with times and LBIs
	crud <- data.frame(label = tree$tip.label)

	# We don't need to select tau based on node values - just tips (but uncomment to change):
	## the node labels have the times; extract these:
	#m<- sapply(tree$node.label, function(z) substr(z, regexpr('=',z)[1]+1, regexpr(',re',z)[1]-1 ) ) 

	#crud2 <- data.frame(label = names(m))

	## the node labels are super clunky, but we need to keep them to match with the tree
	#crud <- rbind(crud,crud2)

	# calculate LBI for the tips and the nodes for several values of bandwidth, tau:
	taulbis <- seq(1,100,length=100)

	# just use LBI at tips:
	for(taulbi in taulbis){crud[,paste0('lbi',taulbi)] <- lbi(tree, tau=taulbi)[1:ntests]}

	# add in a column for the state of the node/tip:
	crud$state <- NA
	crud[grep('IH1',crud$label[1:ntests]),'state'] <- 'Group 1'
	crud[grep('IL1',crud$label[1:ntests]),'state'] <- 'Group 1'
	crud[grep('IH2',crud$label[1:ntests]),'state'] <- 'Group 2'
	crud[grep('IL2',crud$label[1:ntests]),'state'] <- 'Group 2'


#	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], 
#		function(z) substr(z, regexpr("S+",z)[1]+3, regexpr(".[+]=",z)[1]+0))
#
#	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms

	# Calculate F statistics from analysis of variance for each tau:
	Fstats <- apply(crud[,paste0('lbi',taulbis)], 2, function(x) summary(aov(x~state,crud))[[1]][['F value']][1])

	lbidat <- crud[1:ntests,c('label','state')]
	lbidat$lbi <- crud[1:ntests,1+which(Fstats==max(Fstats))] 

	getPredMetrics <- function(thresh){
		# Now, try classifying stuff
		preds <- lbidat$state
		preds[lbidat$lbi > quantile(lbidat$lbi,thresh)] <- 'Group 2'
		preds[lbidat$lbi < quantile(lbidat$lbi,thresh)] <- 'Group 1'

		preds <- as.data.frame(cbind(lbidat$state, preds))
		colnames(preds) <- c('true','pred')

		# calculate true positives, false positives, precision, and accuracy:
		tab <- table(preds)
		#tp <- tab[1,1]
		#fp <- tab[2,1]
		#tn <- tab[2,2]
		#fn <- tab[1,2]

		# We switched around the order - make it consistent:
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

	# Can change these thresholds to see if anything changes:
	thresholds <- seq(0.01,0.999,length=20)

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



####### Case 1, theta1: groups in population equal, homogeneous mixing:
#######

treedir <- 'sims/pophetmodel/case1/theta1/'
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
rocdat_case1_theta1 <- apply(simplify2array(metlist),1:2,mean)
save(rocdat_case1_theta1,file='sims/pophetmodel/roc/rocdat_case1_theta1.Rdata')

####### Case 1, theta2: groups in population equal, preferential mixing:
#######

treedir <- 'sims/pophetmodel/case1/theta2/'
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
rocdat_case1_theta2 <- apply(simplify2array(metlist),1:2,mean)
save(rocdat_case1_theta2,file='sims/pophetmodel/roc/rocdat_case1_theta2.Rdata')



####### Case 2, theta1: groups in population equal, homogeneous mixing:
#######

treedir <- 'sims/pophetmodel/case2/theta1/'
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
rocdat_case2_theta1 <- apply(simplify2array(metlist),1:2,mean)
save(rocdat_case2_theta1,file='sims/pophetmodel/roc/rocdat_case2_theta1.Rdata')

####### Case 2, theta2: groups in population equal, preferential mixing:
#######

treedir <- 'sims/pophetmodel/case2/theta2/'
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
rocdat_case2_theta2 <- apply(simplify2array(metlist),1:2,mean)
save(rocdat_case2_theta2,file='sims/pophetmodel/roc/rocdat_case2_theta2.Rdata')



### We load the rocdat dataframes to create our plots:
#
load('sims/pophetmodel/roc/rocdat_case1_theta1.Rdata')
load('sims/pophetmodel/roc/rocdat_case1_theta2.Rdata')
load('sims/pophetmodel/roc/rocdat_case2_theta1.Rdata')
load('sims/pophetmodel/roc/rocdat_case2_theta2.Rdata')


# We need to calculate AUC from each rocdat by quadrature:

getauc <- function(x,y){
	# given gridpoints at x, y, calculate area under the curve y=y(x)
	# by averaging the left and right-endpoint quadratures:
	left <- sum(diff(c(0,x,1))*c(0,y))
	right <- sum(diff(c(0,x,1))*c(y,1))
	auc <- mean(c(left,right))
	return(auc)
}

auc_case1_theta1 <- with(as.data.frame(rocdat_case1_theta1), getauc(rev(falsepos),rev(truepos)) )
auc_case1_theta2 <- with(as.data.frame(rocdat_case1_theta2), getauc(rev(falsepos),rev(truepos)) )
auc_case2_theta1 <- with(as.data.frame(rocdat_case2_theta1), getauc(rev(falsepos),rev(truepos)) )
auc_case2_theta2 <- with(as.data.frame(rocdat_case2_theta2), getauc(rev(falsepos),rev(truepos)) )






# Save plot:
pdf(file='figures/roc/pophetmodel_roc.pdf',height=7,width=14)
par(mfrow=c(1,2))
plot(truepos~falsepos,rocdat_case1_theta1,type='l',lwd=2,
	xlab='False-positive rate',
	ylab='True-positive rate',
	main='Using LBI to identify risk groups: Case 1',
	cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(truepos~falsepos,rocdat_case1_theta2,lwd=2,lty=2)
legend('bottomright',
	legend=c(bquote('Homogeneous mixing, AUC'==.(round(auc_case1_theta1,2))),
		bquote('Preferential mixing, AUC'==.(round(auc_case1_theta2,2)))
		),
	lty=c(1,2),lwd=2,bty='n',cex=1.5,
	title='Contact structure')
abline(a=0,b=1,lwd=1)

plot(truepos~falsepos,rocdat_case2_theta1,type='l',lwd=2,
	xlab='False-positive rate',
	ylab='True-positive rate',
	main='Using LBI to identify risk groups: Case 2',
	cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(truepos~falsepos,rocdat_case2_theta2,lwd=2,lty=2)
legend('bottomright',
	legend=c(
		bquote('Homogeneous mixing, AUC'==.(round(auc_case2_theta1,2))),
		bquote('Preferential mixing, AUC'==.(round(auc_case2_theta2,2)))
		),
	lty=c(1,2),lwd=2,bty='n',cex=1.5,
	title='Contact structure')
abline(a=0,b=1,lwd=1)

# Turn of plotting:
dev.off()
	

# Somewhat surprising that case 2 theta 2 does so poorly - why is that?	

# use treenms for case2, theta2:
treelist <- list()
for(treenm in treenms){treelist[[treenm]] <- read.tree(treenm)}


gettauopt <- function(tree){

	ntests <- Ntip(tree)

	# want to add in rows for nodes with times and LBIs
	crud <- data.frame(label = tree$tip.label)

	# We don't need to select tau based on node values - just tips (but uncomment to change):
	## the node labels have the times; extract these:
	#m<- sapply(tree$node.label, function(z) substr(z, regexpr('=',z)[1]+1, regexpr(',re',z)[1]-1 ) ) 

	#crud2 <- data.frame(label = names(m))

	## the node labels are super clunky, but we need to keep them to match with the tree
	#crud <- rbind(crud,crud2)

	# calculate LBI for the tips and the nodes for several values of bandwidth, tau:
	taulbis <- seq(1,100,length=100)

	# just use LBI at tips:
	for(taulbi in taulbis){crud[,paste0('lbi',taulbi)] <- lbi(tree, tau=taulbi)[1:ntests]}

	# add in a column for the state of the node/tip:
	crud$state <- NA
	crud[grep('IH1',crud$label[1:ntests]),'state'] <- 'Group 1'
	crud[grep('IL1',crud$label[1:ntests]),'state'] <- 'Group 1'
	crud[grep('IH2',crud$label[1:ntests]),'state'] <- 'Group 2'
	crud[grep('IL2',crud$label[1:ntests]),'state'] <- 'Group 2'


#	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], 
#		function(z) substr(z, regexpr("S+",z)[1]+3, regexpr(".[+]=",z)[1]+0))
#
#	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms

	# Calculate F statistics from analysis of variance for each tau:
	Fstats <- apply(crud[,paste0('lbi',taulbis)], 2, function(x) summary(aov(x~state,crud))[[1]][['F value']][1])
	tauopt <- taulbis[Fstats==max(Fstats)]
	return(tauopt)
}

tauopts <- sapply(treelist, gettauopt)

# set k
# getmainplot(treelist[[k]],tauopts[[k]])

# k=3 is a good example of Group 1 looking higher than Group 2 ( taulbi = 79 )
# k=29 is a good example of when Group 2 has higher LBi than Group 1  (taulbi = 36)
source('getpophetplotter.R')

p1 <- getmainplot(treelist[[3]], taulbi=79, title='Example: negative correlation\nbeween LBI and host susceptibility')[[1]] 
p2 <- getmainplot(treelist[[29]], taulbi=36, title='Example: positive correlation\nbetween LBI and host susceptibility')[[1]]

require(cowplot)

plot_example <- plot_grid(p1,p2,nrow=1)
ggsave(plot_example, file='figures/pophetmodelfigs/case2example.png',height=10,width=16)




###
# Choose some representative trees from the other simulations to make Figure 6:
###

# Inefficient, but calculate the optimal taus for each tree...

# Pick a good example from case 1, theta1:
treedir <- 'sims/pophetmodel/case1/theta1/'
treenms <- paste0('tree',1:30,'.nwk')
treenms <- paste0(treedir,treenms)

treelist <- list()

for(treenm in treenms){
	treelist[[treenm]] <- read.tree(treenm)
}

treelist_case1_theta1 <- treelist

tauopts_case1_theta1 <- sapply(treelist, gettauopt)

# Pick a good example from case 1, theta2:
treedir <- 'sims/pophetmodel/case1/theta2/'
treenms <- paste0('tree',1:30,'.nwk')
treenms <- paste0(treedir,treenms)

treelist <- list()

for(treenm in treenms){
	treelist[[treenm]] <- read.tree(treenm)
}

treelist_case1_theta2 <- treelist

tauopts_case1_theta2 <- sapply(treelist, gettauopt)


# Pick a good example from case 2, theta1:
treedir <- 'sims/pophetmodel/case2/theta1/'
treenms <- paste0('tree',1:30,'.nwk')
treenms <- paste0(treedir,treenms)

treelist <- list()

for(treenm in treenms){
	treelist[[treenm]] <- read.tree(treenm)
}

treelist_case2_theta1 <- treelist

tauopts_case2_theta1 <- sapply(treelist, gettauopt)

# Pick a good example from case 2, theta2:
treedir <- 'sims/pophetmodel/case2/theta2/'
treenms <- paste0('tree',1:30,'.nwk')
treenms <- paste0(treedir,treenms)

treelist <- list()

for(treenm in treenms){
	treelist[[treenm]] <- read.tree(treenm)
}

treelist_case2_theta2 <- treelist

tauopts_case2_theta2 <- sapply(treelist, gettauopt)

# Choose some representatives:
tree_case1_theta1 <- treelist_case1_theta1[[1]]
tree_case1_theta2 <- treelist_case1_theta2[[1]]
tree_case2_theta1 <- treelist_case2_theta1[[1]]
tree_case2_theta2 <- treelist_case2_theta2[[1]]

p11 <- getmainplot(tree_case1_theta1, taulbi=tauopts_case1_theta1[1],
	title='Case 1: Homogeneous mixing')[[1]]

p12 <- getmainplot(tree_case1_theta2, taulbi=tauopts_case1_theta2[1],
	title='Case 1: Assortative mixing')[[1]]


p21 <- getmainplot(tree_case2_theta1, taulbi=tauopts_case2_theta1[1],
	title='Case 2: Homogeneous mixing')[[1]]

p22 <- getmainplot(tree_case2_theta2, taulbi=tauopts_case2_theta2[1],
	title='Case 2: Assortative mixing')[[1]]


fig_case1 <- plot_grid(p11,p12,nrow=1)
fig_case2 <- plot_grid(p21,p22,nrow=1)

ggsave(fig_case1, file='figures/pophetmodelfigs/fig_case1.png',height=14,width=18)
ggsave(fig_case2, file='figures/pophetmodelfigs/fig_case2.png',height=14,width=18)



