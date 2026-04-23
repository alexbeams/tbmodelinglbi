# Calculate the optimal tau for each tree in a given simulation run 


require(ape)

source('lbi.R')


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
	crud[grep('IH',crud$label[1:ntests]),'state'] <- 'IH'
	crud[grep('IL',crud$label[1:ntests]),'state'] <- 'IL'

	# Uncomment this if we want to include node values, rather than just tips:

#	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], 
#		function(z) substr(z, regexpr("S+",z)[1]+3, regexpr(".[+]=",z)[1]+0))
#
#	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms

	# Calculate F statistics from analysis of variance for each tau:
	Fstats <- apply(crud[,paste0('lbi',taulbis)], 2, function(x) summary(aov(x~state,crud))[[1]][['F value']][1])
	tauopt <- taulbis[Fstats==max(Fstats)]
	return(tauopt)
}


# What are the paths for all of the trees we want to load?



####### Case 1, theta1: groups in population equal, homogeneous mixing:
#######

treedir <- 'sims/nullmodel/500/theta1/'
treenms <- paste0('tree',1:30,'.nwk')
treenms <- paste0(treedir,treenms)

# make a big list where we will store all of the matrices to create our ROC plot
# One ROC curve for each configuration of tests (Nyr and ntests)
# ROC curves are averaged over 30 simulations
taulist <- list()

for(treenm in treenms){
	tree <- read.tree(treenm)
	tau <- gettauopt(tree)
	taulist[[treenm]] <- tauopt
}

# average over the simulations:
save(tau_case1_theta1,file='sims/pophetmodel/roc/tau_case1_theta1.Rdata')

####### Case 1, theta2: groups in population equal, preferential mixing:
#######

treedir <- 'sims/pophetmodel/case1/theta2/'
treenms <- paste0('tree',1:30,'.nwk')
treenms <- paste0(treedir,treenms)

# make a big list where we will store all of the matrices to create our ROC plot
# One ROC curve for each configuration of tests (Nyr and ntests)
# ROC curves are averaged over 30 simulations
taulist <- list()

for(treenm in treenms){
	tree <- read.tree(treenm)
	tauopt <- gettauopt(tree)
	taulist[[treenm]] <- tauopt
}

# save the tau values:
save(tau_case1_theta2,file='sims/pophetmodel/roc/tau_case1_theta2.Rdata')



####### Case 2, theta1: groups in population equal, homogeneous mixing:
#######

treedir <- 'sims/pophetmodel/case2/theta1/'
treenms <- paste0('tree',1:30,'.nwk')
treenms <- paste0(treedir,treenms)

# make a big list where we will store all of the matrices to create our ROC plot
# One ROC curve for each configuration of tests (Nyr and ntests)
# ROC curves are averaged over 30 simulations
taulist <- list()

for(treenm in treenms){
	tree <- read.tree(treenm)
	tau <- gettauopt(tree)
	taulist[[treenm]] <- tau	
}

# average over the simulations:
save(tau_case2_theta1,file='sims/pophetmodel/roc/tau_case2_theta1.Rdata')

####### Case 2, theta2: groups in population equal, preferential mixing:
#######

treedir <- 'sims/pophetmodel/case2/theta2/'
treenms <- paste0('tree',1:30,'.nwk')
treenms <- paste0(treedir,treenms)

# make a big list where we will store all of the matrices to create our ROC plot
# One ROC curve for each configuration of tests (Nyr and ntests)
# ROC curves are averaged over 30 simulations
taulist <- list()

for(treenm in treenms){
	tree <- read.tree(treenm)
	tau <- gettauopt(tree)
	taulist[[treenm]] <- tau	
}

# average over the simulations:
tau_case2_theta2 <- apply(simplify2array(taulist),1:2,mean)
save(tau_case2_theta2,file='sims/pophetmodel/roc/tau_case2_theta2.Rdata')


