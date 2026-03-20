rm(list=ls())

# Just dialing up infectiousness independently  - can change later

# This simulates the pophet model under Case 2: The two types of hosts are equally
# represented among infections (and therefore data), but one of the subpopulations
# is smaller than the other, so that it is over-represented due to heightened
# susceptibility relative to the other group. (N1 =\= N2, but I1=I2, L1=L2).

# EDIT THIS: 
# The key parameters that get varied will be 
# 1. f = N2/(N1+N2), the proportion of Group 2 individuals in the population

# 2. w, the parameter describing the level of mixing between the two groups: 
#    0 < w < 1, with w=1 homogeneous mixing and w=0 no cross-group transmission

library(ggnewscale)

#compile the model:
source("compile_pophet_and_superspreading_model.R")

# load in the local branching index (LBI) function
source('lbi.R')

#############
#############

ntests <- 500

# time to simulate:
time <- c(-125,50)

# when should we start to sample?
tm1 <- 45
tms.lin4 <- runif(ntests,min=tm1,max=time[2])
tms.lin4 <- tms.lin4[order(tms.lin4)]


# what is the total population size?
Npop <- 2.0 * 10^5

# Specify the level of preferential mixing, w, and the proportion of 
# Group 2 in the population, f=N2/N. For 0 < f < 0.5, Group 2 is smaller
# but has elevated susceptibility, zeta > 1.
 
gettheta <- function(f,w,fH=0.9,pH=0.1){
	
	N <- Npop

	# parameterize in terms of equilibrium values of S, L, I
	Seq <- 2/3 * Npop
	Ieq <- 0.01 * Npop
	Leq <- Npop-Seq-Ieq

	# what proportion of newly infected immediately become infectious?
	# 5% develop active TB in the first two years; use this:
	pfast <- 0.05
	# This variable corresponds to "p" in our analysis from Maple that gives us zeta, etc.
	p <- pfast

	# Average duration of untreated pulmonary TB, historically, is ~ 1- 3 years:
	gamma <- 1/1.5
	# The solutions fro Maple use this:	
	gam <- gamma

	# The annual risk of developing active TB is 1e-4 - 2e-4 per year:
	sigma <- 1e-4
	
        beta <- (gamma*Ieq-sigma*Leq)/pfast/Seq/Ieq
 
	# determine rho:
        rho <- (beta*Seq*Ieq - gamma*Ieq)/Leq

	# What is the infectiousness of group 2 relative to group 1?
	theta <- 5

	# What is the susceptibility of Group 2, given y=L2/(L1+L2) and a value of w, where w <1
	# describes reduced cross-group contact and w =1 is homogeneous mixing between the groups?
	zeta <- 1 
 
	# Baseline superspreading parameters: infectiousness
	thetaL <- (1-fH)/(1-pH)
	thetaH <- fH/pH

	# specify the transmission rates (but theta into the contact rates b/c the reactions
	#  don't include it separately -- they do include zeta1 and zeta2 already)
	c11 <- beta
	c12 <- w * beta * theta
	c21 <- w * beta
	c22 <- beta * theta
	
	theta <- list(
		zeta1 = 1,
		zeta2 = zeta,
		theta1 = 1,
		theta2 = theta,
		thetaH = thetaH,
		thetaL = thetaL,
		c11 = c11,
		c12 = c12,
		c21 = c21,
		c22 = c22,
		gamma = gamma,
		sigma = sigma,
		rho = rho,
		#pmu = pmu,
		pH = pH,
		pfast=pfast
	)
	return(theta)
}

# create a function to simulate trees given tms and output:
getsimtree <- function(tms,output){
	simulate_tree(
	simuResults=output,
	dates=tms,
	deme=c('IH1','IL1','L1','IH2','IL2','L2'),
	root = 'IH2', # set this to match the initial condition used to produce simuResults=output
	nTrials=50,
	resampling=FALSE,
	addInfos = TRUE)
}


# 1. No preferential mixing, different levels of heightened susceptibility for Group 2:
theta1 <- gettheta(f=0.2,w=1,pH=0.1,fH=0.9)

# 2. No preferential mixing, different levels of heightened susceptibility for Group 2:
theta2 <- gettheta(f=0.2,w=.1,pH=0.1,fH=0.9)


# Make sure these initial conditions match the values we use for Group2 and Group 1:
initialStates <- c(Npop * 0.8,0,0,0,Npop * 0.2,0,1,0)
names(initialStates) <- c('S1','L1','IH1','IL1','S2','L2','IH2','IL2') 

# Folder to save outputs:
case <- 'infectiousness/case2/'

# First, simulate using theta1 (homogeneous mixing):
numouts <- 30
siminds <- 1:numouts

# sub-folder to save outputs:
thetafolder <- 'theta1/'

for(outind in siminds){
 outpath <- paste0('sims/pophetmodel/infectiousness/case2/theta1/',outind,'.txt')

 # Run the sim, save the seed!
 out <-  sir_simu(
    paramValues = as.list(theta1),
    initialStates = initialStates,
    tau = .0001,
    times = time,
    method = "mixed",
    verbose = TRUE,
    nTrials = 100)
 # Save the seed:
 writeLines(as.character(out$seed), paste0('sims/pophetmodel/infectiousness/case2/theta1/seed',outind,'.txt'))

 # Simulate a tree :
 source('tedit_pophet_casex_bigsim_gettrees.R')
} 


# Second, simulate using theta2 (preferential mixing):
numouts <- 30
siminds <- 1:numouts

# sub-folder to save outputs:
thetafolder <- 'theta2/'


for(outind in siminds){
 outpath <- paste0('sims/pophetmodel/infectiousness/case2/theta2/',outind,'.txt')

 # Run the sim, save the seed!
 out <-  sir_simu(
    paramValues = as.list(theta2),
    initialStates = initialStates,
    tau = .0001,
    times = time,
    method = "mixed",
    verbose = TRUE,
    nTrials = 100)
 # Save the seed:
 writeLines(as.character(out$seed), paste0('sims/pophetmodel/infectiousness/case2/theta2/seed',outind,'.txt'))

 # Simulate a tree :
 source('tedit_pophet_casex_bigsim_gettrees.R')
} 



#
## Run the simulations with perfect mixing for different levels of altered suscep
#
#out1 <-  sir_simu(
#   paramValues = as.list(theta1),
#   initialStates = initialStates,
#   tau = .0001,
#   times = time,
#   method = "mixed",
#   verbose = TRUE,
#   nTrials = 100)
#
#out2 <-  sir_simu(
#   paramValues = as.list(theta2),
#   initialStates = initialStates,
#   tau = .0001,
#   times = time,
#   method = "mixed",
#   verbose = TRUE,
#   nTrials = 100)
#
## use a smaller dataframe for plotting (just sample rows)
#traj1 <- out1$traj
#if(dim(traj1)[1]<500){plottraj1<-traj1}else{
#plottraj1 <- traj1[c(1:1500,sample(1:dim(traj1)[1],1000)),]
#plottraj1 <- plottraj1[order(plottraj1$Time),]}
## transform time back to date for plotting:
#plottraj1$date <- as.Date(plottraj1$Time * 365)
#
#traj2 <- out2$traj
#plottraj2 <- traj2[c(1:1500,sample(1:dim(traj2)[1],1000)),]
#plottraj2 <- plottraj2[order(plottraj2$Time),]
## transform time back to date for plotting:
#plottraj2$date <- as.Date(plottraj2$Time * 365)
#
## for cross-sectional sampling similar to lineage 4:
#tmstart <- 45 #min(tms.lin4)
#tmend <- 50 #max(tms.lin4)
#
#
#x1 <- plottraj1[plottraj1$Time < tmend & plottraj1$Time > tmstart,c('Time','IH1','IL1','IH2','IL2')]
#x2 <- plottraj2[plottraj2$Time < tmend & plottraj2$Time > tmstart,c('Time','IH1','IL1','IH2','IL2')]
#
#y1 <- colSums(x1)/sum(colSums(x1))
#y2 <- colSums(x2)/sum(colSums(x2))
#
## randomly sample states according to the probabilities with which they were observed in the simulation:
#s1 <- c('IH1','IL1','IH2','IL2')[t(rmultinom(ntests,1,prob=y1[c('IH1','IL1','IH2','IL2')])) %*% c(1,2,3,4)]
#s2 <- c('IH1','IL1','IH2','IL2')[t(rmultinom(ntests,1,prob=y2[c('IH1','IL1','IH2','IL2')])) %*% c(1,2,3,4)]
#
## create data frames with the tms.lin4 sampling times and the sampled state:
#tms.1 <- data.frame(Date=runif(ntests,tmstart,tmend),Comp=s1) 
#tms.2 <- data.frame(Date=runif(ntests,tmstart,tmend),Comp=s2) 
#
#tms.1 <- tms.1[order(tms.1$Date),]
#tms.2 <- tms.2[order(tms.2$Date),]
#
#
#
#
## create a function to simulate trees given tms and output:
#getsimtree <- function(tms,output){
#	simulate_tree(
#	simuResults=output,
#	dates=tms,
#	deme=c('IH1','IL1','L1','IH2','IL2','L2'),
#	root = 'IH2', # set this to match the initial condition used to produce simuResults=output
#	nTrials=50,
#	resampling=FALSE,
#	addInfos = TRUE)
#}
#
## Make sure that pop2frac here matches the initial conditions above
#tree1 <- getsimtree(tms.1,out1)
#tree2 <- getsimtree(tms.2,out2)
#
#getmainplot <- function(tree,taulbi=4,tauthd=5,taurels=6,tauclust=6,title='title'){
#
#	# want to add in rows for nodes with times and LBIs
#	crud <- data.frame(time = tree$tip.height,
#		label = tree$tip.label)
#
#	# the node labels have the times; extract these:
#	m<- sapply(tree$node.label, function(z) substr(z, regexpr('=',z)[1]+1, regexpr(',re',z)[1]-1 ) ) 
#	crud2 <- data.frame(time=as.numeric(m), 
#			label = names(m))
#	# the node labels are super clunky, but we need to keep them to match with the tree
#	crud <- rbind(crud,crud2)
#
#	# rearrange columns with labels first:
#	crud <- crud[,c(2,1)]
#
#	# calculate LBI for the tips and the nodes:
#	crud$lbi <- lbi(tree, tau=taulbi)
#
#	# add in a column for the state of the node/tip:
#	crud$state <- NA
#	crud[grep('IH1',crud$label[1:ntests]),'state'] <- 'Group 1'
#	crud[grep('IL1',crud$label[1:ntests]),'state'] <- 'Group 1'
#	crud[grep('IH2',crud$label[1:ntests]),'state'] <- 'Group 2'
#	crud[grep('IL2',crud$label[1:ntests]),'state'] <- 'Group 2'
#
#
#	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], function(z) substr(z, regexpr("S+",z)[1]+3, regexpr(".[+]=",z)[1]-0))
#	nodenms[grep('IH1',nodenms)] <- 'Group 1'
#	nodenms[grep('IL1',nodenms)] <- 'Group 1'
#	nodenms[grep('IH2',nodenms)] <- 'Group 2'
#	nodenms[grep('IL2',nodenms)] <- 'Group 2'
#
#
#	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms
#
#	p <- ggtree(tree,layout='rectangular') %<+% crud
#
#	p1 <- p + aes(col=state) + geom_tree(linewidth=0.60) +
#		scale_color_manual(name='Host',
#		values=c('Group 1'='#E1BE6A','Group 2'='#40B0A6')) +
#		theme(legend.position='none') +
#		labs(title=title) + theme(plot.title=element_text(hjust=0.5,face='bold',size=18))
#		#theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face='bold')) +
#		#theme(legend.text=element_text(size=18),legend.title=element_text(size=16,face='bold')) +
#		#theme(legend.position='left')
#
#
#        # calculate tree height:
#        treeheight <- max(node.depth.edgelength(tree))
#
#        # use cophenetic distances to calculate THD and No. of close relatives:
#        x = cophenetic(tree)
#
#        # calculate statistics for the tree (THD, LBI, No. of close relatives):
#
#        # calculate THD from cophenetic distances:
#        dat <- apply(x, 1, function(z) sum(exp(-z/tauthd)))
#
#        # organize into a dataframe:
#        dat <- as.data.frame(dat)
#        colnames(dat) <- 'THD'
#        dat$label = rownames(dat)
#        dat <- dat[,c(2,1)]
#
#        # calculate LBI directly from the tree:
#        dat$LBI <- lbi(tree,tau=taulbi)[1:length(tree$tip.label)]
#
#	# save the raw values in columns:
#	dat$LBIraw <- dat$LBI
#
#	# standardize statistics for ease of comparison (uncomment to show raw stats):
#	dat$LBI <- (dat$LBI-mean(dat$LBI))/sd(dat$LBI)
#
#	# create a ggtree plot:
#       	plin4 <- p1
# 
#
#	lbidat <- as.data.frame(dat[,'LBIraw'])
#	rownames(lbidat) <- rownames(dat)
#	colnames(lbidat) <- 'LBI'
#
#        # use a heatmap to visualize the statistics:
#        heatfig <-  gheatmap(plin4,lbidat,
#                colnames=T, colnames_position="bottom", hjust=0.0,
#                colnames_offset_y=-3,colnames_angle=-45,width=0.1)+
#                scale_fill_continuous(name='Value of\nLBI\nat tips\n(raw)',
#                low='#FEFE62',high='#5D3A9B') +
#                theme(plot.margin=unit(c(1,1,3,1),'cm')) +
#                coord_cartesian(clip = 'off') +
#                ggtitle(title) +
#                theme(plot.title=element_text(hjust=0.5,size=18,face="bold")) + 
#		theme(axis.text=element_text(size=18), 
#			axis.title=element_text(size=18, face='bold')) +
#		theme(legend.text = element_text(size=16), legend.key.size = unit(1.0,'cm'),
#			legend.title=element_text(size=18)) + 
#		guides(color=guide_legend(override.aes=list(linewidth=2)))
#	
#	#reorder the factor levels in crud$state:
#	# If we want to just look at LBI at the tips, we need to just use the first
#	# ntests rows of crud:
#	crud$state <- factor(crud$state, levels=c('Group 1','Group 2'))
#
#	p1.box <- ggplot(crud[1:ntests,]) + geom_boxplot(aes(x=state,y=lbi,fill=state)) +
#		scale_fill_manual(name='Host' ,values=c('Group 1'='#E1BE6A','Group 2'='#40B0A6')) + 
#		theme_classic() + 
#		ylab('LBI (raw)') + 
#		xlab('Subpopulation') + 
#		theme(legend.position='none') +
#		theme(axis.text=element_text(size=18), 
#			axis.title=element_text(size=18, face='bold'))
#	
#
#	# make the tree and boxplot figures w/out the title first:
#	alnd <- align_plots(heatfig,p1.box,align='v',axis='lr')
#	fig.p1 <- plot_grid(alnd[[1]],alnd[[2]],ncol=1, rel_heights=c(2,0.5))
#
#
#
#        return(list(fig.p1,dat))
#}
#
#
