rm(list=ls())

library(ggnewscale)

#compile the model:
source("compile_pophet_and_superspreading_model.R")

# load in the local branching index (LBI) function
source('lbi.R')

#############
#############

ntests <- 500

# time to simulate:
#time <- c(tstart,tmax)/365

time <- c(-125,50)

# when should we start to sample?
tm1 <- 45
tms.lin4 <- runif(ntests,min=tm1,max=time[2])
tms.lin4 <- tms.lin4[order(tms.lin4)]


# what is the total population size?
Npop <- 2.0 * 10^5

# specify the level of preferential mixing (mii, mij), the different
# levels of susceptibility across the groups (zeta1, zeta2), the levels
# of infectiousness across the two grups (theta1, theta2), and the 
# level of superspreading functioning independently of all of this (fH,pH):

gettheta <- function(mii,mij,zeta1,zeta2,theta1,theta2,fH=0.9,pH=0.1){
	
	#of those with active TB, what fraction are superspreaders?
	#pH <- 0.1
	# infectiousness multipliers
	#fH <- 0.9 # fraction of transmission from superspreaders

	# parameterize in terms of equilibrium values of S, L, I
	Seq <- 2/3 * Npop
	Ieq <- 0.01 * Npop
	Leq <- Npop-Seq-Ieq

	# what proportion of newly infected immediately become infectious?
	# 5% develop active TB in the first two years; use this:
	pfast <- 0.05

	# Average duration of untreated pulmonary TB, historically, is ~ 1- 3 years:
	gamma <- 1/1.5

	# The annual risk of developing active TB is 1e-4 - 2e-4 per year:
	sigma <- 1e-4
	 
	# transmission rates (beta=zeta*c*theta):
	beta <- (gamma*Ieq-sigma*Leq)/pfast/Seq/Ieq

	# determine rho:
	rho <- (beta*Seq*Ieq - gamma*Ieq)/Leq

	# susceptibility 
	#zeta1 <- zeta1
	#zeta2 <- zeta2

	# infectiousness
	thetaL <- (1-fH)/(1-pH)
	thetaH <- fH/pH


	# contact rates within/across groups:
	# make c12 and c22 larger to model increased
	# infectiousness from those groups indep of superspreading: 
	c <- beta
	c11 <- beta
	c12 <- beta
	c21 <- beta
	c22 <- beta

	cij <- matrix(c(c11,c12,c21,c22),byrow=T,nrow=2)
	# make the off-diagonals smaller for preferential mixing:
	diag(cij) <- mii * diag(cij)
	cij[1,2] <- mij * cij[1,2]
	cij[2,1] <- mij * cij[2,1]
	#cij <- cij / max(eigen(cij)$values)
	#cij <- cij*c

	c11 <- cij[1,1]	
	c12 <- cij[1,2]	
	c21 <- cij[2,1]	
	c22 <- cij[2,2]	

	# can make c12 and c22 larger to model increased
	# infectiousness from group 2: 
	c11 <- c11 * theta1
	c21 <- c21 * theta1
	c12 <- c12 * theta2
	c22 <- c22 * theta2

	theta <- list(
		beta=beta,
		zeta1 = zeta1,
		zeta2 = zeta2,
		theta1 = theta1,
		theta2 = theta2,
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

# 6 parametrizations: for some of these, we might be able to leave
# the other effect on (superspreading, altsus) but have one 
# uncorrelated with the group status

# The size of the host subpopulations could be important. Let's set it
# to 50% for the "at-risk" group for now and change later if needed
# In all of these, we have dialed the baseline level of superspreading down
# somewhat, to a 20-80 rule. In cases 1 and 3 below, we superimpose some additional
# variation in infectiousness that is correlated with the groups (so that the overall
# level of superpsreading will be more pronounced than a 20-80 rule, but perhaps not
# quite as dramatic as a 10-90 rule).

# 1. No preferential mixing, but altered infectiousness (same as null model)
theta1 <- gettheta(mii=1,mij=1,zeta1=1,zeta2=1,
	theta1=1/2,theta2=2,pH=0.2,fH=0.8)

# 2. No preferential mixing, but altered susceptibility:
theta2 <- gettheta(mii=1,mij=1,zeta1=1/2,zeta2=2,
	theta1=1,theta2=1,pH=0.2,fH=0.8)

# 3. Preferential mixing turned on, and altered infectiousnes correlated w/it:
theta3 <- gettheta(mii=2,mij=.5,zeta1=1,zeta2=1,
	theta1=1/2,theta2=2,pH=0.2,fH=0.8)

# 4. Preferential mixing turned on, and altsus correlated w/it:
theta4 <- gettheta(mii=2,mij=.5,zeta1=1/2,zeta2=2,
	theta1=1,theta2=1,pH=0.2,fH=0.8)


# Set the up the simulation time and the initial conditions:
# set the timestep size:
dT <- 1

# specify initial conditions -- we need to be careful about this now to make
# sure that we have the correct proportion of superspreaders for the
# various cases if we change proportions above:

# 1.
initialStates1 <- c(Npop * 0.5,0,0,0,Npop * 0.5,0,1,0)
names(initialStates1) <- c('S1','L1','IH1','IL1','S2','L2','IH2','IL2') 

# 2.
initialStates2 <- c(Npop * 0.5,0,0,0,Npop * 0.5,0,1,0)
names(initialStates2) <- c('S1','L1','IH1','IL1','S2','L2','IH2','IL2') 

# 3.
initialStates3 <- c(Npop * 0.5,0,0,0,Npop * 0.5,0,1,0)
names(initialStates3) <- c('S1','L1','IH1','IL1','S2','L2','IH2','IL2') 

# 4.
initialStates4 <- c(Npop * 0.5,0,0,0,Npop * 0.5,0,1,0)
names(initialStates4) <- c('S1','L1','IH1','IL1','S2','L2','IH2','IL2') 

out1 <-  sir_simu(
   paramValues = as.list(theta1),
   initialStates = initialStates1,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100)
#   seed=919755)

out2 <-  sir_simu(
   paramValues = as.list(theta2),
   initialStates = initialStates2,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100)
   #seed=280361)

out3 <-  sir_simu(
   paramValues = as.list(theta3),
   initialStates = initialStates3,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100)
   #seed=280361)

out4 <-  sir_simu(
   paramValues = as.list(theta4),
   initialStates = initialStates4,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100)
   #seed=280361)

# use a smaller dataframe for plotting (just sample rows)
traj1 <- out1$traj
if(dim(traj1)[1]<500){plottraj1<-traj1}else{
plottraj1 <- traj1[c(1:1500,sample(1:dim(traj1)[1],1000)),]
plottraj1 <- plottraj1[order(plottraj1$Time),]}
# transform time back to date for plotting:
plottraj1$date <- as.Date(plottraj1$Time * 365)


traj2 <- out2$traj
plottraj2 <- traj2[c(1:1500,sample(1:dim(traj2)[1],1000)),]
plottraj2 <- plottraj2[order(plottraj2$Time),]
# transform time back to date for plotting:
plottraj2$date <- as.Date(plottraj2$Time * 365)

traj3 <- out3$traj
plottraj3 <- traj3[c(1:1500,sample(1:dim(traj3)[1],1000)),]
plottraj3 <- plottraj3[order(plottraj3$Time),]
# transform time back to date for plotting:
plottraj3$date <- as.Date(plottraj3$Time * 365)

traj4 <- out4$traj
plottraj4 <- traj4[c(1:1500,sample(1:dim(traj4)[1],1000)),]
plottraj4 <- plottraj4[order(plottraj4$Time),]
# transform time back to date for plotting:
plottraj4$date <- as.Date(plottraj4$Time * 365)



# plot the trajectories of out1 and out2:

#pdf(file='nullmodelfigs/nullmodeltrajs.pdf',height=5,width=9)
#
#par(mfrow=c(2,2))
#
## low prevalence:
#plot(log10(L1+L2)~date,plottraj1,type='l',col='#005AB5',main='Low prevalence',lwd=2.5,
#	ylab=bquote(Log[10](.('No. of infections'))),xlab='Date',lty='dotted', ylim=c(0,6))
#
#lines(log10(IL1+IL2)~date,plottraj1,type='l',col='#005AB5',main='Active TB',lwd=2.5,
#	ylab='No. of infections',xlab='Date',lty='dashed')
#lines(log10(IH1+IH2)~date,plottraj1,type='l',col='#005AA0',lwd=2.5)
#legend('top',col=c('#005AB5','#005AA0'),
#	lty=c(3,2,1),lwd=2.5,legend=c(bquote(L),bquote(I[1]),bquote(I[2])),
#	cex=1.3)
#
## high prevalence:
#plot(log10(L1+L2)~date,plottraj2,type='l',col='#005AB5',main='High prevalence',lwd=2.5,
#	ylab=bquote(Log[10](.('No. of infections'))),xlab='Date',lty='dotted', ylim=c(0,6))
#
#lines(log10(IL1+IL2)~date,plottraj2,type='l',col='#005AB5',main='Active TB',lwd=2.5,
#	ylab='No. of infections',xlab='Date',lty='dashed')
#lines(log10(IH1+IH2)~date,plottraj2,type='l',col='#005AA0',lwd=2.5)
#
##dev.off()


# use tms.lin4 to generate simulated trees with the same height (approx.) as the
#	empirical lineage 4 tree.
# We want to ensure we are sampling states according to the probabilities they occur in
# the simulation, which will be affected by our parameter choices

# for cross-sectional sampling similar to lineage 4:
tmstart <- 45 #min(tms.lin4)
tmend <- 50 #max(tms.lin4)
ntests <- 2000


x1 <- plottraj1[plottraj1$Time < tmend & plottraj1$Time > tmstart,c('Time','IH1','IL1','IH2','IL2')]
x2 <- plottraj2[plottraj2$Time < tmend & plottraj2$Time > tmstart,c('Time','IH1','IL1','IH2','IL2')]
x3 <- plottraj3[plottraj3$Time < tmend & plottraj3$Time > tmstart,c('Time','IH1','IL1','IH2','IL2')]
x4 <- plottraj4[plottraj4$Time < tmend & plottraj4$Time > tmstart,c('Time','IH1','IL1','IH2','IL2')]

y1 <- colSums(x1)/sum(colSums(x1))
y2 <- colSums(x2)/sum(colSums(x2))
y3 <- colSums(x3)/sum(colSums(x3))
y4 <- colSums(x4)/sum(colSums(x4))

# randomly sample states according to the probabilities with which they were observed in the simulation:
s1 <- c('IH1','IL1','IH2','IL2')[t(rmultinom(ntests,1,prob=y1[c('IH1','IL1','IH2','IL2')])) %*% c(1,2,3,4)]
s2 <- c('IH1','IL1','IH2','IL2')[t(rmultinom(ntests,1,prob=y2[c('IH1','IL1','IH2','IL2')])) %*% c(1,2,3,4)]
s3 <- c('IH1','IL1','IH2','IL2')[t(rmultinom(ntests,1,prob=y3[c('IH1','IL1','IH2','IL2')])) %*% c(1,2,3,4)]
s4 <- c('IH1','IL1','IH2','IL2')[t(rmultinom(ntests,1,prob=y4[c('IH1','IL1','IH2','IL2')])) %*% c(1,2,3,4)]

# create data frames with the tms.lin4 sampling times and the sampled state:
tms.1 <- data.frame(Date=runif(ntests,tmstart,tmend),Comp=s1) 
tms.2 <- data.frame(Date=runif(ntests,tmstart,tmend),Comp=s2) 
tms.3 <- data.frame(Date=runif(ntests,tmstart,tmend),Comp=s3) 
tms.4 <- data.frame(Date=runif(ntests,tmstart,tmend),Comp=s4) 

tms.1 <- tms.1[order(tms.1$Date),]
tms.2 <- tms.2[order(tms.2$Date),]
tms.3 <- tms.3[order(tms.3$Date),]
tms.4 <- tms.4[order(tms.4$Date),]




## longitudinal sampling:
#
## set the minimum time for longitudinal sampling to commence:
#mintime <- -120
#ntests <- 500
#
#tms.1 <- plottraj1[plottraj1$Time > mintime,]
#tms.1 <- tms.1[sort(sample(1:dim(tms.1)[1],ntests)),]
#tms.1$N <- tms.1$IH1+tms.1$IL1+tms.1$IH2+tms.1$IL2
#tms.1$pih1 <- tms.1$IH1/tms.1$N
#tms.1$pil1 <- tms.1$IL1/tms.1$N
#tms.1$pih2 <- tms.1$IH2/tms.1$N
#tms.1$pil2 <- tms.1$IL2/tms.1$N
#
#tms.1$sampled <- apply(tms.1, 1, function(x) c('IH1','IL1','IH2','IL2')[which(t(rmultinom(1,1,prob=x[c('pih1','pil1','pih2','pil2')])) > 0)])
#
#tms.1 <- tms.1[,c('Time','sampled')]
#colnames(tms.1) <- c('Date','Comp')
#
## for Sim 2:
#tms.2 <- plottraj2[plottraj2$Time > mintime,]
#tms.2 <- tms.2[sort(sample(1:dim(tms.2)[1],ntests)),]
#tms.2$N <- tms.2$IH1+tms.2$IL1+tms.2$IH2+tms.2$IL2
#tms.2$pih1 <- tms.2$IH1/tms.2$N
#tms.2$pil1 <- tms.2$IL1/tms.2$N
#tms.2$pih2 <- tms.2$IH2/tms.2$N
#tms.2$pil2 <- tms.2$IL2/tms.2$N
#
#tms.2$sampled <- apply(tms.2, 1, function(x) c('IH1','IL1','IH2','IL2')[which(t(rmultinom(1,1,prob=x[c('pih1','pil1','pih2','pil2')])) > 0)])
#
#tms.2 <- tms.2[,c('Time','sampled')]
#colnames(tms.2) <- c('Date','Comp')
#
## for Sim 3:
#tms.3 <- plottraj3[plottraj3$Time > mintime,]
#tms.3 <- tms.3[sort(sample(1:dim(tms.3)[1],ntests)),]
#tms.3$N <- tms.3$IH1+tms.3$IL1+tms.3$IH2+tms.3$IL2
#tms.3$pih1 <- tms.3$IH1/tms.3$N
#tms.3$pil1 <- tms.3$IL1/tms.3$N
#tms.3$pih2 <- tms.3$IH2/tms.3$N
#tms.3$pil2 <- tms.3$IL2/tms.3$N
#
#tms.3$sampled <- apply(tms.3, 1, function(x) c('IH1','IL1','IH2','IL2')[which(t(rmultinom(1,1,prob=x[c('pih1','pil1','pih2','pil2')])) > 0)])
#
#tms.3 <- tms.3[,c('Time','sampled')]
#colnames(tms.3) <- c('Date','Comp')
#
## for Sim 4:
#tms.4 <- plottraj4[plottraj4$Time > mintime,]
#tms.4 <- tms.4[sort(sample(1:dim(tms.4)[1],ntests)),]
#tms.4$N <- tms.4$IH1+tms.4$IL1+tms.4$IH2+tms.4$IL2
#tms.4$pih1 <- tms.4$IH1/tms.4$N
#tms.4$pil1 <- tms.4$IL1/tms.4$N
#tms.4$pih2 <- tms.4$IH2/tms.4$N
#tms.4$pil2 <- tms.4$IL2/tms.4$N
#
#tms.4$sampled <- apply(tms.4, 1, function(x) c('IH1','IL1','IH2','IL2')[which(t(rmultinom(1,1,prob=x[c('pih1','pil1','pih2','pil2')])) > 0)])
#
#tms.4 <- tms.4[,c('Time','sampled')]
#colnames(tms.4) <- c('Date','Comp')
#



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

# Make sure that pop2frac here matches the initial conditions above
tree1 <- getsimtree(tms.1,out1)
tree2 <- getsimtree(tms.2,out2)
tree3 <- getsimtree(tms.3,out3)
tree4 <- getsimtree(tms.4,out4)


# Make a figure in base R to make sure things are working ok:

##pdf(file='figures/pophetmodelfigs/trees.pdf',width=6,height=9)
#par(mfrow=c(2,2))
#
#tree <- tree1
#plot(tree, type='phylogram',show.tip.label=F,main='Altered infectiousness, no preferential mixing')
#tips_cols <- ifelse(grepl(x=tree$tip.label,pattern = c("1_")),"#005AB5","#DC3220")
#nodes_cols <- ifelse(grepl(x=tree$node.label,pattern="+IH1+"),"#005AB5",
#	ifelse(grepl(x=tree$node.label,pattern="+IL1+"),'#005AB5','#DC3220') )
#tiplabels(pch=20,col=tips_cols)
##nodelabels(pch=20,col=nodes_cols)
#legend('bottomleft',legend=c('Group1','Group2'),col=c('#005AB5','#DC3220'),pch=20)
#
#tree <- tree2
#plot(tree, type='phylogram',show.tip.label=F,main='Altered susceptibility, no preferential mixing')
#tips_cols <- ifelse(grepl(x=tree$tip.label,pattern = c("1_")),"#005AB5","#DC3220")
#nodes_cols <- ifelse(grepl(x=tree$node.label,pattern="+IH1+"),"#005AB5",
#	ifelse(grepl(x=tree$node.label,pattern="+IL1+"),'#005AB5','#DC3220') )
#tiplabels(pch=20,col=tips_cols)
##nodelabels(pch=20,col=nodes_cols)
#
#tree <- tree3
#plot(tree, type='phylogram',show.tip.label=F,main='Altered infectiousness + preferential mixing')
#tips_cols <- ifelse(grepl(x=tree$tip.label,pattern = c("1_")),"#005AB5","#DC3220")
#nodes_cols <- ifelse(grepl(x=tree$node.label,pattern="+IH1+"),"#005AB5",
#	ifelse(grepl(x=tree$node.label,pattern="+IL1+"),'#005AB5','#DC3220') )
#tiplabels(pch=20,col=tips_cols)
##nodelabels(pch=20,col=nodes_cols)
#
#tree <- tree4
#plot(tree, type='phylogram',show.tip.label=F,main='Altered susceptibility + preferential mixing')
#tips_cols <- ifelse(grepl(x=tree$tip.label,pattern = c("1_")),"#005AB5","#DC3220")
#nodes_cols <- ifelse(grepl(x=tree$node.label,pattern="+IH1+"),"#005AB5",
#	ifelse(grepl(x=tree$node.label,pattern="+IL1+"),'#005AB5','#DC3220') )
#tiplabels(pch=20,col=tips_cols)
##nodelabels(pch=20,col=nodes_cols)
#
##dev.off()
#
#
### Want to generate plots of the trees with edges annotated with LBI,
### the variant clade highlighted, and display LBI~time:
#
## function to plot ggtrees colored in with Host subpopulation, and accompanying
## boxplots of LBI~Host subpopulation:
#getfig <- function(tree,title='add a title!',titleadjust=0.45,ptcex=1){
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
#	crud$lbi20 <- lbi(tree, tau=20)
#
#	# add in a column for the state of the node/tip:
#	crud$state <- NA
#	crud[grep('IH1',crud$label[1:ntests]),'state'] <- 'IH1'
#	crud[grep('IL1',crud$label[1:ntests]),'state'] <- 'IL1'
#	crud[grep('IH2',crud$label[1:ntests]),'state'] <- 'IH2'
#	crud[grep('IL2',crud$label[1:ntests]),'state'] <- 'IL2'
#
#	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], function(z) substr(z, regexpr("S+",z)[1]+3, regexpr(".[+]=",z)[1]))
#
#	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms
#
#	# add a column for group:
#	crud$group <- NA
#	crud$group[grep(1,crud$state)] <- 'Group 1'
#	crud$group[grep(2,crud$state)] <- 'Group 2'
#
#
#	# make a ggtree object:
#	p <- ggtree(tree)
#
#	# merge the dataframe with LBI onto it:
#	p <- p %<+% crud
#
#	# plot it!
#	p1 <- p + geom_tippoint(aes(col=group)) + geom_nodepoint(aes(col=group)) + 
#		scale_color_discrete(name='Host subpopulation',
#		type=c('#E1BE6A','#40B0A6'))  +
#		theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face='bold')) +
#		theme(legend.position='none')	
#		#theme(legend.text=element_text(size=18),legend.title=element_text(size=16,face='bold')) +
#		#theme(legend.position='left')
#
#	## make a panel plot
#
#	p1.box <- ggplot(crud) + geom_boxplot(aes(x=group,y=lbi20,fill=group)) +
#		scale_fill_manual(name='Host\nsubpopulation' ,values=c('#E1BE6A','#40B0A6')) + 
#		theme_classic() + 
#		ylab('LBI') + 
#		xlab('Host subpopulation') + 
#		theme(legend.position='none') +
#		theme(axis.text=element_text(size=18), 
#			axis.title=element_text(size=18, face='bold'))
#		
#
#	# make the tree and boxplot figures w/out the title first:
#	alnd <- align_plots(p1,p1.box,align='h',axis='lr')
#	fig.p1 <- plot_grid(alnd[[1]],alnd[[2]],ncol=2, rel_widths=c(1,0.5))
#
#	# create a common title:
#	title <- ggdraw() + 
#	  draw_label(
#	    title,
#	    fontface = 'bold',
#	    x = titleadjust,
#	    hjust = 0.0,
#	    size = 18
#	  ) +
#	  theme(
#	    # add margin on the left of the drawing canvas,
#	    # so title is aligned with left edge of first plot
#	    plot.margin = margin(0, 0, 0, 7)
#	  )
#
#	# make the tree and boxplot figures w/out the title first:
#
#	# add the title:
#	fig.p1 <- plot_grid(title, fig.p1, ncol=1, rel_heights=c(0.1,1))
#
#	return(fig.p1)
#}
#
## create the individual sub-figures:
#fig.1 <- getfig(tree1,'Altered infectiousness',titleadjust=0.40)
#fig.2 <- getfig(tree2,'Altered susceptibility',titleadjust=0.40)
#fig.3 <- getfig(tree3,'Altered infectiousness + preferential mixing',titleadjust=0.30)
#fig.4 <- getfig(tree4,'Altered susceptibility + preferential mixing',titleadjust=0.25)
#
## place them in a combined figure:
#mainfig <- plot_grid(fig.1,fig.2,fig.3,fig.4,byrow=T,nrow=2)
#
#ggsave(mainfig, file='figures/pophetmodelfigs/mainfig_crosssectional.png', dpi=300, width=20,height=10)
#



getmainplot <- function(tree,taulbi=4,tauthd=5,taurels=6,tauclust=6,title='title'){

	# want to add in rows for nodes with times and LBIs
	crud <- data.frame(time = tree$tip.height,
		label = tree$tip.label)

	# the node labels have the times; extract these:
	m<- sapply(tree$node.label, function(z) substr(z, regexpr('=',z)[1]+1, regexpr(',re',z)[1]-1 ) ) 
	crud2 <- data.frame(time=as.numeric(m), 
			label = names(m))
	# the node labels are super clunky, but we need to keep them to match with the tree
	crud <- rbind(crud,crud2)

	# rearrange columns with labels first:
	crud <- crud[,c(2,1)]

	# calculate LBI for the tips and the nodes:
	crud$lbi <- lbi(tree, tau=taulbi)

	# add in a column for the state of the node/tip:
	crud$state <- NA
	crud[grep('IH1',crud$label[1:ntests]),'state'] <- 'Group 1'
	crud[grep('IL1',crud$label[1:ntests]),'state'] <- 'Group 1'
	crud[grep('IH2',crud$label[1:ntests]),'state'] <- 'Group 2'
	crud[grep('IL2',crud$label[1:ntests]),'state'] <- 'Group 2'


	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], function(z) substr(z, regexpr("S+",z)[1]+3, regexpr(".[+]=",z)[1]-0))
	nodenms[grep('IH1',nodenms)] <- 'Group 1'
	nodenms[grep('IL1',nodenms)] <- 'Group 1'
	nodenms[grep('IH2',nodenms)] <- 'Group 2'
	nodenms[grep('IL2',nodenms)] <- 'Group 2'


	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms

	p <- ggtree(tree,layout='rectangular') %<+% crud

	p1 <- p + geom_nodepoint(aes(col=state)) + geom_tree(linewidth=0.60) +
		scale_color_manual(name='Host',
		values=c('Group 1'='#E1BE6A','Group 2'='#40B0A6')) +
		theme(legend.position='none') +
		labs(title=title) + theme(plot.title=element_text(hjust=0.5,face='bold',size=18))
		#theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face='bold')) +
		#theme(legend.text=element_text(size=18),legend.title=element_text(size=16,face='bold')) +
		#theme(legend.position='left')


        # calculate tree height:
        treeheight <- max(node.depth.edgelength(tree))

        # use cophenetic distances to calculate THD and No. of close relatives:
        x = cophenetic(tree)

        # calculate statistics for the tree (THD, LBI, No. of close relatives):

        # calculate THD from cophenetic distances:
        dat <- apply(x, 1, function(z) sum(exp(-z/tauthd)))

        # organize into a dataframe:
        dat <- as.data.frame(dat)
        colnames(dat) <- 'THD'
        dat$label = rownames(dat)
        dat <- dat[,c(2,1)]

        # calculate LBI directly from the tree:
        dat$LBI <- lbi(tree,tau=taulbi)[1:length(tree$tip.label)]

	# save the raw values in columns:
	dat$LBIraw <- dat$LBI

	# standardize statistics for ease of comparison (uncomment to show raw stats):
	dat$LBI <- (dat$LBI-mean(dat$LBI))/sd(dat$LBI)

	# create a ggtree plot:
       	plin4 <- p1
 

	lbidat <- as.data.frame(dat[,'LBIraw'])
	rownames(lbidat) <- rownames(dat)
	colnames(lbidat) <- 'LBI'

        # use a heatmap to visualize the statistics:
        heatfig <-  gheatmap(plin4,lbidat,
                colnames=T, colnames_position="bottom", hjust=0.0,
                colnames_offset_y=-3,colnames_angle=-45,width=0.1)+
                scale_fill_continuous(name='Value of\nLBI\nat tips\n(raw)',
                low='#FEFE62',high='#5D3A9B') +
                theme(plot.margin=unit(c(1,1,3,1),'cm')) +
                coord_cartesian(clip = 'off') +
                ggtitle(title) +
                theme(plot.title=element_text(hjust=0.5,size=18,face="bold")) + 
		theme(axis.text=element_text(size=18), 
			axis.title=element_text(size=18, face='bold')) +
		theme(legend.text = element_text(size=16), legend.key.size = unit(1.0,'cm'),
			legend.title=element_text(size=18)) + 
		guides(color=guide_legend(override.aes=list(linewidth=2)))
	
	#reorder the factor levels in crud$state:
	# If we want to just look at LBI at the tips, we need to just use the first
	# ntests rows of crud:
	crud$state <- factor(crud$state, levels=c('Group 1','Group 2'))

	p1.box <- ggplot(crud[1:ntests,]) + geom_boxplot(aes(x=state,y=lbi,fill=state)) +
		scale_fill_manual(name='Host' ,values=c('Group 1'='#E1BE6A','Group 2'='#40B0A6')) + 
		theme_classic() + 
		ylab('LBI (raw)') + 
		xlab('Subpopulation') + 
		theme(legend.position='none') +
		theme(axis.text=element_text(size=18), 
			axis.title=element_text(size=18, face='bold'))
	

	# make the tree and boxplot figures w/out the title first:
	alnd <- align_plots(heatfig,p1.box,align='v',axis='lr')
	fig.p1 <- plot_grid(alnd[[1]],alnd[[2]],ncol=1, rel_heights=c(2,0.5))



        return(list(fig.p1,dat))
}


# 1. No preferential mixing, but altered infectiousness (same as null model)
fig1 <- getmainplot(tree1,title='Altered infectiousness',taulbi=2)[[1]]

# 2. No preferential mixing, but altered susceptibility:
fig2 <- getmainplot(tree2, title='Altered susceptibility')[[1]]

# 3. Preferential mixing turned on, and altered infectiousnes correlated w/it:
fig3 <- getmainplot(tree3, title='Altered infectiousness +\npreferential mixing',taulbi=2)[[1]]

# 4. Preferential mixing turned on, and altsus correlated w/it:
fig4 <- getmainplot(tree4, title='Altered susceptibility +\npreferential mixing',taulbi=2)[[1]]

fig <- plot_grid(fig1,fig2,fig3,fig4,nrow=1)

ggsave(fig, file='figures/pophetmodelfigs/pophetmodeltrees.png', dpi=600, height=16, width=24)

gettauplot <- function(tree,title='add a title!'){

	# want to add in rows for nodes with times and LBIs
	crud <- data.frame(time = tree$tip.height,
		label = tree$tip.label)

	# the node labels have the times; extract these:
	m<- sapply(tree$node.label, function(z) substr(z, regexpr('=',z)[1]+1, regexpr(',re',z)[1]-1 ) ) 
	crud2 <- data.frame(time=as.numeric(m), 
			label = names(m))
	# the node labels are super clunky, but we need to keep them to match with the tree
	crud <- rbind(crud,crud2)

	# rearrange columns with labels first:
	crud <- crud[,c(2,1)]

	# calculate LBI for the tips and the nodes for lots of tau values:

	# tau vals for plotting:
	tauvals <- seq(-4,2, length=30)
	tauvals <- 10^tauvals

	for(tau in tauvals) crud[,paste0('lbi',tau)] <- lbi(tree, tau=tau)

	# add in a column for the state of the node/tip:
	crud$state <- NA
	crud[grep('IH1',crud$label[1:ntests]),'state'] <- 'Group 1'
	crud[grep('IL1',crud$label[1:ntests]),'state'] <- 'Group 1'
	crud[grep('IH2',crud$label[1:ntests]),'state'] <- 'Group 2'
	crud[grep('IL2',crud$label[1:ntests]),'state'] <- 'Group 2'


	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], function(z) substr(z, regexpr("S+",z)[1]+2, regexpr(".[+]=",z)[1]+0))

	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms

	ratios <- apply(crud[,paste0('lbi',tauvals)], 2, function(x) mean(x[crud$state=='Group 2'])/mean(x[crud$state=='Group 1'])  )

	plot(tauvals, ratios, xlab=bquote(tau),ylab='LBI ratio (Group 2/Group 1)',  
		type='l',lwd=3, ylim=c(0, max(ratios)*1.5),
		cex.axis=1.3,cex.lab=1.3, main=title)
}

pdf(file='figures/pophetmodelfigs/tauplot_altinf.pdf',height=6,width=6)
gettauplot(tree1, title='Altered infectiousness')
dev.off()

pdf(file='figures/pophetmodelfigs/tauplot_altsus.pdf',height=6,width=6)
gettauplot(tree2, title='Altered susceptibility')
dev.off()

pdf(file='figures/pophetmodelfigs/tauplot_altinf_pref.pdf',height=6,width=6)
gettauplot(tree3, title='Altered infectiousness + preferential mixing')
dev.off()

pdf(file='figures/pophetmodelfigs/tauplot_altsus_pref.pdf',height=6,width=6)
gettauplot(tree4, title='Altered susceptibility + preferential mixing')
dev.off()



